(ns org.nfrac.xylo.align
  (:require [clojure.spec.alpha :as s])
  (:import (parasail Aligner
                     Matrix
                     Cigar
                     Result)))

(def cigar-op
  {(int \=) :match
   (int \X) :mismatch
   (int \D) :deletion
   (int \I) :insertion})

(s/def ::match-weight pos?)
(s/def ::mismatch-weight neg?)
(s/def ::gap-open-weight pos?)
(s/def ::gap-ext-weight pos?)
(s/def ::alphabet string?)

(s/def ::opts
  (s/keys :opt-un [::match-weight
                   ::mismatch-weight
                   ::gap-open-weight
                   ::gap-ext-weight
                   ::alphabet
                   ::w-matrix]))

(s/def ::match
  (s/keys :req-un [::match-path
                   ::score]))

(defn cigar->path
  [cigar query-start ref-start global?]
  (let [q0 query-start
        r0 ref-start]
    (loop [cigar cigar
           prev-q (dec q0)
           prev-r (dec r0)
           ans []
           first? true]
      (if-let [[n op] (first cigar)]
        (let [qrs (case op
                    (:match :mismatch)
                    (for [i (range n)] [(+ prev-q i 1) (+ prev-r i 1)])
                    :deletion
                    (for [i (range n)] [prev-q (+ prev-r i 1)])
                    :insertion
                    (for [i (range n)] [(+ prev-q i 1) prev-r]))
              ]
          (recur (rest cigar)
                 (int (+ prev-q (if (= op :deletion) 0 n)))
                 (int (+ prev-r (if (= op :insertion) 0 n)))
                 (if (and first? (not= op :match)
                          (not global?))
                   ans ;; discard spurious leading mismatch
                   (into ans qrs))
                 false))
        ;; done
        ans))))

(defn cigar->clj
  [^Cigar cigar]
  (mapv (fn [cop] [(Cigar/decode_len cop)
                   (cigar-op (Cigar/decode_op cop))])
        (take (.getLen cigar) (.getSeq cigar))))

(defn align
  "Finds the best match between the strings with either local (default)
  or global alignment.

  Remember that either a substitution matrix must be supplied or the
  alphabet specified."
  [s1 s2 {:as opts
          :keys [match-weight
                 mismatch-weight
                 gap-open-weight
                 gap-ext-weight
                 global?
                 w-matrix
                 alphabet]
          :or {match-weight 2
               mismatch-weight -1
               gap-open-weight 5
               gap-ext-weight 1
               global? false
               alphabet "acgt"}}]
  (let [w-matrix (or w-matrix
                     (Matrix/create alphabet match-weight mismatch-weight))
        result (if global?
                 (Aligner/nw_trace_striped_sat s1 s2 gap-open-weight gap-ext-weight w-matrix)
                 (Aligner/sw_trace_striped_sat s1 s2 gap-open-weight gap-ext-weight w-matrix))
        cigar-raw (.getCigar result s1 s2 w-matrix)
        cigar (cigar->clj cigar-raw)
        ans {:score (.getScore result)
             :end-query (.getEndQuery result)
             :end-ref (.getEndRef result)
             :match-path (cigar->path cigar (.getBegQuery cigar-raw)
                                      (.getBegRef cigar-raw) global?)
             :cigar cigar
             }]
    (.delete cigar-raw)
    (.delete result)
    ans))

(s/fdef align
        :args (s/cat :s1 string?
                     :s2 string?
                     :opts ::opts)
        :ret ::match)

(defn mask-out
  "Replace substring from from (inclusive) to to (exclusive) with char."
  [s from to char]
  (let [s (str s)]
    (str (.substring s 0 from)
         (apply str (repeat (- to from) char))
         (.substring s to (.length s)))))

(defn multi-align
  "Iteratively finds the best matches between the strings, while masking
  out previously matches sections of s1. Stops when the best remaining
  match score is below min-score. Return is ordered by score decreasing.

  Remember that either a substitution matrix must be supplied or the
  alphabet specified."
  [s1 s2 min-score max-sites opts]
  (loop [s1-mask s1
         ans []]
    (let [a (align s1-mask s2 opts)]
      (if (>= (:score a) min-score)
        (let [[q0 _] (first (:match-path a))
              q1 (:end-query a)
              new-ans (conj ans a)]
          (if (= (count new-ans) max-sites)
            new-ans
            (recur (mask-out s1-mask q0 q1 \_)
                   new-ans)))
        ;; done
        ans))))

(s/fdef multi-align
        :args (s/cat :s1 string?
                     :s2 string?
                     :min-score pos?
                     :opts ::opts)
        :ret (s/every ::match, :kind vector?))
