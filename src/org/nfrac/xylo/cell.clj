(ns org.nfrac.xylo.cell
  (:require [org.nfrac.xylo.dna :as dna]
            [org.nfrac.str-alignment.core :as ali]
            [clojure.spec.alpha :as s]
            [clojure.test.check.random :as random]))

(def match-weight 2)
(def min-binding-score (* match-weight (+ 1 dna/codon-length)))
(def min-template-bases (* 2 dna/codon-length))
(def max-ops 36)
(def max-energy 8)

(s/def ::dna
  (s/coll-of ::dna/base, :min-count 1))

(s/def ::dna-open?
  (s/coll-of boolean?, :min-count 1, :kind vector?))

(s/def ::energy nat-int?)

(s/def ::products
  (s/map-of ::dna nat-int?))

(s/def ::cell
  (s/keys :req [::dna
                ::dna-open?
                ::products
                ::energy]))

(defn new-cell
  [dna]
  (let []
    {::dna dna
     ::dna-open? (vec (repeat (count dna) true))
     ::products {}
     ::energy 1}))

(defn seed-cell
  []
  (->>
   (concat (first dna/terminators)
           (:sun dna/fixed-stimuli)
           (dna/op->codon 'to-stimulus)
           (dna/op->codon 'clone)
           (dna/op->codon 'bond-form)
           (dna/op->codon 'stop-reaction)
           (first dna/terminators)
           (:ground dna/fixed-stimuli)
           (first dna/terminators)
           (last dna/no-ops)
           (last dna/no-ops)
           (dna/complement (last dna/no-ops))
           (dna/complement (last dna/no-ops))
           (dna/op->codon 'energy-test)
           (first dna/no-ops)
           (first (drop 2 dna/no-ops))
           (first dna/terminators)
           (dna/op->codon 'stop-reaction)
           (first dna/no-ops)
           (first (drop 2 dna/no-ops))
           (first dna/terminators)
           (dna/op->codon 'sugar-start)
           (dna/op->codon 'stop-reaction)
           )
   (apply str)
   (new-cell)))

(defn vector-subset
  [xs ?s]
  (loop [xs xs
         ?s ?s
         out (transient [])]
    (if-let [x (first xs)]
      (recur (rest xs)
             (rest ?s)
             (if (first ?s)
               (conj! out x)
               out))
      ;; done
      (persistent! out))))

(defn open-dna
  [dna dna-open?]
  (apply str (vector-subset dna dna-open?)))

(defn binding-sites
  "Returns a collection of possible binding sites of vs-dna onto open-dna.
  Only sites with a sufficient score from Smith-Waterman matching are
  included. Keys bind-begin and bind-end-x give the start (inclusive)
  and end (exclusive) indexes to open-dna of the match, aligned to
  codon boundaries."
  [open-dna vs-dna min-score]
  (let [amat (ali/alignments open-dna vs-dna {:match-weight match-weight})
        matches (ali/distinct-local-matches amat min-score)
        nnn dna/codon-length]
    (->> matches
         (map (fn [locs]
                (let [loc (last locs)
                      score (:score (get amat loc))
                      begin-loc (first locs)]
                  {:bind-begin-base (first loc)
                   :vs-end-base (second loc)
                   :bind-begin (let [x (first begin-loc)] (- x (mod x nnn)))
                   :bind-end-x (let [x (first loc)] (+ nnn (- x (mod x nnn))))
                   :score score}))))))

(defn all-binding-sites
  "Returns a collection of maps {:offset :score :path}, where offset is
  an index into open-dna. :path describes the source of each match
  as [kind i j]; kind is :stimuli or :products, j is the index into
  its dna-like sequence. Ordered by score decreasing."
  [open-dna products stimuli min-score]
  (->>
   (concat (map vector (repeat :products) (range) products)
           (map vector (repeat :stimuli) (range) stimuli))
   (reduce (fn [out [kind kind-i vs-dna]]
             (->> (binding-sites open-dna vs-dna min-score)
                  (map (fn [m]
                         (assoc m :path [kind kind-i (:vs-end-base m)])))
                  (into out)))
           ())
   (sort-by :score >)))

(defn select-binding-site
  ""
  [cell stimuli time-step]
  (let [odna (open-dna (::dna cell) (::dna-open? cell))
        products (::products cell)
        binds (all-binding-sites odna products stimuli min-binding-score)
        cum (reductions + (map :score binds))
        sum (last cum)
        t' (mod time-step sum)
        bind-index (count (take-while #(< % t') cum))]
    (nth binds bind-index)))

(defmulti reaction-op
  "Args are
  * op - the symbol or keyword interpreted from current codon.
  * cell - current cell state.
  * open-dna - dna sequence omitting any silenced parts.
  * offset - current offset _after_ current op codon.

  Reaction operations can change:

  _NOW_
  * current offset into open DNA (read head)
  * stop current reaction
  * cell energy
  * cell orientation

  _LATER_
  * cell DNA silencing
  * cell products
  * world - creating new cells - with initialisation for specialisation
  * world - forming bonds
  * world - force
  * world - sugar channels

  Therefore, this function returns keys
  * :offset -- the next offset to read, or :stop-reaction; nil means no change.
  * :cell-immediate -- to be merged into cell
    - :energy
    - :orientation
    - :dna-open?
  * :delayed
    - :product
    - :new-cell
      - :in-direction
      - :init-offset
    - :new-bond
      - :in-direction
    - :break-bond
      - :in-direction
    - :new-channel
    - :break-channel
  "
  (fn [op cell open-dna offset]
    op))

(defmethod reaction-op :no-op
  [op cell open-dna offset]
  {})

(defmethod reaction-op :terminator
  [op cell open-dna offset]
  {})

(defmethod reaction-op 'stop-reaction
  [op cell open-dna offset]
  {:offset :stop-reaction})

(defn read-template
  [open-dna offset]
  (loop [dna (drop offset open-dna)
         tem []]
    (if (seq dna)
      (let [codon (vector (take dna/codon-length dna))
            op (dna/codon->op codon)]
        (case op
          :terminator
          tem
          'stop-reaction
          tem
          ;; otherwise, append and continue
          (recur (drop dna/codon-length dna) (into tem codon))))
      ;; end of DNA
      tem)))

(defmethod reaction-op 'product
  [op cell open-dna offset]
  (let [cost 1
        energy (:energy cell)]
    (if (>= energy cost)
      (let [tem (read-template open-dna offset)
            terminator dna/codon-length]
        (if (>= (count tem) min-template-bases)
          {:delayed {:product (map dna/translate tem)}
           :offset (+ offset (count tem) terminator)
           :cell-immediate {:energy (- energy cost)}}
          ;; template too short
          {:offset (+ offset (count tem) terminator)}))
      ;; not enough energy
      {:offset :stop-reaction})))

(defmethod reaction-op 'silence
  [op cell open-dna offset]
  {})

(defmethod reaction-op 'unsilence
  [op cell open-dna offset]
  {})

(defmethod reaction-op 'goto
  [op cell open-dna offset]
  (let [tem (read-template open-dna offset)
        terminator dna/codon-length]
    (if (>= (count tem) min-template-bases)
      ;; find best matching site
      (let [matches (binding-sites open-dna tem min-binding-score)]
        (if (seq matches)
          (let [match (apply max-key :score matches)]
            {:offset (:bind-end-x match)})
          ;; no match
          {:offset (+ offset (count tem) terminator)}))
      ;; template too short
      {:offset (+ offset (count tem) terminator)})))

(defmethod reaction-op 'energy-test
  [op cell open-dna offset]
  (let [e-threshold (* max-energy (/ offset (count open-dna)))]
    (if (> (:energy cell) e-threshold)
      (reaction-op 'goto cell open-dna offset)
      ;; test failed, skip goto
      (let [tem (read-template open-dna offset)
            terminator dna/codon-length]
        {:offset (+ offset (count tem) terminator)}))))

(defmethod reaction-op 'to-sun
  [op cell open-dna offset]
  {:cell-immediate {:orientation 0.0}})

(defmethod reaction-op 'about-face
  [op cell open-dna offset]
  (let [angle (:orientation cell)]
    {:cell-immediate {:orientation (+ angle Math/PI)}}))


(defn react-at-site
  [cell open-dna bind-offset]
  (loop [offset bind-offset
         cell cell
         delayed []
         counter 0]
    (cond
      (> counter max-ops)
      {:cell cell :delayed delayed}
      (>= offset (count open-dna))
      {:cell cell :delayed delayed}
      :else
      (let [codon (vector (take dna/codon-length (drop offset open-dna)))
            _ (assert (= dna/codon-length (count codon)))
            op (dna/codon->op codon)
            next-off* (+ offset dna/codon-length)
            re (reaction-op op cell open-dna next-off*)
            next-off (or (:offset re) next-off*)
            ]
        (if (or (= next-off :stop-reaction)
                (>= next-off (count open-dna)))
          {:cell cell :delayed delayed}
          (recur next-off
                 (merge cell (:cell-immediate re))
                 (conj delayed (:delayed re))
                 (inc counter)))))))

;; TODO instrument fn specs
;; invariant: dna / open-dna should be a multiple of codon-length
