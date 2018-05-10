(ns org.nfrac.xylo.dna
  (:require [org.nfrac.xylo.align :as ali]
            [clojure.test.check.random :as random]
            [clojure.spec.gen.alpha :as gen]
            [clojure.spec.alpha :as s])
  (:refer-clojure :exclude [complement bases rand-int rand-nth]))

(def align-options
  {:match-weight 1
   :mismatch-weight -1.5
   :gap-open-weight 4
   :gap-ext-weight 1.5})

(def bases
  (seq "acgt"))

(def base? (set bases))

(s/def ::base base?)

(def codon-length 3)

(def n-codons (int (Math/pow (count bases) codon-length)))

(def codons
  (for [a bases
        b bases
        c bases]
    [a b c]))

(def complement (zipmap bases (drop 2 (cycle bases))))
(def translate (zipmap bases (drop 1 (cycle bases))))
(def untranslate (zipmap (vals translate) (keys translate)))

(defn codon-boundary? [n] (zero? (mod n codon-length)))

(s/def ::dna
  (->
   (s/and string?
          #(every? base? %)
          #(codon-boundary? (count %)))
   (s/with-gen
     (fn []
       (->> (s/gen (s/every ::base, :min-count (* 3 codon-length)))
            (gen/fmap (fn [es]
                        (let [extra (mod (count es) codon-length)]
                          (apply str (drop extra es))))))))))

(s/def ::dna-open?
  (s/every boolean?, :min-count 1, :kind vector?))

(def op-names
  '[to-sun
    about-face
    rot-left
    rot-right
    clone
    sex
    push
    bond-form
    bond-break
    sugar-start
    sugar-stop
    ;; non-physical - control flow
    product
    silence
    unsilence
    goto
    energy-test
    stop-reaction])

(def n-ops (count op-names))
(def n-silence-terminators 6)
(def n-terminators 12)

;; excludes terminators and no-ops (as ambiguous)
(def op->codon
  (zipmap op-names codons))

(def codon->op
  (zipmap codons (concat op-names
                         (repeat n-terminators :terminator)
                         (repeat n-silence-terminators :silence-terminator)
                         (repeat :no-op))))

(s/def ::op-code (set (vals codon->op)))

(def terminators (keep (fn [[codon op]] (when (= op :terminator) codon)) codon->op))
(def no-ops (keep (fn [[codon op]] (when (= op :no-op) codon)) codon->op))

(def fixed-stimuli
  (zipmap [:sun :ground :birth-default]
          (map #(apply str (apply concat %))
               (partition 2 2 no-ops))))

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
  "Cuts out any silenced parts of the DNA."
  [dna dna-open?]
  (apply str (vector-subset dna dna-open?)))

(s/fdef open-dna
        :args (s/cat :dna ::dna :dna-open? ::dna-open?)
        :ret ::dna)

(defn offset-into-full-dna
  [offset-into-open-dna dna-open?]
  (->> dna-open?
       (map #(if % 1 0))
       (reductions + 0)
       (take-while #(<= % offset-into-open-dna))
       (count)
       (dec)))

(defn offset-into-open-dna
  "Truncates to most recent open base (towards start of dna). Returns -1
  if the offset is before any open dna."
  [offset-into-full-dna dna-open?]
  (->> (take (inc offset-into-full-dna) dna-open?)
       (filter true?)
       (count)
       (dec)))

(s/def ::rng (-> #(satisfies? random/IRandom %)
                 (s/with-gen #(gen/fmap random/make-random (gen/int)))))

(defn rand-int
  "Uniform integer between lower (inclusive) and upper (exclusive)."
  [rng lower upper]
  (-> (random/rand-double rng)
      (* (- upper lower))
      (+ lower)
      (Math/floor)
      (long)))

(defn rand-nth
  [rng xs]
  (nth xs (rand-int rng 0 (count xs))))

(def mutation-point-prob (double (/ 1 100)))
(def mutation-delete-prob (double (/ 1 20)))
(def mutation-shift-prob (double (/ 1 20)))
(def mutation-dup-prob (double (/ 1 20)))
(def crossover-freq-bases (* 50 codon-length))

(defn mutate-pointwise
  [[dna dna-open?] rng]
  (loop [i 0
         ds (transient (vec dna))
         rng rng]
    (if (< i (count ds))
      (let [[rng rng*] (random/split rng)]
        (if (< (random/rand-double rng*) mutation-point-prob)
          (let [[rng rng*] (random/split rng)
                x (rand-nth rng* bases)]
            (recur (inc i) (assoc! ds i x) rng))
          (recur (inc i) ds rng)))
      ;; done
      [(apply str (persistent! ds)) dna-open?])))

(defn mutate-delete
  [[dna dna-open?] rng]
  (let [[r1 r2] (random/split-n rng 2)
        i (rand-int r1 0 (- (count dna) codon-length))
        j* (rand-int r2 0 (- (count dna) i))
        j (-> j* (quot codon-length) (* codon-length))]
    [(->> (concat (take i dna)
                  (drop j dna))
          (apply str))
     (->> (concat (take i dna-open?)
                  (drop j dna-open?))
          (vec))]
    ))

(defn coll-mutate-shift
  "Delete a section between i and j, and insert it at k."
  [coll i j k]
  (let [coll* (concat (take i coll)
                      (drop j coll))]
    (concat (take k coll*)
            (take (- j i) (drop i coll))
            (drop k coll*))))

(defn mutate-shift
  [[dna dna-open?] rng]
  (let [[r1 r2 r3] (random/split-n rng 3)
        ;; delete a section between i and j, and insert it at k
        i (rand-int r1 0 (- (count dna) codon-length))
        j* (rand-int r2 0 (- (count dna) i))
        j (-> j* (quot codon-length) (* codon-length))
        k (rand-int r3 0 (count dna))]
    [(->> (coll-mutate-shift dna i j k)
          (apply str))
     (->> (coll-mutate-shift dna-open? i j k)
          (vec))]
    ))

(defn coll-mutate-dup
  "At k, insert a copy of the section between i and j."
  [coll i j k]
  (concat (take k coll)
          (take (- j i) (drop i coll))
          (drop k coll)))

(defn mutate-dup
  [[dna dna-open?] rng]
  (let [[r1 r2 r3] (random/split-n rng 3)
        ;; take a section between i and j and insert a copy of it to k
        i (rand-int r1 0 (- (count dna) codon-length))
        j* (rand-int r2 0 (- (count dna) i))
        j (-> j* (quot codon-length) (* codon-length))
        k (rand-int r3 0 (count dna))]
    [(->> (coll-mutate-dup dna i j k)
          (apply str))
     (->> (coll-mutate-dup dna-open? i j k)
          (vec))]
    ))

(defn mutate
  [dna dna-open? rng]
  (let [[r1 r2 r3 r4 r5 r6 r7] (random/split-n rng 7)]
    (cond-> [dna dna-open?]
      true
      (mutate-pointwise r1)
      (< (random/rand-double r2) mutation-dup-prob)
      (mutate-dup r3)
      (< (random/rand-double r4) mutation-shift-prob)
      (mutate-shift r5)
      (< (random/rand-double r6) mutation-delete-prob)
      (mutate-delete r7))))

(s/fdef mutate
        :args (s/cat :dna ::dna
                     :dna-open? ::dna-open?
                     :rng ::rng)
        :ret (s/tuple ::dna ::dna-open?))

(defn crossover
  "Finds the best global alignment between the two DNA sequences,
  chooses a number of crossover points, and alternates the sequences
  between those points."
  [dna1 dna2 rng]
  (let [match (ali/align dna1 dna2 (assoc align-options
                                          :global? true))
        mpath (:match-path match)
        n-cross (inc (quot (max (count dna1) (count dna2))
                           crossover-freq-bases))
        rngs (random/split-n rng n-cross)
        xlocs (sort (for [r rngs] (rand-nth r mpath)))]
    (loop [xlocs xlocs
           [prev1 prev2] [0 0]
           dna []
           use-2? false]
      (if-let [[i1 i2] (first xlocs)]
        (let [n1 (- i1 prev1)
              n2 (- i2 prev2)]
          (recur (rest xlocs)
                 [i1 i2]
                 (concat dna (if use-2?
                               (take n2 (drop prev2 dna2))
                               (take n1 (drop prev1 dna1))))
                 (not use-2?)))
        ;; done
        (let [dna (concat dna (if use-2?
                                (drop prev2 dna2)
                                (drop prev1 dna1)))]
          (->> dna
               (take (-> (count dna) (quot codon-length) (* codon-length)))
               (apply str)))))))

(s/fdef crossover
        :args (s/cat :dna1 ::dna
                     :dna2 ::dna
                     :rng ::rng)
        :ret ::dna)
