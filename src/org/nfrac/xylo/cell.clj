(ns org.nfrac.xylo.cell
  (:require [org.nfrac.xylo.dna :as dna]
            [org.nfrac.str-alignment.core :as ali]
            [clojure.spec.alpha :as s]
            [clojure.test.check.random :as random]))

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
  (new-cell "x1y2z3 qtrusv abcd=.1234?="))

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
  "Returns a collection of maps {:offset :score :vs-offset}, where
  offset is an index into open-dna and the score is from
  Smith-Waterman matching."
  [open-dna vs-dna]
  (let [amat (ali/alignments open-dna vs-dna {:match-weight 2})
        matches (ali/distinct-local-matches amat 4)]
    (->> matches
         (map (fn [locs]
                (let [loc (last locs)
                      score (:score (get amat loc))]
                  {:offset (first loc)
                   :vs-offset (second loc)
                   :score score}))))))

(defn all-binding-sites
  "Returns a collection of maps {:offset :score :path}, where offset is
  an index into open-dna. :path describes the source of each match
  as [kind i j]; kind is :stimuli or :products, j is the index into
  its dna-like sequence. Ordered by score decreasing."
  [open-dna products stimuli]
  (->>
   (concat (map vector (repeat :products) (range) products)
           (map vector (repeat :stimuli) (range) stimuli))
   (reduce (fn [out [kind kind-i vs-dna]]
             (->> (binding-sites open-dna vs-dna)
                  (map (fn [m]
                         (assoc m :path [kind kind-i (:vs-offset m)])))
                  (into out)))
           ())
   (sort-by :score >)))

(defn select-binding-site
  ""
  [cell stimuli time-step]
  (let [odna (open-dna (::dna cell) (::dna-open? cell))
        products (::products cell)
        binds (all-binding-sites odna products stimuli)
        cum (reductions + (map :score binds))
        sum (last cum)
        t' (mod time-step sum)
        bind-index (count (take-while #(< % t') cum))]
    (nth binds bind-index)))

(defmulti reaction-op
  "Reaction operations can change:

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
  * :offset -- the next offset to read, or :stop-reaction; nil means increment.
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

(defmethod reaction-op 'product
  [op cell open-dna offset]
  (let [cost 1
        energy (:energy cell)]
    (if (>= energy cost)
      {:delayed {:product nil}
       :cell-immediate {:energy (- energy cost)}}
      ;; not enough energy
      {:offset :stop-reaction})))

(defmethod reaction-op 'silence
  [op cell open-dna offset]
  {:offset :terminate})

(defmethod reaction-op 'to-sun
  [op cell open-dna offset]
  {:cell-immediate {:orientation 0.0}})

(defmethod reaction-op 'about-face
  [op cell open-dna offset]
  (let [angle (:orientation cell)]
    {:cell-immediate {:orientation (+ angle Math/PI)}}))


(defn react-at-site
  [cell open-dna bind-offset]
  )
