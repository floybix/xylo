(ns org.nfrac.xylo.cell
  (:require [org.nfrac.xylo.dna :as dna]
            [org.nfrac.xylo.physics :as phys]
            [org.nfrac.xylo.physics-grid :refer [in-pi-pi]]
            [org.nfrac.xylo.align :as ali]
            [clojure.spec.alpha :as s]
            [clojure.test.check.random :as random]))

(def max-binding-sites 4)
(def baseline-score 5)
(def weight-power 1)
(def min-template-bases (* 2 dna/codon-length))
(def max-ops 36)
;; note cell metabolism requirement is fixed at 1 unit per time step
(def energies
  {:max 10
   :sun 3
   :sugar 1
   :bond 1
   :clone 4
   :product 1
   :push 2})
(def starvation-steps 8)

(s/def ::product-counts
  (s/map-of ::dna/dna pos-int?))

(s/def ::energy (s/int-in 0 (inc (:max energies))))

(s/def ::orientation (s/double-in :min (- -0.01 Math/PI)
                                  :max (+ 0.01 Math/PI) :NaN? false))

(s/def ::starvation (s/int-in 0 (inc starvation-steps)))

(s/def ::birthstep nat-int?)

(s/def ::cell
  (->
   (s/keys :req-un [::dna/dna
                    ::dna/dna-open?
                    ::product-counts
                    ::orientation
                    ::energy
                    ::starvation
                    ::birthstep])
   (s/and #(= (count (:dna %)) (count (:dna-open? %))))))

(defn new-cell
  [dna time-step]
  (->
    {:dna dna
     :dna-open? (vec (repeat (count dna) true))
     :product-counts {}
     :orientation 0.0
     :energy 1
     :starvation 0
     :birthstep time-step}
    (vary-meta merge
               {::open-dna-cache (atom nil)})))

(s/fdef new-cell
        :args (s/cat :dna ::dna/dna
                     :time-step nat-int?)
        :ret ::cell)

(defn example-dna
  []
  (->>
   (concat (first dna/terminators)
           (dna/fixed-stimuli :sun)
           (dna/op->codon 'rot-left)
           (dna/op->codon 'clone)
           (dna/op->codon 'bond-form)
           (dna/op->codon 'stop-reaction)
           (first dna/terminators)
           (dna/fixed-stimuli :ground)
           (first dna/terminators)
           (dna/op->codon 'stop-reaction)
           (last dna/no-ops)
           (second dna/no-ops)
           (map dna/complement (last dna/no-ops))
           (map dna/complement (second dna/no-ops))
           (dna/op->codon 'energy-test)
           (first dna/no-ops)
           (first (drop 2 dna/no-ops))
           (first dna/terminators)
           (dna/op->codon 'sugar-stop)
           (dna/op->codon 'stop-reaction)
           (first dna/no-ops)
           (first (drop 2 dna/no-ops))
           (first dna/terminators)
           (dna/op->codon 'sugar-start)
           (dna/op->codon 'stop-reaction)
           (first dna/sterminators)
           (dna/op->codon 'rot-right)
           (dna/op->codon 'rot-right)
           (dna/op->codon 'sex)
           (dna/op->codon 'product)
           (first dna/no-ops)
           (second dna/sterminators)
           (second dna/no-ops)
           (last dna/terminators)
           (dna/op->codon 'silence)
           (dna/op->codon 'sugar-start)
           (dna/op->codon 'stop-reaction)
           (last dna/sterminators)
           (dna/op->codon 'unsilence)
           (dna/op->codon 'energy-test)
           (dna/op->codon 'bond-form)
           (dna/op->codon 'energy-test)
           (last dna/sterminators)
           (dna/op->codon 'rot-left)
           (dna/op->codon 'sex)
           (dna/op->codon 'about-face)
           (dna/op->codon 'push)
           (first dna/terminators)
           )
   (apply str)))

(defn get-open-dna
  [cell]
  (let [od-cache (::open-dna-cache (meta cell))
        dna (:dna cell)
        dna-open? (:dna-open? cell)
        hash-key [(hash dna) (hash dna-open?)]]
    (if (= hash-key (:hash-key @od-cache))
      (:val @od-cache)
      (:val (reset! od-cache {:hash-key hash-key
                              :val (dna/open-dna dna dna-open?)})))))

(s/fdef get-open-dna
        :args (s/cat :cell ::cell)
        :ret ::dna/dna)

(s/def ::bind-begin (s/and nat-int? dna/codon-boundary?))
(s/def ::bind-end-x (s/and nat-int? dna/codon-boundary?))
(s/def ::score (s/double-in :min 0 :NaN? false))

(s/def ::bind
  (s/keys :req-un [::bind-begin
                   ::bind-end-x
                   ::score]))

(s/def ::path (s/tuple #{:products :stimuli} nat-int? nat-int?))

(s/def ::bind-with-path
  (s/and ::bind
         (s/keys :req-un [::path])))

(defn binding-sites
  "Finds possible binding sites of vs-dna onto dna, with scores.
  Only sites with a sufficient score from Smith-Waterman matching are
  included. Keys bind-begin and bind-end-x give the start (inclusive)
  and end (exclusive) indexes to dna of the match, aligned to codon
  boundaries."
  [dna vs-dna min-score max-sites]
  (let [matches (ali/multi-align dna vs-dna min-score max-sites dna/align-options)
        nnn dna/codon-length]
    (->> matches
         (map (fn [m]
                (let [{:keys [match-path score]} m
                      [start-q start-r] (first match-path)
                      [end-q end-r] (last match-path)]
                  (assoc m
                         :bind-begin (let [x start-q] (- x (mod x nnn)))
                         :bind-end-x (let [x end-q] (+ nnn (- x (mod x nnn))))
                         :bind-begin-base start-q
                         :bind-end-base end-q
                         :vs-end-base end-r)))))))

(s/fdef binding-sites
        :args (s/cat :dna ::dna/dna
                     :vs-dna ::dna/dna
                     :min-score pos-int?
                     :max-sites pos-int?)
        :ret (s/coll-of ::bind))

(defn possible-binding-sites
  "Finds all possible binding sites of a cell's DNA onto internal
  products or external stimuli. The results include a key :path which
  describes the source of each match as [kind i j]; kind is :stimuli
  or :products, i indicates which stimulus or product, and j is the
  index into its dna-like sequence. Ordered by score decreasing, and
  then by DNA offset."
  [dna products stimuli min-score max-sites]
  (->>
   (concat (map vector (repeat :products) (range) products)
           (map vector (repeat :stimuli) (range) stimuli))
   (reduce (fn [out [kind kind-i vs-dna]]
             (->> (binding-sites dna vs-dna min-score max-sites)
                  (map (fn [m]
                         (assoc m :path [kind kind-i (:vs-end-base m)])))
                  (into out)))
           ())
   (sort-by (juxt (comp - :score) :bind-end-x :vs-end-base))
   (take max-sites)))

(s/fdef possible-binding-sites
        :args (s/cat :dna ::dna/dna
                     :products (s/nilable (s/every ::dna/dna))
                     :stimuli (s/nilable (s/every ::dna/dna))
                     :min-score pos-int?
                     :max-sites pos-int?)
        :ret (s/coll-of ::bind-with-path))

(defn binding-site-weight
  [score baseline-score weight-power]
  (->
   (- score baseline-score)
   (Math/pow weight-power)))

(defn select-binding-site
  "Selects one binding site on the cell's open DNA from any matching
  products or stimuli. The result includes a :path indicating the
  match source. Selection probability is according to bind
  score (transformed). However, it is not stochastic, but simply
  divided proportionally."
  [cell stimuli time-step]
  (let [odna (get-open-dna cell)
        products (keys (:product-counts cell))
        binds (possible-binding-sites odna products stimuli (inc baseline-score)
                                      max-binding-sites)
        cum (reductions + (map #(binding-site-weight (:score %) baseline-score
                                                     weight-power)
                               binds))
        sum (last cum)]
    (when (pos? sum)
      (let [t' (mod time-step sum)
            bind-index (count (take-while #(<= % t') cum))]
        (nth binds bind-index)))))

(s/fdef select-binding-site
        :args (s/cat :cell ::cell
                     :stimuli (s/nilable (s/every ::dna/dna))
                     :time-step nat-int?)
        :ret (s/nilable ::bind-with-path))

(s/def ::offset (s/and nat-int? dna/codon-boundary?))

(s/def ::next-offset (s/or :i ::offset
                           :kw #{:stop-reaction}))

(s/def ::cell-immediate
  (s/keys :opt-un [::energy
                   ::orientation]))

(s/def ::in-direction ::orientation)
(s/def ::target ::phys/part-id)
(s/def ::sex-id ::phys/part-id)
(s/def ::child-product ::dna/dna)

(s/def ::product (s/keys :req-un [::dna/dna]))
(s/def ::clone (s/keys :req-un [::in-direction ::child-product]))
(s/def ::sex (s/keys :req-un [::in-direction ::child-product ::sex-id]))
(s/def ::bond-form (s/keys :req-un [::target]))
(s/def ::bond-break (s/keys :req-un [::target]))
(s/def ::sugar-start (s/keys :req-un [::target]))
(s/def ::sugar-stop (s/keys :req-un [::target]))
(s/def ::push (s/keys :req-un [::in-direction]))

(s/def ::delayed
  (s/keys :req-un [(or ::product
                       ::clone
                       ::sex
                       ::bond-form
                       ::bond-break
                       ::sugar-start
                       ::sugar-stop
                       ::push)]))

(s/def ::reaction-step
  (s/keys :opt-un [::next-offset
                   ::cell-immediate
                   ::delayed]))

(defmulti reaction-op*
  (fn [op cell offset cell-id phy]
    op))

(defn reaction-op
  "Computes a reaction from a single codon of DNA, interpreted into a
  symbolic op code. The offset argument here is an index into DNA
  _after_ the current op codon."
  [op cell offset cell-id phy]
  (reaction-op* op cell offset cell-id phy))

(s/fdef reaction-op
        :args (s/cat :op ::dna/op-code
                     :cell ::cell
                     :offset ::offset
                     :cell-id ::phys/part-id
                     :phy ::phys/physics)
        :ret ::reaction-step)

(defmethod reaction-op* :no-op
  [op cell offset cell-id phy]
  {})

(defmethod reaction-op* :terminator
  [op cell offset cell-id phy]
  {})

(defmethod reaction-op* :silence-terminator
  [op cell offset cell-id phy]
  {})

(defmethod reaction-op* 'stop-reaction
  [op cell offset cell-id phy]
  {:next-offset :stop-reaction})

(defn read-template
  ([dna offset]
   (read-template dna offset #(= % :terminator)))
  ([dna offset terminator?]
   (loop [dna (drop offset dna)
          tem ""]
     (if (seq dna)
       (let [codon (vec (take dna/codon-length dna))
             op (dna/codon->op codon)]
         (if (terminator? op)
           tem
           ;; otherwise, append and continue
           (recur (drop dna/codon-length dna) (apply str tem codon))))
       ;; end of DNA
       tem))))

(defmethod reaction-op* 'product
  [op cell offset cell-id phy]
  (let [open-dna (get-open-dna cell)
        cost (:product energies)
        energy (:energy cell)]
    (if (>= energy cost)
      (let [tem (read-template open-dna offset)
            terminator dna/codon-length
            tem-idx [offset (+ offset (count tem))]
            next-offset (+ offset (count tem) terminator)]
        (if (>= (count tem) min-template-bases)
          {:delayed {:product {:dna (apply str (map dna/translate tem))}}
           :next-offset next-offset
           :cell-immediate {:energy (- energy cost)}
           :read-idx {:template tem-idx}}
          ;; template too short
          {:next-offset next-offset
           :read-idx {:template tem-idx}}))
      ;; not enough energy
      {:next-offset :stop-reaction})))

(defn silence-target
  "Finds best site matching tem against dna, then reads the following
  dna until a silence-terminator. Returns this as a range of indexes
  to silence or unsilence: start inclusive, end exclusive. Returns nil
  if there was no match to the template, or if the target dna segment
  is too short."
  [dna dna-open? tem]
  ;; NOTE: use dna not open-dna because that is idempotent,
  ;; i.e. stably repeatable after silencing has taken effect
  (let [matches (binding-sites dna tem baseline-score 1)
        terminator dna/codon-length]
    (if (seq matches)
      (let [match (apply max-key :score matches)
            ;; NOTE matched template is silenced too:
            silence-start (:bind-begin match)
            target-dna (read-template dna silence-start
                                      #(= % :silence-terminator))
            silence-end (min (+ silence-start (count target-dna) terminator)
                             (count dna))]
        (if (seq target-dna)
          {:start silence-start
           :end silence-end
           :match-idx [(:bind-begin-base match) (:bind-end-base match)]}
          ;; empty target
          nil))
      ;; no match
      nil)))

(s/def :silence/start nat-int?)
(s/def :silence/end nat-int?)

(s/fdef silence-target
        :args (s/cat :dna ::dna/dna
                     :dna-open ::dna/dna-open?
                     :tem ::dna/dna)
        :ret (s/nilable (s/keys :req-un [:silence/start
                                         :silence/end])))

(defmethod reaction-op* 'silence
  [op cell offset cell-id phy]
  (let [open-dna (get-open-dna cell)
        tem (read-template open-dna offset)
        terminator dna/codon-length
        tem-idx [offset (+ offset (count tem))]
        next-offset (-> (+ offset (count tem) terminator) (min (count open-dna)))
        dna (:dna cell)
        dna-open (:dna-open? cell)]
    (if (>= (count tem) min-template-bases)
      ;; find best matching site
      (if-let [{:keys [start end match-idx]} (silence-target dna dna-open tem)]
        (let [new-dna-open (dna/set-vector-range dna-open start end false)]
          (if (>= (count (filter true? new-dna-open))
                  (* dna/min-open-dna-length-in-codons dna/codon-length))
            (let [next-offset-full (dna/offset-into-full-dna next-offset dna-open)
                  next-offset-adj (dna/offset-into-open-dna next-offset-full new-dna-open)]
              {:cell-immediate {:dna-open? new-dna-open}
               :next-offset (if (<= start next-offset-full end)
                              :stop-reaction ;; read head was silenced
                              next-offset-adj)
               :read-idx {:template tem-idx
                          :match-full match-idx
                          :silenced [start end]}
               })
            ;; open dna would be too short; abort silencing
            {:next-offset next-offset
             :read-idx {:template tem-idx}}))
        ;; no match/target
        {:next-offset next-offset
         :read-idx {:template tem-idx}})
      ;; template too short
      {:next-offset next-offset
       :read-idx {:template tem-idx}})))

(defmethod reaction-op* 'unsilence
  [op cell offset cell-id phy]
  (let [open-dna (get-open-dna cell)
        tem (read-template open-dna offset)
        terminator dna/codon-length
        tem-idx [offset (+ offset (count tem))]
        next-offset (-> (+ offset (count tem) terminator) (min (count open-dna)))]
    (if (>= (count tem) min-template-bases)
      ;; find best matching site
      (let [dna (:dna cell)
            dna-open (:dna-open? cell)]
        (if-let [{:keys [start end match-idx]} (silence-target dna dna-open tem)]
          (let [new-dna-open (dna/set-vector-range dna-open start end true)
                next-offset-full (dna/offset-into-full-dna next-offset dna-open)
                next-offset-adj (dna/offset-into-open-dna next-offset-full new-dna-open)]
            {:cell-immediate {:dna-open? new-dna-open}
             :next-offset next-offset-adj
             :read-idx {:template tem-idx
                        :match-full match-idx
                        :unsilenced [start end]}})
          ;; no match/target
          {:next-offset next-offset
           :read-idx {:template tem-idx}}))
      ;; template too short
      {:next-offset next-offset
       :read-idx {:template tem-idx}})))

(defmethod reaction-op* 'goto
  [op cell offset cell-id phy]
  (let [open-dna (get-open-dna cell)
        tem (read-template open-dna offset)
        terminator dna/codon-length
        tem-idx [offset (+ offset (count tem))]
        next-offset (+ offset (count tem) terminator)]
    (if (>= (count tem) min-template-bases)
      ;; find best matching site
      (let [matches (binding-sites open-dna tem baseline-score 1)]
        (if (seq matches)
          (let [match (apply max-key :score matches)]
            {:next-offset (:bind-end-x match)
             :read-idx {:template tem-idx
                        :match [(:bind-begin-base match) (:bind-end-base match)]}})
          ;; no match
          {:next-offset next-offset
           :read-idx {:template tem-idx}}))
      ;; template too short
      {:next-offset next-offset
       :read-idx {:template tem-idx}})))

(defmethod reaction-op* 'energy-test
  [op cell offset cell-id phy]
  (let [open-dna (get-open-dna cell)
        e-threshold (* (:max energies) (/ offset (count open-dna)))]
    (if (> (:energy cell) e-threshold)
      (reaction-op 'goto cell offset cell-id phy)
      ;; test failed, skip goto
      (let [open-dna (get-open-dna cell)
            tem (read-template open-dna offset)
            terminator dna/codon-length
            next-offset (+ offset (count tem) terminator)]
        {:next-offset next-offset}))))

(defmethod reaction-op* 'to-sun
  [op cell offset cell-id phy]
  {:cell-immediate {:orientation (/ Math/PI 2)}})

(defmethod reaction-op* 'about-face
  [op cell offset cell-id phy]
  (let [angle (:orientation cell)]
    {:cell-immediate {:orientation (in-pi-pi (+ angle Math/PI))}}))

(defmethod reaction-op* 'rot-left
  [op cell offset cell-id phy]
  (let [angle (:orientation cell)]
    {:cell-immediate {:orientation (in-pi-pi (+ angle (/ Math/PI 8)))}}))

(defmethod reaction-op* 'rot-right
  [op cell offset cell-id phy]
  (let [angle (:orientation cell)]
    {:cell-immediate {:orientation (in-pi-pi (- angle (/ Math/PI 8)))}}))

(defmethod reaction-op* 'push
  [op cell offset cell-id phy]
  (let [cost (:push energies)
        energy (:energy cell)]
    (if (>= energy cost)
      {:delayed {:push {:in-direction (:orientation cell)}}
       :cell-immediate {:energy (- energy cost)}}
      ;; not enough energy
      {:next-offset :stop-reaction})))

(defmethod reaction-op* 'sugar-start
  [op cell offset cell-id phy]
  (let [dir (:orientation cell)]
    (if-let [targ-id (phys/part-in-direction phy cell-id dir)]
      {:delayed {:sugar-start {:target targ-id}}}
      ;; no target cell
      {})))

(defmethod reaction-op* 'sugar-stop
  [op cell offset cell-id phy]
  (let [dir (:orientation cell)]
    (if-let [targ-id (phys/part-in-direction phy cell-id dir)]
      {:delayed {:sugar-stop {:target targ-id}}}
      ;; no target cell
      {})))

(defmethod reaction-op* 'bond-form
  [op cell offset cell-id phy]
  (let [cost (:bond energies)
        energy (:energy cell)
        dir (:orientation cell)]
    (if (>= energy cost)
      (if-let [targ-id (phys/part-in-direction phy cell-id dir)]
        {:delayed {:bond-form {:target targ-id}}
         :cell-immediate {:energy (- energy cost)}}
        ;; no target cell
        {})
      ;; not enough energy
      {:next-offset :stop-reaction})))

(defmethod reaction-op* 'bond-break
  [op cell offset cell-id phy]
  (let [dir (:orientation cell)]
    (if-let [targ-id (phys/part-in-direction phy cell-id dir)]
      {:delayed {:bond-break {:target targ-id}}}
      ;; no target cell
      {})))

(defmethod reaction-op* 'clone
  [op cell offset cell-id phy]
  (let [cost (:clone energies)
        energy (:energy cell)
        dir (:orientation cell)]
    (if (>= energy cost)
      (if (not (phys/part-in-direction phy cell-id dir))
        (let [prod-re (reaction-op 'product cell offset cell-id phy)
              prod (or (get-in prod-re [:delayed :product :dna])
                       (dna/fixed-stimuli :birth-default))]
          {:delayed {:clone {:in-direction (:orientation cell)
                             :child-product prod}}
           :next-offset (:next-offset prod-re)
           :cell-immediate {:energy (- energy cost)}
           :read-idx (:read-idx prod-re)})
        ;; blocked by existing cell
        {})
      ;; not enough energy
      {:next-offset :stop-reaction})))

(defmethod reaction-op* 'sex
  [op cell offset cell-id phy]
  (let [re (reaction-op 'clone cell offset cell-id phy)
        dir (:orientation cell)]
    (if-let [clone (get-in re [:delayed :clone])]
      (if-let [targ-id (phys/part-in-direction phy cell-id dir)]
        (-> re
            (assoc-in [:delayed :sex] (assoc clone :sex-id targ-id))
            (update :delayed dissoc :clone))
        ;; no sexual partner
        {})
      re)))

(defn react-at-site
  "Runs a reaction on cell's open DNA starting at the given
  offset. Interprets each codon of DNA as a symbolic operation. Each
  codon's reaction step may control the next read offset; otherwise we
  continue along the DNA. The reaction stops when one of:

  * a reaction step returns a stop instruction (e.g. stop-reaction codon);
  * we reach the end of the DNA;
  * the instruction counter reaches a limit."
  [cell bind-offset cell-id phy]
  (let [open-dna (get-open-dna cell)]
    (loop [offset bind-offset
           counter 0
           stop? false
           ret {:cell cell
                :effects []
                :reaction-log []}]
      (cond
        stop?
        ret
        (> counter max-ops)
        ret
        (>= offset (dec (count open-dna)))
        ret
        :else
        (let [codon (vec (take dna/codon-length (drop offset open-dna)))
              _ (assert (= dna/codon-length (count codon))
                        (str "Re-step " counter " offset " offset " open-dna " (count open-dna)))
              op (dna/codon->op codon)
              next-off* (+ offset dna/codon-length)
              r (reaction-op op (:cell ret) next-off* cell-id phy)
              next-off (or (:next-offset r) next-off*)
              _ (assert (or (= next-off :stop-reaction)
                            (dna/codon-boundary? next-off))
                        (str "Non codon aligned offset from " op ":" next-off
                             "/" offset ". " r))]
          (recur next-off
                 (inc counter)
                 (= next-off :stop-reaction)
                 (-> ret
                     (update :cell merge (:cell-immediate r))
                     (update :effects conj (:delayed r))
                     (update :reaction-log conj [op offset r]))))))))

(s/def ::effects (s/every (s/nilable ::delayed)))

(s/def ::reaction-log (s/every (s/tuple ::dna/op-code ::offset ::reaction-step)))

(s/def ::reaction (s/keys :req-un [::cell
                                   ::effects
                                   ::reaction-log]))

(s/fdef react-at-site
        :args (s/cat :cell ::cell
                     :bind-offset ::offset
                     :cell-id ::phys/part-id
                     :phy ::phys/physics)
        :ret ::reaction)
