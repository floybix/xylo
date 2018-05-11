(ns org.nfrac.xylo.fancy-generators
  (:require [org.nfrac.xylo.sim :as sim]
            [org.nfrac.xylo.cell :as cell]
            [org.nfrac.xylo.dna :as dna]
            [org.nfrac.xylo.physics :as phys]
            [org.nfrac.xylo.physics-grid :as phys-g]
            [clojure.spec.alpha :as s]
            [clojure.spec.gen.alpha :as gen]
            [com.gfredericks.test.chuck :as chuck]
            [com.gfredericks.test.chuck.generators :as gen']))

(defn to-codon-boundary [i] (- i (mod i dna/codon-length)))

(defn cell-gen
  []
  (gen'/for [dna (s/gen ::dna/dna)
             :let [cell* (cell/new-cell dna)
                   n (count dna)]
             prods (s/gen ::cell/product-counts)
             energy (s/gen ::cell/energy)
             silence-from (->> (s/gen (s/int-in 0 (dec n)))
                               (gen/fmap to-codon-boundary))
             silence-to (->> (s/gen (s/int-in silence-from n))
                             (gen/fmap to-codon-boundary))
             :let [dopen (dna/set-vector-range (:dna-open? cell*) silence-from
                                               silence-to false)]]
            (assoc cell*
                   :product-counts prods
                   :energy energy
                   :dna-open? dopen)))

(defn world-gen
  []
  (gen'/for [cell1 (cell-gen)
             cell2 (cell-gen)
             rng (s/gen ::dna/rng)
             :let [phy (phys-g/init 10 5)
                   [id1 phy] (phys/create-part phy [1 0])
                   [id2 phy] (phys/create-part phy [1 1])]]
            {:physics phy
             :cell-pop {id1 cell1
                        id2 cell2}
             :sugar-from-to {}
             :rng rng}))

(defn mutate-args-gen
  []
  (gen'/for [dna (s/gen ::dna/dna)
             rng (s/gen ::dna/rng)
             :let [dna-open? (vec (repeat (count dna) true))]]
            [[dna dna-open?] rng]))

(def fancy-gens
  {::cell/cell cell-gen
   ::sim/world world-gen
   ::dna/mutate-args mutate-args-gen})
