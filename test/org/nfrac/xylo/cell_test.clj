(ns org.nfrac.xylo.cell-test
  (:require [org.nfrac.xylo.cell :as cell]
            [org.nfrac.xylo.dna :as dna]
            [org.nfrac.xylo.physics :as phys]
            [org.nfrac.xylo.physics-grid :as phys-g]
            [clojure.test :as t
             :refer (is deftest testing run-tests)]))

(deftest seed-cell-test
  (let []
    (testing ""
      (let [cell (-> (cell/seed-dna) (cell/new-cell) (assoc :energy 4))
            [cell-id phy] (-> (phys-g/init 10 10) (phys/create-part [1 1]))
            dna (:dna cell)
            odna (cell/get-open-dna cell)
            cdna (apply str (map dna/complement dna))
            stim [(dna/fixed-stimuli :ground) cdna]
            binds (cell/all-binding-sites odna [] stim (inc cell/baseline-score))
            bind (cell/select-binding-site cell stim 1)]
        (doseq [b binds] (println b))
        (println "selected:")
        (println bind)
        (let [ans (cell/react-at-site cell (:bind-end-x bind) cell-id phy)]
          (println (:reaction-log ans)))))))
