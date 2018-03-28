(ns org.nfrac.xylo.physics-grid-test
  (:require [org.nfrac.xylo.physics :as phys]
            [org.nfrac.xylo.physics-grid :as phys-g]
            [clojure.test :as t
             :refer (is deftest testing run-tests)]))

(deftest physics-grid-test
  (testing "Toy physics 4 parts"
    (let [w (phys-g/init 16 12)
          [id1 w] (phys/create-part w [1.0 0.5])
          [id2 w] (phys/create-part w [5.0 1.0])
          [id3 w] (phys/create-part w [5.0 2.0])
          [id4 w] (phys/create-part w [5.5 2.0])
          w (phys/bond-form w id3 id4)]
      (println "initial world")
      (println (sort (:parts w)))
      (let [ws (iterate #(phys/step % 0.5) w)]
        (println "step 0.5")
        (println (sort (:parts (nth ws 1))))
        (println "step 0.5")
        (println (sort (:parts (nth ws 2))))
        (println "step 0.5")
        (println (sort (:parts (nth ws 3))))
        (is (= (set (vals (:parts (nth ws 3))))
               #{[1.0 0.5] [5.0 0.5] [5.0 1.0] [5.5 1.0]})
            "Parts with bonds fall into expected places.")
        ))))
