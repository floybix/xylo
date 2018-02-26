(ns org.nfrac.xylo.sim
  (:require [org.nfrac.xylo.physics :as phys]
            [org.nfrac.xylo.physics-grid :as phys-g]))

(defn test-physics-grid
  []
  (let [w (phys-g/init 16 12)
        [id1 w] (phys/create-cell w [1.0 0.5])
        [id2 w] (phys/create-cell w [5.0 1.0])
        [id3 w] (phys/create-cell w [5.0 2.0])
        [id4 w] (phys/create-cell w [5.5 2.0])
        w (phys/create-bond w id3 id4)]
    (println "initial world")
    (println (sort (:cells w)))
    (let [ws (iterate #(phys/step % 0.5) w)]
      (println "step 0.5")
      (println (sort (:cells (nth ws 1))))
      (println "step 0.5")
      (println (sort (:cells (nth ws 2))))
      (println "step 0.5")
      (println (sort (:cells (nth ws 3))))
      (= (set (vals (:cells (nth ws 3))))
         #{[1.0 0.5] [5.0 0.5] [5.0 1.0]}))))


;; react to: ground, sunlight, sugar, touching cells, products
;; (stochastically)
;; cache match score per location

;; generate energy
;; stop reaction
;; jump to matching DNA
;; conditional branching on energy level

;; silence DNA

;;

;; create product
;; products decay

;; clone
;; mutate?
;; specialisation: silencing via control branching
;; form bond

;;
