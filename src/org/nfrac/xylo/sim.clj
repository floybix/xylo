(ns org.nfrac.xylo.sim
  (:require [org.nfrac.xylo.physics :as phys]
            [org.nfrac.xylo.physics-grid :as phys-g]
            [org.nfrac.xylo.cell :as cell]
            [org.nfrac.xylo.dna :as dna]
            [clojure.spec.alpha :as s]))

(defn test-physics-grid
  []
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
      (= (set (vals (:parts (nth ws 3))))
         #{[1.0 0.5] [5.0 0.5] [5.0 1.0] [5.5 1.0]}))))

(defn test-seed-cell
  [time-step]
  (let [cell (-> (cell/seed-cell) (assoc :energy 4))
        dna (:dna cell)
        odna (cell/open-dna (:dna cell) (:dna-open? cell))
        cdna (apply str (map dna/complement dna))
        stim [(dna/fixed-stimuli :ground) cdna]
        binds (cell/all-binding-sites odna [] stim (inc cell/baseline-score))
        bind (cell/select-binding-site cell stim time-step)]
    (doseq [b binds] (println b))
    (println "selected:")
    (println bind)
    (let [ans (cell/react-at-site cell odna (:bind-end-x bind))]
      (:reaction-log ans))
    ))


(s/def ::cell-pop
  (s/map-of ::phys/cell-id ::cell/cell))

(s/def ::physics #(satisfies? phys/PPhysics %))

(s/def ::sugar-from-to
  (s/every-kv ::phys/part-id (s/every ::phys/part-id, :kind set?)))

(s/def ::world
  (s/keys :req-un [::physics
                   ::cell-pop
                   ::sugar-from-to]))

(defn find-stimuli
  [world cell-id]
  (let [cell-pop (:cell-pop world)
        phy (:physics world)
        cell (get cell-pop cell-id)
        others (phys/touching phy cell-id)
        others-dna (map (fn [id]
                          (get-in cell-pop [id :dna]))
                        others)]
    (cond-> []
      (phys/touching-ground? phy cell-id)
      (conj (dna/fixed-stimuli :ground))
      (phys/in-sunlight? phy cell-id)
      (conj (dna/fixed-stimuli :sun))
      (seq others)
      (into (map (fn [dna]
                   (apply str (map dna/complement dna)))
                 others-dna)))))

(s/fdef find-stimuli
        :args (s/cat :world ::world
                     :cell-id ::phys/part-id)
        :ret (s/every ::dna))

(defn reaction
  "Selects a reaction partner -- either external stimulus or an internal
  product -- and a binding site on the cell. If a product is selected,
  then the returned cell has that product decremented or removed."
  [world cell-id time-step]
  (let [cell-pop (:cell-pop world)
        phy (:physics world)
        cell (get cell-pop cell-id)
        dna (:dna cell)
        odna (cell/open-dna (:dna cell) (:dna-open? cell))
        stim (find-stimuli world cell-id)
        bind (cell/select-binding-site cell stim time-step)]
    (when bind
      (let [ret (cell/react-at-site cell odna (:bind-end-x bind))
            [kind kind-i _] (:path bind)]
        (if (= kind :products)
          (let [[product n] (-> cell :product-count (nth kind-i))]
            (if (> n 1)
              (update-in ret [:cell :product-count product] dec)
              (update-in ret [:cell :product-count] dissoc product))))))))

(s/fdef reaction
        :args (s/cat :world ::world
                     :cell-id ::phys/part-id
                     :time-step nat-int?)
        :ret (s/nilable ::cell/reaction))

(defn apply-reaction-effect
  [world cell-id effect]
  (let [cell-pop (:cell-pop world)
        phy (:physics world)
        cell (get cell-pop cell-id)
        [fx fx-args] (first effect)]
    (case fx

      :product
      (let [prod fx-args]
        (update-in world [:cell-pop cell-id :product-count prod]
                   (fnil inc 0)))

      :clone
      (let [dir (:in-direction fx-args)
            [x y] (phys/position phy cell-id)
            [dx dy] (phys-g/polar-xy 1.0 dir)
            new-pos [(+ x dx) (+ y dy)]
            child-product (:child-product fx-args)
            [new-id phy] (phys/create-part phy new-pos)
            dna (:dna cell)
            new-dna dna ;; TODO mutation
            new-dna-open? (:dna-open? cell)
            new-cell (-> (cell/new-cell new-dna)
                         (assoc :dna-open? new-dna-open?)
                         (assoc :orientation (:orientation cell))
                         (assoc-in [:product-counts child-product] 1))]
        (-> world
            (assoc :physics phy)
            (update :cell-pop assoc new-id new-cell)))

      :bond-form
      (let [dir (:in-direction fx-args)]
        (if-let [targ-id (phys/part-in-direction phy cell-id dir)]
          (update world :physics phys/bond-form cell-id targ-id)
          world))

      :bond-break
      (let [dir (:in-direction fx-args)]
        (if-let [targ-id (phys/part-in-direction phy cell-id dir)]
          (update world :physics phys/bond-break cell-id targ-id)
          world))

      :push
      (let [dir (:in-direction fx-args)]
        (update world :physics phys/force-in-direction cell-id dir))

      :sugar-start
      (let [dir (:in-direction fx-args)]
        (if-let [targ-id (phys/part-in-direction phy cell-id dir)]
          (update-in world [:sugar-from-to cell-id] #(conj (or % #{}) targ-id))
          world))

      :sugar-stop
      (let [dir (:in-direction fx-args)]
        (if-let [targ-id (phys/part-in-direction phy cell-id dir)]
          (update-in world [:sugar-from-to cell-id] disj targ-id)
          world))

      )))

(defn world-step
  [world time-step]
  (let [cell-pop (:cell-pop world)]
    ;; TODO sunlight / sugar
    (reduce (fn [world cell-id]
              (if-let [re (reaction world cell-id time-step)]
                (reduce (fn [w effect]
                          (apply-reaction-effect w cell-id effect))
                        (assoc-in world [:cell-pop cell-id] (:cell re))
                        (remove nil? (:effects re)))
                world))
            world
            (sort (keys cell-pop)))))

(defn new-world
  [width height]
  (let [phy (phys-g/init width height)
        cell (cell/seed-cell)
        [cell-id phy] (phys/create-part phy [(quot width 2) 0])]
    {:physics phy
     :cell-pop {cell-id cell}
     :sugar-from-to {}}))

;; cache match score per location

;; silence DNA

;; products decay

;; mutate
;; specialisation: silencing via control branching
