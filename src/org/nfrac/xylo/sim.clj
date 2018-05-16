(ns org.nfrac.xylo.sim
  (:require [org.nfrac.xylo.physics :as phys]
            [org.nfrac.xylo.physics-grid :as phys-g]
            [org.nfrac.xylo.cell :as cell]
            [org.nfrac.xylo.dna :as dna]
            [clojure.test.check.random :as random]
            [clojure.spec.alpha :as s]))

(def physics-dt 1.0)

(s/def ::cell-pop
  (s/map-of ::phys/part-id ::cell/cell))

(s/def ::sugar-from-to
  (s/every-kv ::phys/part-id (s/every ::phys/part-id, :kind set?)))

(s/def ::world
  (s/keys :req-un [::phys/physics
                   ::cell-pop
                   ::sugar-from-to
                   ::dna/rng]))

(defn find-stimuli
  [world cell-id touch-ids]
  (let [cell-pop (:cell-pop world)
        phy (:physics world)
        cell (get cell-pop cell-id)
        [x y] (phys/position phy cell-id)]
    (cond-> []
      (phys/touching-ground? phy cell-id)
      (conj {:ident :ground
             :dna (dna/fixed-stimuli :ground)
             :orientation (- (/ Math/PI 2))})
      (phys/in-sunlight? phy cell-id)
      (conj {:ident :sun
             :dna (dna/fixed-stimuli :sun)
             :orientation (/ Math/PI 2)})
      (seq touch-ids)
      (into (map (fn [id]
                   (let [cell-i (get cell-pop id)
                         odna (cell/get-open-dna cell-i)
                         cdna (apply str (map dna/complement odna))
                         [xi yi] (phys/position phy id)
                         angle (phys-g/v-angle [(- xi x) (- yi y)])]
                     {:ident id
                      :dna cdna
                      :orientation angle}))
                 touch-ids)))))

(s/def ::stimulus
  (s/keys :req-un [::dna/dna
                   ::cell/orientation]))

(s/fdef find-stimuli
        :args (s/cat :world ::world
                     :cell-id ::phys/part-id
                     :touch-ids (s/every ::phys/part-id))
        :ret (s/every ::stimulus))

(defn cell-reaction
  "Selects a reaction partner -- either external stimulus or an internal
  product -- and a binding site on the cell. If a product is selected,
  then the returned cell has that product decremented or removed."
  [world cell-id touch-ids time-step]
  (let [cell-pop (:cell-pop world)
        phy (:physics world)
        cell (get cell-pop cell-id)
        stim (find-stimuli world cell-id touch-ids)
        bind (cell/select-binding-site cell (map :dna stim) time-step)]
    (when bind
      (let [[kind kind-i _] (:path bind)
            cell (cond-> cell
                   (= kind :stimuli)
                   (assoc :orientation (:orientation (nth stim kind-i))))
            ret (-> (cell/react-at-site cell (:bind-end-x bind) cell-id phy)
                    (assoc :stimuli stim)
                    (assoc ::cell/path (:path bind)))]
        (if (= kind :products)
          (let [[product n] (-> cell :product-counts (seq) (nth kind-i))]
            (if (> n 1)
              (update-in ret [:cell :product-counts product] dec)
              (update-in ret [:cell :product-counts] dissoc product)))
          ret)))))

(s/fdef cell-reaction
        :args (s/cat :world ::world
                     :cell-id ::phys/part-id
                     :touch-ids (s/every ::phys/part-id)
                     :time-step nat-int?)
        :ret (s/nilable ::cell/reaction))

(defn- abs [x] (if (neg? x) (- x) x))

(defn intruding
  [phy part-id]
  (let [within 0.4
        [x y] (phys/position phy part-id)
        parts (:parts phy)]
    (->> (dissoc parts part-id)
         (keep (fn [[id [ix iy]]]
                 (when (and (<= (abs (- x ix)) within)
                            (<= (abs (- y iy)) within))
                   id))))))

(defn apply-reaction-effect
  [world cell-id effect]
  (let [cell-pop (:cell-pop world)
        phy (:physics world)
        cell (get cell-pop cell-id)
        [fx fx-args] (first effect)]
    (case fx

      :product
      (let [prod (:dna fx-args)]
        (update-in world [:cell-pop cell-id :product-counts prod]
                   (fnil inc 0)))

      :clone
      (let [dir (:in-direction fx-args)
            [x y] (phys/position phy cell-id)
            [dx dy] (phys-g/polar-xy 1.0 dir)
            new-pos [(+ x dx) (+ y dy)]
            child-product (:child-product fx-args)
            [new-id phy] (phys/create-part phy new-pos)
            [rng rng*] (random/split (:rng world))
            [new-dna new-dna-open] (dna/mutate [(:dna cell) (:dna-open? cell)] rng*)
            new-cell (-> (cell/new-cell new-dna)
                         (assoc :dna-open? new-dna-open)
                         (assoc :orientation (:orientation cell))
                         (assoc-in [:product-counts child-product] 1))]
        ;; a new cell might have been created in the same place; if so, abort
        (if (seq (intruding phy new-id))
          world
          (-> world
              (assoc :physics phy)
              (assoc :rng rng)
              (update :cell-pop assoc new-id new-cell))))

      :sex
      (let [dir (:in-direction fx-args)
            sex-id (:sex-id fx-args)
            [x y] (phys/position phy cell-id)
            [dx dy] (phys-g/polar-xy 1.0 (+ dir Math/PI))
            new-pos [(+ x dx) (+ y dy)]
            child-product (:child-product fx-args)
            [new-id phy] (phys/create-part phy new-pos)
            sex-cell (get cell-pop sex-id)
            dna (:dna cell)
            sex-dna (:dna sex-cell)
            [rng rng*] (random/split (:rng world))
            new-dna (dna/crossover dna sex-dna rng*)
            new-cell (-> (cell/new-cell new-dna)
                         (assoc :orientation (:orientation cell))
                         (assoc-in [:product-counts child-product] 1))]
        ;; a new cell might have been created in the same place; if so, abort
        (if (seq (intruding phy new-id))
          world
          (-> world
             (assoc :physics phy)
             (assoc :rng rng)
             (update :cell-pop assoc new-id new-cell))))

      :bond-form
      (let [targ-id (:target fx-args)]
        (update world :physics phys/bond-form cell-id targ-id))

      :bond-break
      (let [targ-id (:target fx-args)]
        (update world :physics phys/bond-break cell-id targ-id))

      :push
      (let [dir (:in-direction fx-args)]
        (update world :physics phys/force-in-direction cell-id dir))

      :sugar-start
      (let [targ-id (:target fx-args)]
        (update-in world [:sugar-from-to cell-id] #(conj (or % #{}) targ-id)))

      :sugar-stop
      (let [targ-id (:target fx-args)]
        (update-in world [:sugar-from-to cell-id] #(disj (or % #{}) targ-id)))

      )))

(defn init-step
  "Builds caches shared by sub step functions."
  [world]
  (let [phy (:physics world)]
    (assoc world ::touching-cache
           (into {}
                 (map (fn [cell-id]
                        [cell-id (set (phys/touching phy cell-id))]))
                 (keys (:cell-pop world))))))

(s/fdef init-step
        :args (s/cat :world ::world)
        :ret ::world)

(defn sun-step
  "Accumulates solar energy in cells receiving direct sunlight."
  [world]
  (let [phy (:physics world)
        max-energy (:max cell/energies)]
    (reduce (fn [world cell-id]
              (cond-> world
                (phys/in-sunlight? phy cell-id)
                (update-in [:cell-pop cell-id :energy]
                           #(min max-energy (+ % (:sun cell/energies))))))
            world
            (keys (:cell-pop world)))))

(s/fdef sun-step
        :args (s/cat :world ::world)
        :ret ::world)

(defn sugar-step
  "Transfers energy between cells along sugar channels. Channels are
  only open when the cells are also touching. The amount transferred
  is limited by availability on the supply side and capacity on the
  demand side."
  [world]
  (let [sugar-from-to (:sugar-from-to world)
        touching (::touching-cache world)
        max-energy (:max cell/energies)]
    (reduce (fn [world cell-id]
              (let [sugar-to (->> (sugar-from-to cell-id)
                                  (filter (touching cell-id)))]
                (loop [sugar-to sugar-to
                       cell-pop (:cell-pop world)
                       my-e (get-in world [:cell-pop cell-id :energy])]
                  (if-let [to (first sugar-to)]
                    (let [to-e (get-in cell-pop [to :energy])
                          diff-e (-> (:sugar cell/energies)
                                     (min my-e)
                                     (min (- max-energy to-e)))]
                      (recur (rest sugar-to)
                             (assoc-in cell-pop [to :energy] (+ to-e diff-e))
                             (- my-e diff-e)))
                    ;; done
                    (assoc world :cell-pop
                           (assoc-in cell-pop [cell-id :energy] my-e))))))
            world
            (->> (:cell-pop world)
                 (sort-by #(:energy (val %)) >)
                 (map key)))))

(s/fdef sugar-step
        :args (s/cat :world ::world)
        :ret ::world)

(defn metabolism-step
  "Consumes 1 energy in all cells representing metabolic requirements of
  life. Increments the starvation counter on cells without enough
  energy. If a cell has been starving for long enough it dies and is
  removed from the simulation, along with any bonds and sugar
  channels. When a starving cell has energy again, it uses as much as
  it can to pay off the metabolic deficit incurred by starvation."
  [world]
  (reduce (fn [world cell-id]
            (let [cell (get-in world [:cell-pop cell-id])
                  starv (:starvation cell)]
              (if (pos? (:energy cell))
                (let [e-use (-> (:starvation cell) (min (:energy cell)) (max 1))]
                  (-> world
                      (update-in [:cell-pop cell-id :energy] - e-use)
                      (update-in [:cell-pop cell-id :starvation] #(max 0 (- % e-use)))))
                ;; starving
                (if (< starv cell/starvation-steps)
                  (update-in world [:cell-pop cell-id :starvation] inc)
                  ;; death
                  (-> world
                      (update :physics phys/delete-part cell-id)
                      (update :cell-pop dissoc cell-id)
                      (update :sugar-from-to dissoc cell-id))))))
          world
          (keys (:cell-pop world))))

(s/fdef metabolism-step
        :args (s/cat :world ::world)
        :ret ::world)

(defn pre-reaction-steps
  [world]
  (-> world
      (update :physics phys/step physics-dt)
      (init-step)
      (sun-step)
      (sugar-step)
      (metabolism-step)
      (init-step) ;; update for any destroyed cells
      ))

(defn reaction-step
  "Selects a reaction for each cell, then applies all the reaction effects."
  [world]
  (let [touching (::touching-cache world)
        time-step (inc (:time-step world 0))
        res (map (fn [cell-id]
                   [cell-id
                    (cell-reaction world cell-id (touching cell-id) time-step)])
                 (keys (:cell-pop world)))]
    (reduce (fn [world [cell-id re]]
              (if re
                (reduce (fn [w effect]
                          (apply-reaction-effect w cell-id effect))
                        (assoc-in world [:cell-pop cell-id] (:cell re))
                        (remove nil? (:effects re)))
                world))
            (assoc world :time-step time-step)
            res)))

(s/fdef reaction-step
        :args (s/cat :world ::world)
        :ret ::world)

(defn world-step
  [world]
  (-> world
      (pre-reaction-steps)
      (reaction-step)))

(s/fdef world-step
        :args (s/cat :world ::world)
        :ret ::world)

(defn new-world
  [width height sim-seed init-dna]
  (let [phy (phys-g/init width height)
        cell (cell/new-cell init-dna)
        [cell-id phy] (phys/create-part phy [(quot width 2) 0])
        rng (random/make-random sim-seed)]
    {:physics phy
     :cell-pop {cell-id cell}
     :sugar-from-to {}
     :rng rng}))

(s/fdef new-world
        :args (s/cat :width pos?
                     :height pos?
                     :sim-seed integer?
                     :init-dna ::dna/dna)
        :ret ::world)

;; cache match score per location

;; silence DNA

;; products decay

;; mutate
;; specialisation: silencing via control branching
