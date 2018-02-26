(ns org.nfrac.xylo.physics-grid
  (:require [org.nfrac.xylo.physics :as phys]
            [loom.graph :as graph]
            [loom.alg :as alg]
            [clojure.spec.alpha :as s]))

(s/def ::pos (s/tuple (s/double-in ::min 0 ::NaN? false)
                      (s/double-in ::min 0 ::NaN? false)))
(s/def ::cells (s/every-kv ::phys/cell-id ::pos))
(s/def ::width (s/double-in ::min 1 ::NaN? false))
(s/def ::height (s/double-in ::min 1 ::NaN? false))
(s/def ::bonds (s/every-kv ::phys/cell-id (s/every ::phys/cell-id)))

(s/def ::world
  (s/keys :req-un [::width
                   ::height
                   ::cells
                   ::bonds]))

(defn- abs [x] (if (neg? x) (- x) x))

(defn grounding
  "Returns connected components -- collections of cell ids bonded
  together -- in a tuple of two lists. Those in the first list are
  supported, those in the second list aren't."
  [ccs supp-id?]
  (loop [ccs ccs
         supp ()
         un ()]
    (if-let [ids (first ccs)]
      (if (some supp-id? ids)
        (recur (rest ccs) (conj supp ids) un)
        (recur (rest ccs) supp (conj un ids)))
      ;; done
      [supp un])))

(defn on-ground?
  [[x y]]
  (< y 1.0))

(defn- sitting-on?
  [[x y] supp-xys]
  (->> supp-xys
       (some (fn [[ix iy]]
               (and (< (abs (- ix x)) 1.0)
                    (< 0.0 (- y iy) 1.0))))))

(defn unsupported-ids
  [cells bonds]
  (let [singletons (keys (apply dissoc cells (keys bonds)))
        g (apply graph/graph bonds singletons)
        ccs (alg/connected-components g)]
    (loop [[supp un] (grounding ccs #(on-ground? (cells %)))]
      (if (seq un)
        ;; look for components sitting directly on top of supported ones
        (let [supp-xys (map cells (apply concat supp))
              [newly-supp new-un] (grounding un #(sitting-on? (cells %) supp-xys))
              ]
          (if (seq newly-supp)
            (recur [(into supp newly-supp) new-un])
            ;; nothing changed so we are done
            (apply concat un)))
        ;; nothing unsupported
        nil))))

(defrecord GridWorld [width
                      height
                      cells
                      bonds]
  phys/PPhysics

  (step
   [this dt]
   (let [uns (unsupported-ids cells bonds)]
     (reduce (fn [this id]
               (update-in this [:cells id 1] - dt))
             this
             uns)))

  (create-cell
   [this pos]
   (let [id (gensym "cell")]
     [id (assoc-in this [:cells id] pos)]))

  (delete-cell
   [this cell-id]
   (let [bonded (get bonds cell-id)]
     (-> (assoc this :bonds (apply dissoc bonds bonded))
         (update :bonds dissoc cell-id)
         (update :cells dissoc cell-id))))

  (position
   [this cell-id]
   (get cells cell-id))

  (touching
   [this cell-id]
   (when-let [[x y] (get cells cell-id)]
     (->> (dissoc cells cell-id)
          (keep (fn [[id [ix iy]]]
              (when (and (< (abs (- x ix)) 1.0)
                         (< (abs (- y iy)) 1.0))
                id))))))

  (touching-ground?
   [this cell-id]
   (when-let [[x y] (get cells cell-id)]
     (< y 1.0)))

  (in-sunlight?
   [this cell-id]
   (when-let [[x y] (get cells cell-id)]
     (->> (dissoc cells cell-id)
          (vals)
          (not-any? (fn [[ix iy]]
                      (and (< (abs (- x ix)) 1.0)
                           (> iy y)))))))

  (cell-in-direction
   [this cell-id angle]
   (let [[x y] (get cells cell-id)
         dx (* 0.5 (Math/cos angle))
         dy (* 0.5 (Math/sin angle))
         [ax ay] [(+ x dx) (+ y dy)]]
     (->> (dissoc cells cell-id)
          (keep (fn [[id [ix iy]]]
                  (let [d (+ (Math/pow (- ix ax) 2)
                             (Math/pow (- iy ay) 2))]
                    (when (< d 0.5)
                      [id d]))))
          (sort-by second)
          (ffirst))))

  (create-bond
   [this from-id to-id]
   (-> this
       (update-in [:bonds from-id] conj to-id)
       (update-in [:bonds to-id] conj from-id)))

  (force-in-direction
   [this cell-id angle]
   this)
  )

(defn init
  [width height]
  (map->GridWorld
   {:width width
    :height height
    :cells {}
    :bonds {}}))
