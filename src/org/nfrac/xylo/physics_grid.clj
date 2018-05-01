(ns org.nfrac.xylo.physics-grid
  (:require [org.nfrac.xylo.physics :as phys]
            [loom.graph :as graph]
            [loom.alg :as alg]
            [clojure.spec.alpha :as s]))

(s/def ::pos (s/tuple (s/double-in ::min 0 ::NaN? false)
                      (s/double-in ::min 0 ::NaN? false)))
(s/def ::parts (s/every-kv ::phys/part-id ::pos))
(s/def ::width (s/double-in ::min 1 ::NaN? false))
(s/def ::height (s/double-in ::min 1 ::NaN? false))
(s/def ::bonds (s/every-kv ::phys/part-id (s/every ::phys/part-id, :kind set?)))
(s/def ::parts-list (s/every ::phys/part-id))

(s/def ::world
  (s/keys :req-un [::width
                   ::height
                   ::parts
                   ::bonds]))

(def PI Math/PI)
(def TWOPI (* 2.0 PI))

(defn in-pi-pi
  "Returns the angle expressed in the range -pi to pi."
  [angle]
  (cond
    (> angle PI) (in-pi-pi (- angle TWOPI))
    (< angle (- PI)) (in-pi-pi (+ angle TWOPI))
    :else angle))

(defn polar-xy
  "Convert polar coordinates (magnitude, angle) to cartesian
   coordinates (x, y)."
  [mag angle]
  [(* mag (Math/cos angle))
   (* mag (Math/sin angle))])

(defn v-angle
  "Angle of a 2d geometric vector in radians in range -pi to pi."
  [[x y]]
  (Math/atan2 y x))

(defn- abs [x] (if (neg? x) (- x) x))

(defn grounding
  "Returns connected components -- collections of part ids bonded
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

(s/fdef grounding
        :args (s/cat :ccs ::parts-list
                     :supp-id? (s/fspec :args (s/cat :id ::phys/part-id)))
        :ret (s/cat :supp ::parts-list
                    :un ::parts-list))

(defn on-ground?
  [[x y]]
  (< y 1.0))

(defn- sitting-on?
  [[x y] supp-xys]
  (->> supp-xys
       (some (fn [[ix iy]]
               (and (<= (abs (- ix x)) 1.0)
                    (<= 0.0 (- y iy) 1.0))))))

(defn unsupported-ids
  [parts bonds]
  (let [singletons (keys (apply dissoc parts (keys bonds)))
        g (apply graph/graph (cond-> singletons
                               (seq bonds) (conj bonds)))
        ccs (alg/connected-components g)]
    (loop [[supp un] (grounding ccs #(on-ground? (parts %)))]
      (if (seq un)
        ;; look for components sitting directly on top of supported ones
        (let [supp-xys (map parts (apply concat supp))
              [newly-supp new-un] (grounding un #(sitting-on? (parts %) supp-xys))
              ]
          (if (seq newly-supp)
            (recur [(into supp newly-supp) new-un])
            ;; nothing changed so we are done
            (apply concat un)))
        ;; nothing unsupported
        nil))))

(s/fdef unsupported-ids
        :args (s/cat :parts ::parts
                     :bonds ::bonds)
        :ret (s/nilable ::parts-list))

(defrecord GridWorld [width
                      height
                      parts
                      bonds]
  phys/PPhysics

  (step
   [this dt]
   (let [uns (unsupported-ids parts bonds)]
     (reduce (fn [this id]
               (update-in this [:parts id 1] - dt))
             this
             uns)))

  (create-part
   [this pos]
   (let [id (gensym "p")]
     [id (assoc-in this [:parts id] pos)]))

  (delete-part
    [this part-id]
    (let [bonded (get bonds part-id)
          bonds (reduce (fn [bonds other-id]
                          (update bonds other-id disj part-id))
                        bonds
                        bonded)]
      (-> (assoc this :bonds bonds)
          (update :bonds dissoc part-id)
          (update :parts dissoc part-id))))

  (position
   [this part-id]
   (get parts part-id))

  (touching
   [this part-id]
   (when-let [[x y] (get parts part-id)]
     (->> (dissoc parts part-id)
          (keep (fn [[id [ix iy]]]
              (when (and (<= (abs (- x ix)) 1.0)
                         (<= (abs (- y iy)) 1.0))
                id))))))

  (touching-ground?
   [this part-id]
   (when-let [[x y] (get parts part-id)]
     (< y 1.0)))

  (in-sunlight?
   [this part-id]
   (when-let [[x y] (get parts part-id)]
     (->> (dissoc parts part-id)
          (vals)
          (not-any? (fn [[ix iy]]
                      (and (< (abs (- x ix)) 1.0)
                           (> iy y)))))))

  (part-in-direction
   [this part-id angle]
   (let [[x y] (get parts part-id)
         _ (assert x (str "Part ID not found " part-id))
         [dx dy] (polar-xy 0.5 angle)
         [ax ay] [(+ x dx) (+ y dy)]]
     (->> (dissoc parts part-id)
          (keep (fn [[id [ix iy]]]
                  (let [d (Math/sqrt (+ (Math/pow (- ix ax) 2)
                                        (Math/pow (- iy ay) 2)))]
                    (when (<= d 0.5)
                      [id d]))))
          (sort-by second)
          (ffirst))))

  (bond-form
    [this from-id to-id]
    (-> this
        (update-in [:bonds from-id] #(conj (or % #{}) to-id))
        (update-in [:bonds to-id] #(conj (or % #{}) from-id))))

  (bond-break
    [this from-id to-id]
    (-> this
        (update-in [:bonds from-id] disj to-id)
        (update-in [:bonds to-id] disj from-id)))

  (force-in-direction
   [this part-id angle]
   this)
  )

(defn init
  [width height]
  (map->GridWorld
   {:width width
    :height height
    :parts {}
    :bonds {}}))

(s/fdef init
        :args (s/cat :width pos?
                     :height pos?)
        :ret ::phys/physics)
