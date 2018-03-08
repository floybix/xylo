(ns org.nfrac.xylo.physics
  (:require [clojure.spec.alpha :as s]))

(defprotocol PPhysics
  (step [this dt])
  (create-part [this pos])
  (delete-part [this part-id])
  (position [this part-id])
  (touching [this part-id])
  (touching-ground? [this part-id])
  (in-sunlight? [this part-id])
  (part-in-direction [this part-id angle])
  (create-bond [this from-id to-id])
  (force-in-direction [this part-id angle]))
