(ns org.nfrac.xylo.physics
  (:require [clojure.spec.alpha :as s]))

(defprotocol PPhysics
  (step [this dt])
  (create-cell [this pos])
  (delete-cell [this cell-id])
  (position [this cell-id])
  (touching [this cell-id])
  (touching-ground? [this cell-id])
  (in-sunlight? [this cell-id])
  (cell-in-direction [this cell-id angle])
  (create-bond [this from-id to-id])
  (force-in-direction [this cell-id angle]))
