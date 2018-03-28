(ns org.nfrac.xylo.sim-test
  (:require [org.nfrac.xylo.sim :as sim]
            [org.nfrac.xylo.cell :as cell]
            [clojure.test :as t
             :refer (is deftest testing run-tests)]))

(deftest world-test
  (let [world (sim/new-world 10 10)]
    (testing "First step"
      (let [w1 (sim/world-step world 1)]
        ))
    (testing "Second step"
      (let [w1 (sim/world-step world 1)
            w2 (sim/world-step w1 2)]
        ))
    ))
