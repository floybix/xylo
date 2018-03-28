(ns org.nfrac.xylo.spec-tests
  (:require [org.nfrac.xylo.sim :as sim]
            [org.nfrac.xylo.cell :as cell]
            [clojure.spec.test.alpha :as stest]
            [clojure.test :as t
             :refer (is deftest testing run-tests)]))

(stest/instrument
 (concat
  (stest/enumerate-namespace 'org.nfrac.xylo.sim)
  (stest/enumerate-namespace 'org.nfrac.xylo.cell)
  (stest/enumerate-namespace 'org.nfrac.xylo.dna)))

(run-tests 'org.nfrac.xylo.cell-test
           'org.nfrac.xylo.sim-test
           'org.nfrac.xylo.physics-grid-test
           )

(stest/unstrument)
