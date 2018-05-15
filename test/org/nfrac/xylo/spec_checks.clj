(ns org.nfrac.xylo.spec-checks
  (:require [org.nfrac.xylo.sim :as sim]
            [org.nfrac.xylo.cell :as cell]
            [org.nfrac.xylo.dna :as dna]
            [org.nfrac.xylo.fancy-generators :refer [fancy-gens]]
            [clojure.spec.test.alpha :as stest]
            [clojure.test.check.clojure-test :as ctcc]
            [clojure.test :as t
             :refer (is deftest testing run-tests)]))

(alter-var-root #'ctcc/*report-shrinking* (constantly true))
(alter-var-root #'ctcc/*report-trials* (constantly ctcc/trial-report-periodic))

(alias 'stc 'clojure.spec.test.check)
(def opts {::stc/opts {:num-tests 5000}
           :gen fancy-gens
           })

(def instr-syms
  (concat
   (stest/enumerate-namespace 'org.nfrac.xylo.sim)
   (stest/enumerate-namespace 'org.nfrac.xylo.cell)
   (stest/enumerate-namespace 'org.nfrac.xylo.dna)
   ))

(deftest cell-fns-test
  (stest/instrument instr-syms)
  (-> `[dna/offset-into-open-dna
        dna/mutate-pointwise
        dna/mutate-delete
        dna/mutate-shift
        dna/mutate-dup
        dna/crossover
        sim/world-step
        ]
      (stest/check opts)
      (stest/summarize-results))
  (stest/unstrument))
