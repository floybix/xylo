(ns org.nfrac.xylo.cell
  (:require [org.nfrac.xylo.dna :as dna]
            [clojure.spec.alpha :as s]
            [clojure.test.check.random :as random]))

(s/def ::dna
  (s/coll-of ::dna/base, :min-count 1))

(s/def ::dna-open?
  (s/coll-of boolean?, :min-count 1, :kind vector?))

(s/def ::energy nat-int?)

(s/def ::products
  (s/map-of ::dna nat-int?))

(s/def ::cell
  (s/keys :req [::dna
                ::dna-open?
                ::products
                ::energy]))

(defn new-cell
  [dna]
  (let []
    {::dna dna
     ::dna-open? (vec (repeat (count dna) true))
     ::products {}
     ::energy 1}))

(defn seed-cell
  []
  (new-cell "x1y2z3 qtrusv abcd=.1234?="))
