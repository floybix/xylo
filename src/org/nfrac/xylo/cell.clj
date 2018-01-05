(ns org.nfrac.xylo.cell
  "A xylo cell is fundamentally a sequence of bases.
  In that respect it is actually more like an RNA molecule.  Sections
  of the sequence can be in a silenced state, like epigenetic marks on
  DNA, and this restricts the available binding sites and reactions."
  (:require [org.nfrac.xylo.dna :as dna]
            [clojure.spec :as s]
            [clojure.test.check.random :as random])))

(s/def ::dna
  (s/coll-of ::dna/base, :min-count 1))

(s/def ::dna-open
  (s/coll-of boolean?, :min-count 1, :kind vector?))))

(s/def ::energy nat-int?)

(s/def ::cell
  (s/keys :req [::dna
                ::dna-open
                ::energy]))

(defn seed-cell
  []
  (let [dna "x1y2z3 qtrusv abcd=.1234?="]
    {::dna dna
     ::dna-open (vec (repeat (count dna) true))
     ::energy 1}))
