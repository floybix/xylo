(ns org.nfrac.xylo.dna
  (:require [org.nfrac.str-alignment.core :as ali]
            [clojure.spec.alpha :as s])
  (:refer-clojure :exclude [complement bases]))

(def bases
  '(a c g t))

(s/def ::base (set bases))

(def codon-length 3)

(def n-codons (Math/pow (count bases) codon-length))

(def codons
  (for [a bases
        b bases
        c bases]
    [a b c]))

(def complement (zipmap bases (drop 2 (cycle bases))))
(def translate (zipmap bases (drop 1 (cycle bases))))
(def untranslate (zipmap (vals translate) (keys translate)))

(def op-names
  '[to-stimulus
    to-sun
    about-face
    rot-left
    rot-right
    clone
    sex
    push
    bond-form
    bond-break
    sugar-start
    sugar-stop
    ;; non-physical - control flow
    product
    silence
    unsilence
    goto
    energy-test
    stop-reaction])

(def n-ops (count op-names))

(def n-terminators 12)

;; excludes terminators and no-ops (as ambiguous)
(def op->codon
  (zipmap op-names codons))

(def codon->op
  (zipmap codons (concat op-names
                         (repeat (- n-codons n-ops n-terminators) :no-op)
                         (repeat n-terminators :terminator))))
