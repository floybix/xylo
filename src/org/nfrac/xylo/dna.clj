(ns org.nfrac.xylo.dna
  (:require [org.nfrac.str-alignment.core :as ali]
            [clojure.spec.alpha :as s])
  (:refer-clojure :exclude [complement bases]))

(def bases
  (seq "acgt"))

(def base? (set bases))

(s/def ::base base?)

(def codon-length 3)

(def n-codons (int (Math/pow (count bases) codon-length)))

(def codons
  (for [a bases
        b bases
        c bases]
    [a b c]))

(def complement (zipmap bases (drop 2 (cycle bases))))
(def translate (zipmap bases (drop 1 (cycle bases))))
(def untranslate (zipmap (vals translate) (keys translate)))

(def op-names
  '[to-sun
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
(def n-silence-terminators 6)
(def n-terminators 12)

;; excludes terminators and no-ops (as ambiguous)
(def op->codon
  (zipmap op-names codons))

(def codon->op
  (zipmap codons (concat op-names
                         (repeat n-terminators :terminator)
                         (repeat n-silence-terminators :silence-terminator)
                         (repeat :no-op))))

(s/def ::op-code (set (vals codon->op)))

(def terminators (keep (fn [[codon op]] (when (= op :terminator) codon)) codon->op))
(def no-ops (keep (fn [[codon op]] (when (= op :no-op) codon)) codon->op))

(def fixed-stimuli
  (zipmap [:sun :ground :birth-default]
          (map #(apply str (apply concat %))
               (partition 2 2 no-ops))))
