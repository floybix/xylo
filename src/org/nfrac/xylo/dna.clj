(ns org.nfrac.xylo.dna
  (:require [org.nfrac.str-alignment.core :as ali]
            [clojure.spec.alpha :as s]))

(def template-bases
  "abcdefghijklmnopqrstuvwxyz12345")

(def op-bases
  "+-=%:><$*@_^&?.#!")

(s/def ::base (into (set template-bases) op-bases))

(def op-names
  {\+ 'form-bond
   \- 'break-bond
   \= 'clone
   \% 'sex
   \: 'about-face
   \< 'rot-left
   \> 'rot-right
   \$ 'sugar-start
   \* 'sugar-stop
   \@ 'reference
   \_ 'silence
   \^ 'unsilence
   \& 'product
   \? 'energy-test
   \. 'end
   \# 'goto
   \! 'push
   })

(def rock-sig
  "12345")

(def sun-sig
  "abcde")

(def sugar-sig
  "5pqrs")

;; Note that the operation codes are all complements to no-op template
;; codes. The reason is that we want to be able to bind to functional
;; genes, so as to reliably recognise parasites, for example, on a
;; signal that can't be trivially changed to avoid detection. That is
;; how immune systems work.

(def complements
  (let [c0 (zipmap op-bases template-bases)
        tt (drop (count op-bases) template-bases)
        n2 (quot (count tt) 2)
        tt-rot (concat (drop n2 tt) (take n2 tt))]
    (reduce merge c0 [(zipmap (vals c0) (keys c0))
                      (zipmap tt tt-rot)
                      (zipmap tt-rot tt)])))
