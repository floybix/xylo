(ns org.nfrac.xylo.viz
  (:require [org.nfrac.xylo.sim :as sim]
            [org.nfrac.xylo.physics :as phys]
            [org.nfrac.xylo.cell :as cell]
            [org.nfrac.xylo.dna :as dna]
            [org.nfrac.str-alignment.core :as ali])
  )

(defn with-vega-meta
  [m]
  (merge
   {"$schema" "https://vega.github.io/schema/vega/v3.0.json"
    :width   600
    :height  400
    :padding {:top 50, :left 50, :bottom 50, :right 50}
    }
   m))

(defn space-viz*
  [cell-data bonds-data sugar-data [width height] width-px]
  {:width width-px
   :height (* width-px (/ height width))
   :scales [{:name   "x",
             :type   "linear",
             :range  "width",
             :domain [0 width],
             :nice   true
             :round  true}
            {:name   "y",
             :type   "linear",
             :range  "height",
             :domain [0 height],
             :nice   true
             :round  true}
            {:name   "fill",
             :type   "quantize",
             :domain [0 (:max cell/energies)],
             :range  {:scheme "reds", :count (min 10 (int (:max cell/energies)))}}]
   :axes   [{:scale "x", :orient "bottom"}
            {:scale "y", :orient "left"}]
   :data   [{:name "cells"
             :values cell-data}
            {:name "bonds"
             :values bonds-data}]
   :marks [{:type "arc",
            :from {:data "cells"},
            :encode
            {:enter
             {:x           {:scale "x", :field "x"},
              :y           {:scale "y", :field "y"},
              :fill        {:scale "fill", :field "energy"}
              :stroke      {:value "black"}
              :startAngle  {:field "orientation"}
              :endAngle    {:field "lifeFrac" :mult (* 1.99 Math/PI)
                            :offset {:field "orientation"} }
              :outerRadius {:scale "x", :value 0.5}
              :fillOpacity {:value 1}},
             :hover {:fillOpacity {:value 1}}}}
           {:type "text"
            :from {:data "cells"}
            :encode
            {:enter
             {:x {:scale "x", :field "x"}
              :y {:scale "y", :field "y"}
              :text {:field "id"}
              :fontSize {:scale "x", :value 0.3}
              :align {:value "center"}}}}]
   :legends [{:fill "fill"
              :title "Energy"}]
   })

(defn space-viz
  [world width-px]
  (let [phy (:physics world)
        cell-data (for [[id cell] (sort (:cell-pop world))]
                    (let [[x y] (phys/position phy id)
                          z (:starvation cell)]
                      (-> (select-keys cell [:energy
                                             :starvation
                                             :orientation])
                          (assoc :x x
                                 :y y
                                 :lifeFrac (- 1.0 (/ z cell/starvation-steps))
                                 :id (name id)))))
        bonds-data (for [[from tos] (sort (:bonds (:physics world)))
                         to (sort tos)]
                     {:from from
                      :to to})
        sugar-data (for [[from tos] (sort (:sugar-from-to world))
                         to (sort tos)]
                     {:from from
                      :to to})]
    (space-viz* (vec cell-data) (vec bonds-data) (vec sugar-data)
              [(:width phy) (:height phy)]
              width-px)))

(defn silenced-runs
  [dna-open?]
  (loop [dna-open? dna-open?
         i 0
         run-start nil
         runs []]
    (if (seq dna-open?)
      (let [o? (first dna-open?)
            start? (and o? (nil? run-start))
            end? (and (not o?) run-start)]
        (recur (rest dna-open?)
               (inc i)
               (cond start? i
                     end? nil
                     :else run-start)
               (cond-> runs
                 end? (conj [run-start i]))))
      ;; done
      (cond-> runs
        run-start (conj [run-start i])))))

(defn basic-dna-viz*
  [dna-data silence-data row-scale-name]
  {:scales [{:name   "base",
             :type   "linear",
             :range  "width",
             :domain [0 (count dna-data)]}
            {:name   "codon",
             :type   "linear",
             :range  "width",
             :domain [0 (float (/ (count dna-data) dna/codon-length))]}]
   :axes [{:scale "codon"
           :orient "top"
           :grid true
           :ticks false
           :offset {:scale row-scale-name, :band 1}}
          {:scale "codon"
           :orient "top"
           :grid true
           :ticks false
           :labels false}]
   :data [{:name "dna"
           :values dna-data}
          {:name "codon-bars"
           :values (vec (for [i (range 0 (/ (count dna-data) dna/codon-length))
                              :when (odd? i)]
                          {:codon i}))}
          {:name "silence"
           :values silence-data}]
   :marks [{:type "rect"
            :from {:data "silence"}
            :zindex 1
            :encode
            {:enter
             {:x           {:scale "base", :field "begin"}
              :x2          {:scale "base", :field "end"}
              :y           {:value 0}
              :height      "height"
              :fill        {:value "#aaaaaa"}}}}
           {:type "rect"
            :from {:data "codon-bars"}
            :zindex 0
            :encode
            {:enter
             {:x           {:scale "codon", :field "codon"}
              :width       {:scale "codon", :value 1.0}
              :y           {:scale row-scale-name, :band -1}
              :height      {:signal "height", :offset {:scale row-scale-name, :band 1}}
              :fill        {:value "black"}
              :fillOpacity {:value 0.05}}}}
           {:type "text",
            :from {:data "dna"},
            :zindex 3
            :encode
            {:enter
             {:x           {:scale "base", :field "at"}
              :y           {:scale row-scale-name, :band -0.5}
              :text        {:field "base"}
              :font        {:value "monospace"}
              :fontSize    {:scale row-scale-name, :band 0.4}
              :baseline    {:value "middle"}
              }}}]})

(defn basic-dna-viz
  [dna dna-open? row-scale-name]
  (let [odna (dna/open-dna dna dna-open?)
        dna-data (for [[i base] (map vector (range) dna)]
                   {:at i
                    :base base})
        silence-data (for [[begin end] (silenced-runs dna-open?)]
                       {:begin begin, :end end})]
    (-> (basic-dna-viz* dna-data silence-data row-scale-name)
        (with-vega-meta))))

(defn binds-viz*
  [dna-vega binds-data align-data [width-px height-px]]
  (merge-with
   into
   (assoc dna-vega
          :width width-px
          :height height-px
          :padding {:left 20, :right 150, :top 20, :bottom 20})
   {:scales [{:name   "bind",
              :type   "band"
              :range  "height",
              :domain {:data "binds", :field "name"
                       :sort {:field "bidx", :order "ascending"}}
              :padding 0.1}
             {:name   "bind-by-index",
              :type   "band"
              :range  "height",
              :domain {:data "binds", :field "bidx"
                       :sort {:field "bidx", :order "ascending"}}
              :padding 0.1}]
    :axes [{:scale "bind"
            :orient "left"}]
    :data [{:name "binds"
            :values binds-data}
           {:name "align"
            :values align-data}]
    :marks [{:type "rect"
             :from {:data "binds"}
             :zindex 2
             :encode
             {:enter
              {:x           {:scale "base", :field "begin"}
               :x2          {:scale "base", :field "end"}
               :y           {:scale "bind", :field "name"}
               :height      {:scale "bind", :band 1.0}
               :fill        {:value "#9999ff"}
               :fillOpacity {:value 0.5}}}}
            {:type "text",
             :from {:data "align"},
             :zindex 3
             :encode
             {:enter
              {:x           {:scale "base", :field "at"}
               :y           {:scale "bind-by-index", :field "bidx", :band 0.5}
               :text        {:field "base"}
               :font        {:value "monospace"}
               :fontSize    {:scale "bind", :band 0.4}
               :baseline    {:value "middle"}
               }}}
            ]}))

(defn binds-viz
  [world cell-id]
  (let [cell (get-in world [:cell-pop cell-id])
        dna-vega (basic-dna-viz (:dna cell) (:dna-open? cell) "bind")
        touch-ids (phys/touching (:physics world) cell-id)
        stimuli (sim/find-stimuli world cell-id touch-ids)
        dna (:dna cell)
        open? (:dna-open? cell)
        odna (cell/get-open-dna cell)
        products (keys (:product-counts cell))
        binds (->> (cell/all-binding-sites odna products (map :dna stimuli)
                                           (inc cell/baseline-score))
                   (map-indexed (fn [bind-i bind]
                                  (let [[kind kind-i j] (:path bind)
                                        vs-dna (case kind
                                                 :stimuli (:dna (nth stimuli kind-i))
                                                 :products (nth products kind-i))
                                        ident (case kind
                                               :stimuli (name (:ident (nth stimuli kind-i)))
                                               :products (format "%h" (hash vs-dna)))]
                                    (assoc bind :bind-index bind-i
                                           :vs-dna vs-dna
                                           :name (str ident "@" j))))
                                ))
        binds-data (for [m binds]
                     {:bidx (:bind-index m)
                      :name (:name m)
                      :score (:score m)
                      :begin (dna/offset-into-full-dna (:bind-begin-base m) open?)
                      :end (dna/offset-into-full-dna (inc (:bind-end-base m)) open?)})
        ok #(max 0 (dec %))
        align-data (for [m binds
                         :let [vs-dna (:vs-dna m)
                               amat (ali/alignments odna vs-dna cell/alignment-options)
                               anchor-loc [(inc (:bind-end-base m))
                                           (inc (:vs-end-base m))]]
                         [i vs-i] (ali/match-path anchor-loc amat false)]
                     {:bidx (:bind-index m)
                      :at (dna/offset-into-full-dna (ok i) open?)
                      :base (nth vs-dna (ok vs-i))})
        ]
    (-> (binds-viz* dna-vega binds-data align-data
                    [900 (* 50 (count binds))])
        (with-vega-meta))))

(defn reaction-viz*
  [dna-vega re-data [width-px height-px]]
  (merge-with
   into
   (assoc dna-vega
          :width width-px
          :height height-px
          :padding {:left 20, :right 20, :top 20, :bottom 20})
   {:scales [{:name   "re-step",
              :type   "band"
              :range  "height",
              :domain {:data "re", :field "step-name"
                       :sort {:field "step-index", :order "ascending"}}
              :padding 0.1}
             {:name "read-type-cols"
              :type "ordinal"
              :domain ["operation" "template" "match" "match-full" "silenced" "unsilenced"]
              :range {:scheme "set1"}}]
    :axes [{:scale "re-step"
            :orient "left"}]
    :data [{:name "re"
            :values re-data}]
    :marks [{:type "rect"
             :from {:data "re"}
             :zindex 2
             :encode
             {:enter
              {:x           {:scale "base", :field "begin"}
               :x2          {:scale "base", :field "end"}
               :y           {:scale "re-step", :field "step-name"}
               :height      {:scale "re-step", :band 1.0}
               :fill        {:scale "read-type-cols" :field "read-type"}
               :fillOpacity {:value 0.5}
               :stroke      {:value "black"}}}}
            {:type "text"
             :from {:data "re"}
             :zindex 3
             :encode
             {:enter
              {:x {:scale "base", :field "end"}
               :y {:scale "re-step", :field "step-name", :band 0.5}
               :text {:field "post-energy"}
               :baseline {:value "middle"}
               :align {:value "right"}}}}
            ]
    :legends [{:orient "bottom"
               :fill "read-type-cols"
               :encode
               {:symbols
                {:enter
                 {:shape {:value "square"}
                  :size {:value 400}
                  :fillOpacity {:value 0.5}
                  :stroke {:value "black"}
                  :strokeWidth {:value 1}}}}}]}))

(defn reaction-viz
  [world cell-id time-step]
  (let [cell (get-in world [:cell-pop cell-id])
        dna-vega (basic-dna-viz (:dna cell) (:dna-open? cell) "re-step")
        dna (:dna cell)
        open? (:dna-open? cell)
        touch-ids (phys/touching (:physics world) cell-id)
        re (sim/cell-reaction world cell-id touch-ids time-step)
        re-data (for [[step-i [op-code offset re-step]]
                      (map-indexed vector (:reaction-log re))
                      :let [no (:next-offset re-step)
                            next-off (cond (= no :stop-reaction) nil
                                           (nil? no) (+ offset dna/codon-length)
                                           :else no)
                            read-idx (merge {:operation [offset (+ offset dna/codon-length)]}
                                            (:read-idx re-step))]
                      [read-type [read-begin read-end]] read-idx]
                  {:step-index step-i
                   :step-name (str (name op-code) " ." step-i)
                   :op-code (name op-code)
                   :post-energy (when (= read-type :operation)
                                  (when-let [e (get-in re-step [:cell-immediate :energy])]
                                    (str "e=" e)))
                   ;:next-offset next-off
                   :read-type (name read-type)
                   :begin (cond-> read-begin
                            (#{:template :match} read-type)
                            (dna/offset-into-full-dna open?))
                   :end (cond-> read-end
                          (#{:template :match} read-type)
                          (dna/offset-into-full-dna open?))
                   })]
    (-> (reaction-viz* dna-vega re-data
                       [900 (* 50 (count (:reaction-log re)))])
        (assoc :title {:text (str "Reaction of cell " cell-id
                                  " t=" time-step
                                  " e=" (:energy cell))})
        (with-vega-meta))))

(comment
  (require '[gorilla-repl.core :as gorilla])
  (gorilla/run-gorilla-server {:port 8990})

  ;; http://127.0.0.1:8990/worksheet.html?filename=examples/worksheets/first.clj
  (require '[gorilla-repl.vega :refer [vega-view]])

  (require '[oz.core :as oz])
  (require '[org.nfrac.xylo.viz :as viz])
  (require '[org.nfrac.xylo.sim :as sim])
  )
