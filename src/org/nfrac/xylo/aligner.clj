(ns org.nfrac.xylo.aligner
  "copied from https://github.com/ship561/str_alignment/blob/master/src/str_alignment/aligner.clj

   Changed align function to return offset positions of best match.")

(def ^{:dynamic true :private true} global)
(def ^{:dynamic true :private true} match)
(def ^{:dynamic true :private true} mis)
(def ^{:dynamic true :private true} gap-open)
(def ^{:dynamic true :private true} gap-ext)
(def ^{:dynamic true :private true} *anchor-left*)
(def ^{:dynamic true :private true} *anchor-right*)

(defn- array-keys
  "positions of the array to work on"

  [s1 s2]
  (for [i (range (count s1)) ;initialize scoring array. similar to a sparse matrix
        j (range (count s2))]
    [i j]))

(defn- init-array
  "initial conditions of the array to fill in for the gapped row/col"

  [s1 s2 locations]
  (->> locations
       (filter #(or (zero? (first %))
                    (zero? (second %))))
       (reduce (fn [m [i j]]
                 (assoc m [i j]
                        (cond
                         (= i j 0) [0 :stop \- \-]
                         (= i 0) [(if global (* gap-open j) 0) :l \- (.charAt s2 j)]
                         (= j 0) [(if global (* gap-open i) 0) :u (.charAt s1 i) \-])))
               {})))

(defn- dir->coord [dir i j]
  (case dir
    :d [(dec i) (dec j)]
    :u [(dec i) j]
    :l [i (dec j)]
    :stop))

(defn- get-score [m dir i j] (->> (dir->coord dir i j) m))

(defn- maxa
  "Determine max score and direction for position [i j]"
  [coll]
  (->> coll
       (filter #(= (apply max (map second coll)) (second %)))
       (sort-by first)
       first))


(defn- gapfn
  "For gap extention penalties, the direction of the previous cell must also be checked"

  [recur-from]
  (case recur-from
    :d gap-open
    :u gap-ext
    :l gap-ext))

(defn- fill-array

  [locations s1 s2]
  (reduce (fn [m [i j]];;score array format
            (let [[d dfrom] (get-score m :d i j) ;;score match/mismatch (diagonal)
                  [u ufrom] (get-score m :u i j) ;;score deletion (above)
                  [l lfrom] (get-score m :l i j) ;;score insertion (left)
                  aa1 (.charAt s1 i) ;;current char in s1
                  aa2 (.charAt s2 j) ;;current char in s2
                  ;;chooses from d, u, l and scores associated with it.
                  [from score] (maxa [(if (= aa1 aa2 )
                                        [:d (+ d match)]
                                        [:d (+ d mis)])
                                      [:u (+ u (gapfn ufrom))]
                                      [:l (+ l (gapfn lfrom))]])]
              (assoc m [i j] ;;insertion of the best score into the matrix
                     (case from
                       :d [score :d aa1 aa2]
                       :u [score :u aa1 \-]
                       :l [score :l \- aa2]))))
          (init-array s1 s2 locations)
          (remove #(or (zero? (first %))
                       (zero? (second %))) locations)))

(defn- trace-back
  "Traceback to to the beginning of the matrix from a starting
  location. If anchor-left is true then the traceback continues to the
  top-left corner."

  [score start-loc H]
  (loop [loc start-loc
         aln_s1 (list)
         aln_s2 (list)]
    (let [[cscore dir a1 a2] (get H loc) ;;stores the next location [score[i j] to go to in H]
          next-coord (apply dir->coord dir loc)]
      (if (or (pos? cscore)
              (and *anchor-left* (not= :stop next-coord)))
        (recur next-coord
               (cons a1 aln_s1) ;;builds strings up from the right to left
               (cons a2 aln_s2))
        (if (= \- a1 a2)
          [score (apply str aln_s1) (apply str aln_s2)]
          [score (apply str a1 aln_s1) (apply str a2 aln_s2)])))))

(defn align
  "Returns [i1 i2 score]
  where
  i1 is index of end of match on seq1.
  i2 is index of end of match on seq2.
  score is the SWa match score.
  If there are multiple equally good alignments, returns first."
  [seq1 seq2 & {:keys [global
                       anchor-right
                       anchor-left
                       match-weight
                       mismatch-weight
                       gap-open-weight
                       gap-ext-weight]
                :or {global false
                     anchor-right false
                     anchor-left false
                     match-weight 2
                     mismatch-weight -1
                     gap-open-weight -1
                     gap-ext-weight -1}}]
  (binding [global global
            match match-weight  ;match
            mis mismatch-weight ;mismatch
            gap-open gap-open-weight
            gap-ext gap-ext-weight
            *anchor-left* (if global true anchor-left)
            *anchor-right* anchor-right]
    (let [s1 (str "-" seq1)
          s2 (str "-" seq2)
          H (fill-array (array-keys s1 s2) s1 s2) ;creates score matrix
          start (cond global
                      (get H (mapv dec [(count s1) (count s2)]));bottom right corner
                      *anchor-right*
                      (->> (filter #(= (ffirst %) (dec (count s1))) H)
                           (sort-by #(-> % second first) >)
                           first
                           second)
                      :else
                      (-> (sort-by #(-> % second first) H) ;finds highest value in matrix
                          last
                          second))
          start-locs (map first (filter #(= start (val %)) H)) ;;starts traceback from this highest value
          [i1 i2] (first start-locs)]
      [i1 i2 (first start)]
      ;;(map (bound-fn [loc] (trace-back (first start) loc H)) start-locs);thread-bindings
                                        ;(mapv #(trace-back (first start) % H) start-locs))));implicit doall. returns vector
      )))
