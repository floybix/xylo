(defproject org.nfrac/xylo "0.1.0-SNAPSHOT"
  :description "Artificial life simulation of plant-like cells"
  :url "https://github.com/floybix/xylo"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.9.0"]
                 [org.clojure/core.async "0.4.474"]
                 [org.clojure/test.check "0.9.0"]
                 [aysylu/loom "1.0.1"]
                 ]
  :resource-paths ["lib/parasail.jar"]
  :jvm-opts ["-Djava.library.path=lib/"]

  :profiles {:dev {:dependencies [[gorilla-repl "0.4.1-SNAPSHOT"]
                                  [com.gfredericks/test.chuck "0.2.9"]]}}
)
