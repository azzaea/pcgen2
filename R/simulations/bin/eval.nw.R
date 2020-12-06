# Loading dependencies
pacman::p_load("pcalg",
               "Matrix",
               "Hmisc",# needed for using the %nin% operator
               "gaston",# for type A tests
               "Rgraphviz", # for plotting resulting networks
               "sommer", # for calling pcRes
               "clusterGeneration")
wd <- "~/github_repos/pcgen_simulations/"
for (file.name in list.files(paste0(wd,'simulation_functions'), full.names=T)) {source(file.name)}

if ("G" %in% nodes(testdag@graph)) {
  fitted.subgraph <- subGraph(snodes = setdiff(
                                          nodes(testdag@graph),
                                          c('G', names(data)[grep(pattern='qtl_', x=names(data))])),
                                testdag@graph)
    gen.results <- compareGeneticEffects(QTL_to_G(testdag@graph), gsem$gdag)
} else {
  fitted.subgraph <- testdag@graph
  gen.results <- NULL
}

dag.results <- vector("numeric", 4)
dag.results[1:3] <- compareGraphs(fitted.subgraph, gsem$dag)
dag.results[4] <- shd(fitted.subgraph, gsem$dag)
names(dag.results) <- c('tpr', 'fpr', 'tdr', 'shd')

