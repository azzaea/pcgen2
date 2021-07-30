devtools::load_all(path = "/home/azza/github_repos/pcgen2/")
library(pcalg)
library(graph)
library(GraphPhenoSimulator)


cor.residuals <- cor(getResiduals(data, cov.method = "uni", K = K))

learnt.nw <- pcgen(data, K = K, res.cor = cor.residuals, use.res = T, verbose = F, alpha = 1e-2)

gen.results <- compareGeneticEffects(QTL_to_G(learnt.nw@graph), gsem$gdag)

fitted.subgraph <-
  subGraph(snodes = setdiff(nodes(learnt.nw@graph),
                            c('G', names(data)[grep(pattern='qtl_', x=names(data))])),
           learnt.nw@graph)

dag.results <- c(compareGraphs.order(fitted.subgraph, gsem$dag),
                 "shd" = shd.order(fitted.subgraph, gsem$dag))
