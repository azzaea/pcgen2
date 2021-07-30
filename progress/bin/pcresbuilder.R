devtools::load_all(path = "/home/azza/github_repos/pcgen2/")
library(pcalg)

learnt.nw <- pcRes(data, K = K,
                    cov.method = "uni", m.max = Inf, verbose = F, alpha = 1e-2)

dag.results <- c(compareGraphs.order(learnt.nw@graph, gsem$dag),
                 "shd" = shd.order(learnt.nw@graph, gsem$dag))

gen.results <- rep(NA, 3)
names(gen.results) <- c('tpr', 'fpr', 'tdr')

cor.residuals <- NA
