devtools::load_all(path = "/home/azza/github_repos/pcgen2/")
library(pcalg)

learnt.nw <- pcalg::pc(suffStat = list(C = cor(data[,-1]), n = nrow(data)),
                       indepTest = pcalg::gaussCItest, alpha = 0.01,
                       labels = names(data[,-1]), m.max = Inf, # for pcRes
                       skel.method = "stable", u2pd = "relaxed", NAdelete = FALSE,
                       conservative = FALSE, maj.rule = TRUE, solve.confl = TRUE)

dag.results <- c(compareGraphs.order(learnt.nw@graph, gsem$dag),
                 "shd" = shd.order(learnt.nw@graph, gsem$dag))

gen.results <- rep(NA, 3)
names(gen.results) <- c('tpr', 'fpr', 'tdr')

cor.residuals <- NA

