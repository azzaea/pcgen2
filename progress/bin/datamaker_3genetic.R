library(GraphPhenoSimulator)

if (target == "mouse") {
  results.dir <- "/home/azza/github_repos/pcgen2/results/mouse"
  K <- as.matrix(data.table::fread(file.path(results.dir, "kinship.cXX.txt"),
                                   sep = "\t", header = F))
} else if (target == "human") {

} else {
  load("~/github_repos/GraphPhenoSimulator/data/simdata1.RData")
  X <- as.matrix(dm[, -c(1:3)])
  K <- X %*% t(X) / sqrt(ncol(X))
  rm(list = c("d", "dm", "X"))
}

ng <- nrow(K)    # number of genotypes


if (randomizeSamples) {
  genotype.sample <-  sample(1:ng, size = ng) # if looking for randomness!
  K <- K[genotype.sample, genotype.sample]
}
p.edge <- 1 * min(1, 1/(ntraits-1)) # probability of an edge, in the randomDAG


dag <- matrix(0, ntraits, ntraits)
dag[1, 2] <- dagcoeff
rownames(dag) <- colnames(dag) <- paste("Trait", 1: ntraits)

if (diffVar == TRUE) {
  geneticCovMatrix <- clusterGeneration::genPositiveDefMat(covMethod = "onion",
                                                           dim=ntraits, alphad=1,
                                                           rangeVar=c(0, rangeGenVar))$Sigma
} else
  geneticCovMatrix <- clusterGeneration::genPositiveDefMat(covMethod = "onion",
                                                           dim=ntraits, alphad=1,
                                                           rangeVar=c(rangeGenVar, rangeGenVar))$Sigma

rownames(geneticCovMatrix) <- colnames(geneticCovMatrix) <- colnames(dag)

traits.without.direct.gen.effect <- rownames(dag)[1: 2]

# if (pgen < 1) {
#   if (fixngen) {
#     if (topsort == TRUE) {
#       traits.without.direct.gen.effect <-
#         rev(match(names(topological.sort(igraph.from.graphNEL(dag))),
#                   names(edgeL(dag)))) [1:(round((1 - pgen)*ntraits))]
#     } else { # topsort == FALSE
#       traits.without.direct.gen.effect <-
#         sample(1:ntraits, size = round((1 - pgen)*ntraits))
#     }
#   } else { # fixngen == F
#     traits.without.direct.gen.effect <-
#       which(rbinom(p=1-pgen, size=1, n=ntraits)==1)
#   }
# } else { # pgen >=1
#   traits.without.direct.gen.effect <- integer()
# }

geneticCovMatrix[traits.without.direct.gen.effect, ] <- 0
geneticCovMatrix[ ,traits.without.direct.gen.effect] <- 0

if (diffVar == TRUE) {
  errorVars <- runif(min = 0, max = rangeEnvVar, n = ntraits)
} else
  errorVars <- runif(min = rangeEnvVar, max = rangeEnvVar, n = ntraits)

names(errorVars) <- names(dag)

gsem <- create.gsem(dag.adj = dag, genCovMatrix = geneticCovMatrix,
                    errorVars = errorVars)


test <- rmvDAG_H2(ng = ng, r = reps, gsem = gsem, simulate.QTLs = simulateQTLs,
                  lambda.QTL = lambdaQTL, fix.number = fixnumber,
                  alpha.QTL = alphaQTL, p.QTL = pQTL,
                  K = K, use.K.names = FALSE, use.marker.matrix = T)

nw.data <- test$fulldata[, 1:(ntraits+1)]

cor.traits <- cor(nw.data[, -1])["Trait.1", "Trait.2"]
