library(GraphPhenoSimulator)

if (target == "mouse") {
  results.dir <- "/home/azza/github_repos/pcgen2/results/mouse"
  K <- as.matrix(data.table::fread(file.path(results.dir, "kinship.cXX.txt"),
                                   sep = "\t", header = F) )
  use.K.names <- FALSE
} else if (target == "human") {

} else {
  load("~/github_repos/GraphPhenoSimulator/data/simdata1.RData")
  X <- as.matrix(dm[, -c(1:3)])
  K <- X %*% t(X) / sqrt(ncol(X))
  rm(list = c("d", "dm", "X"))
  use.K.names <- TRUE
}

K <- K + 0.005 *diag(nrow(K))
ng <- nrow(K)    # number of genotypes

if (randomizeSamples) {
  genotype.sample <-  sample(1:ng, size = ng) # if looking for randomness!
  K <- K[genotype.sample, genotype.sample]
}
p.edge <- 1 * min(1, 1/(ntraits-1)) # probability of an edge, in the randomDAG


# dag <- matrix(0, ntraits, ntraits)
# dag[1, 2] <- dagcoeff
# rownames(dag) <- colnames(dag) <- paste("Trait", 1: ntraits)

dag <- randomDAG2(n=ntraits, prob=p.edge, lB = dagcoeff[1], uB = dagcoeff[2],
                  V= paste0('Trait', 1:ntraits), negative = negative)


if (diffVar == TRUE) {
  geneticCovMatrix <- clusterGeneration::genPositiveDefMat(covMethod = "onion",
                                                           dim=ntraits, alphad=1,
                                                           rangeVar=c(0, rangeGenVar))$Sigma
} else
  geneticCovMatrix <- clusterGeneration::genPositiveDefMat(covMethod = "onion",
                                                           dim=ntraits, alphad=1,
                                                           rangeVar=c(rangeGenVar[1], rangeGenVar[2]))$Sigma

rownames(geneticCovMatrix) <- colnames(geneticCovMatrix) <- graph::nodes(dag)

if (pgen < 1) {
  if (fixngen) {
    if (topsort == TRUE) {
      traits.without.direct.gen.effect <-
        rev(match(names(igraph::topological.sort(igraph::igraph.from.graphNEL(dag))),
                  names(graph::edgeL(dag)))) [1:(round((1 - pgen)*ntraits))]
    } else { # topsort == FALSE
      traits.without.direct.gen.effect <-
        sample(1:ntraits, size = round((1 - pgen)*ntraits))
    }
  } else { # fixngen == F
    traits.without.direct.gen.effect <-
      which(rbinom(p=1-pgen, size=1, n=ntraits)==1)
  }
} else { # pgen >=1
  traits.without.direct.gen.effect <- integer()
}

geneticCovMatrix[traits.without.direct.gen.effect, ] <- 0
geneticCovMatrix[ ,traits.without.direct.gen.effect] <- 0

if (diffVar == TRUE) {
  errorVars <- runif(min = 0, max = rangeEnvVar, n = ntraits)
} else
  errorVars <- runif(min = rangeEnvVar[1], max = rangeEnvVar[2], n = ntraits)

names(errorVars) <- graph::nodes(dag)

gsem <- create.gsem(dag.adj = dag, genCovMatrix = geneticCovMatrix,
                    errorVars = errorVars)

test <- rmvDAG_H2(ng = ng, r = reps, gsem = gsem, simulate.QTLs = simulateQTLs,
                  lambda.QTL = lambdaQTL, fix.number = fixnumber,
                  alpha.QTL = alphaQTL, p.QTL = pQTL,
                  K = K, use.K.names = use.K.names, use.marker.matrix = T)

nw.data <- test$fulldata[, 1:(ntraits+1)]

if (noise == T) {
  for (j in (2:(ntraits+1)))
    nw.data[,j] <- nw.data[,j] + rnorm(n = nrow(nw.data),
                           sd = sqrt(extranoise * gsem$Ve[j-1,j-1]))
}


cor.traits <- cor(nw.data[, -1])


