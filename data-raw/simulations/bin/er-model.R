# Loading dependencies
pacman::p_load("pcalg",
#               "Matrix",
#               "Hmisc",# needed for using the %nin% operator
#               "graph",
#               "gaston",# for type A tests
#               "simplePHENOTYPES", # for simulating different phenotypes
#               "microbenchmark", # for timing
#               "Rgraphviz", # for plotting resulting networks
#               "sommer", # for calling pcRes
#               "igraph",
               "clusterGeneration")
#source('~/github_repos/pcgen/R/pcgen.R')
#source('~/github_repos/pcgen/R/skeleton2.R')
#source('~/github_repos/pcgen/R/pcgenTest.R')
#source('~/github_repos/pcgen/R/gen.var.test.R')
#source('~/github_repos/pcgen/R/res.covar.test.R')
#source('~/github_repos/pcgen/R/make.Z.matrix.R')
#source('~/github_repos/pcgen/R/fitEM.R')
#source('~/github_repos/pcgen/R/construct.block.R')
#source('~/github_repos/pcgen/R/deviance2.R')
#source('~/github_repos/pcgen/R/pc.cons.intern2.R')
#source('~/github_repos/pcgen/R/checkTriple2.R')
#source('~/github_repos/pcgen/R/udag2pdagRelaxed2.R')
#source("~/github_repos/pcgen/R/pcRes.R")
#source("~/github_repos/pcgen/R/getResiduals.R")

set.seed(seed)
wd <- "~/github_repos/pcgen_simulations/" 

for (file.name in list.files(paste0(wd,'simulation_functions'), full.names=T)) {source(file.name)}


# General settings ####################################
if (useK) {
  load(file = paste0(wd,'kinship/', Kname))
  K <- K + 0.001 *diag(nrow(K))
} else {
  K <- NULL
}
ng <- nrow(K)    # number of genotypes 
K.pop <- K
if (!is.null(K.pop)) {
  nk <- nrow(K.pop) # number of individuals in the kinship matrix
  stopifnot(nk >= ng)
  genotype.sample <- 1:nk # sample(1:nk, size = ng) # if looking for randomness!
  K <- K.pop[genotype.sample, genotype.sample]
} else {
  K <- NULL
}

p.edge <- 1 * min(1, 1/(p-1)) # probability of an edge, in the randomDAG 
# function Kalisch & Buhlmann (2007) 
# E(N) = p.edge / (p-1) => Average neighbourhood size

######################################################################
 
dag <- randomDAG2(n=p, prob=p.edge, lB = dagcoeff[1], uB = dagcoeff[2],
                  V= paste0('V', 1:p), negative = negative)

plot(dag, main = paste("Simulated Traits only network- Erdoes-Renyi: p.edge =", 
                       p.edge, "weights range = ", dagcoeff[1], dagcoeff[2]))
    
geneticCovMatrix <- genPositiveDefMat("onion", dim = p, alphad=1, 
                                       rangeVar=rangeGenVar)$Sigma
rownames(geneticCovMatrix) <- 
  colnames(geneticCovMatrix) <- paste0('V', 1:p)

if (pgen < 1) {
  if (fixngen) {
    if (topsort == TRUE) {
      traits.without.direct.gen.effect <- 
        rev(match(names(topological.sort(igraph.from.graphNEL(dag))),
                  names(edgeL(dag)))) [1:(round((1 - pgen)*p))]
    } else { # topsort == FALSE
      traits.without.direct.gen.effect <- 
        sample(1:p, size = round((1 - pgen)*p))
    }
  } else { # fixngen == F
    traits.without.direct.gen.effect <-
      which(rbinom(p=1-pgen, size=1, n=p)==1)
  }
} else { # pgen >=1
  traits.without.direct.gen.effect <- integer()
}
geneticCovMatrix[traits.without.direct.gen.effect,] <- 0
geneticCovMatrix[,traits.without.direct.gen.effect] <- 0
errorVars <- runif(min=rangeEnvVar[1], max=rangeEnvVar[2], n=p)
names(errorVars) <- paste0('V', 1:p)
 
gsem <- create.gsem(dag.adj=dag, genCovMatrix=geneticCovMatrix,
                    errorVars = errorVars)
plot(gsem[['gdag']], main = "Simulated  genetic network")

test <- rmvDAG_H2(ng=ng, r=r, gsem=gsem, simulate.QTLs = simulateQTLs,
                  lambda.QTL = lambdaQTL, fix.number = fixnumber,
                  alpha.QTL = alphaQTL, p.QTL = pQTL, 
                  K = K, use.marker.matrix = T, M = M)

d <- test$fulldata[, 1:(p+1)]

if (noise == T) {
  for (j in (2:(p+1))) 
    d[,j] <- d[,j] + rnorm(n = nrow(d),
                           sd = sqrt(extranoise * gsem$Ve[j-1,j-1]))
}


