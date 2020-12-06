# Loading dependencies
pacman::p_load("pcalg",
               "Matrix",
               "Hmisc",# needed for using the %nin% operator
               "graph",
               "gaston",# for type A tests
               "simplePHENOTYPES", # for simulating different phenotypes
               "microbenchmark", # for timing
               "Rgraphviz", # for plotting resulting networks
               "sommer", # for calling pcRes
               "igraph",
               "clusterGeneration"
)

rm(list=ls()); gc()

source('~/github_repos/pcgen/R/pcgen.R')
source('~/github_repos/pcgen/R/skeleton2.R')
source('~/github_repos/pcgen/R/pcgenTest.R')
source('~/github_repos/pcgen/R/gen.var.test.R')
source('~/github_repos/pcgen/R/res.covar.test.R')
source('~/github_repos/pcgen/R/make.Z.matrix.R')
source('~/github_repos/pcgen/R/fitEM.R')

# The codes below have not been altered- They are called by the above functions though:
source('~/github_repos/pcgen/R/construct.block.R')
source('~/github_repos/pcgen/R/deviance2.R')
source('~/github_repos/pcgen/R/pc.cons.intern2.R')
source('~/github_repos/pcgen/R/checkTriple2.R')
source('~/github_repos/pcgen/R/udag2pdagRelaxed2.R')

# These functions only learn the network between traits- not modified
source("~/github_repos/pcgen/R/pcRes.R")
source("~/github_repos/pcgen/R/getResiduals.R")

set.seed(1279)
wd <- "~/github_repos/pcgen_simulations/" 
result.file <- 'sim_Azza.RData'

setwd(paste0(wd,'simulation_functions'))
for (file.name in list.files()) {source(file.name)}
setwd(wd)


# General settings ####################################
r  <- 1      # number of (genetically identical) replicates of each genotype
extra.noise <- 0.2
n.sim <- 2  # number of simulations
use.K <- TRUE # use a marker-based kinship matrix
K.name <- 'regmap_100snps.RData' #'drops50k.RData' # 'regmap_kinship.RData'
if (use.K) {
  load(file = paste0(wd,'kinship/', K.name))
  K <- K + 0.001 *diag(nrow(K))
} else {
  K <- NULL
}
ng <- nrow(K)    # number of genotypes 
nTr <- 700    # number of samples in each simulation- moot now
K.pop <- K
method.names <- c('pcgen', 'pcgen_with_res', 'pcRes')
method.numbers <- 1 # chosen method numbers for evaluation
n.method <- length(method.names)
######################### QTLs

simulate.QTLs <- T # simulate QTLs, explaining part of the genetic variance ?
lambda.QTL <- 2 # Poisson parameter. For each simulated data set, the total 
# number of QTLs is drawn from a Poisson distribution
# (and put to 1, in case of a zero outcome)
# The QTLs scores (0,1) are simulated as well, and always have
# frequency 0.5
fix.number <- T # When fix.number is TRUE, we don't draw from a Poisson
# distribution, but instead fix the number of QTLs to
# round(lambda.QTL)
alpha.QTL <- 0.05   # The proportion of genetic variance explained by the QTLs
# together. More precisely, for each trait with direct genetic
# effects, a proportion of (1-alpha.QTL) of the genetic variance
# is still explained by a random genetic effect, while a
# proportion of alpha.QTL is now explained by some of the QTLs.
# Given that there are n.QTL QTLs in total (for all the traits),
# we randomly decide if a given QTL affects a given trait,
# with Bernoulli trials (p.QTL); minimum one QTL per trait
p.QTL <- 0.5

###############################

run.example <- F  # run a specific example (true DAG) ?
# Otherwise, generate a DAG each time, using the parameters below
# Choose an example from the /examples folder, if run.example==T
example.file <- 'rep1.R'

######################
# random GDAG generation, if run.example==F
p  <- 5    # number of traits
p.edge <- 1 * min(1, 1/(p-1)) # probability of an edge, in the randomDAG 
# function Kalisch & Buhlmann (2007) 
# E(N) = p.edge / (p-1) => Average neighbourhood size

rangeGenVar <- c(1,2) # Initially, variances of the direct genetic effects are 
# drawn from a uniform distribution on this support
p.gen       <- 0.5   # Then, with probability (1-p.gen), they are set to zero
fix.n.gen <- T # If this is TRUE, a fixed number of round((1-p.gen) * p) is set to 0
top.sort <- T  # If TRUE, nodes with highest topological order are set to zero
rangeEnvVar <- c(1,1) # Environmental variances are drawn from a uniform 
# distribution on this support

dag.coeff   <- c(0.5, 1) # Lower and upper limites of DAG edge weights
negative    <- FALSE # Create a DAG with negative edge weights?

# Parameter used in the genPositiveDefMat function, used to randomly generate
# the genetic covariance matrix.
# For a diagonal matrix, put Inf
eta <- 1 # for uniform

######################################################################

gsem.list <- list()
sim.result.list <- list()
for (i in 1:n.method) sim.result.list[[i]] <- list()
names(sim.result.list) <- method.names
sim.data.list.without.noise <- sim.data.list <- list()
(start.date <- date())

for (sim in 1:n.sim) {
  #sim=1
  if (!is.null(K.pop)) {
    nk <- nrow(K.pop) # number of individuals in the kinship matrix
    stopifnot(nk >= ng)
    genotype.sample <- 1:nk # sample(1:nk, size = ng) # if looking for randomness!
    K <- K.pop[genotype.sample, genotype.sample]
    rownames(K) <- colnames(K) <- paste0('g',1:ng)
  } else {
    K <- NULL
  }
  
  cat('########','\n','# Simulation: ',sim,'\n','#############','\n')
  if (!run.example) {
    dag <- randomDAG2(n=p, prob=p.edge, lB = dag.coeff[1], uB = dag.coeff[2],
                      V= paste0('V', 1:p), negative = negative)
    #names(topological.sort(igraph.from.graphNEL(dag)))
    #edgeWeights(dag)
    plot(dag, main = paste("Simulated Traits only network- Erdoes-Renyi: p.edge =", 
                           p.edge, "weights range = ", dag.coeff[1], dag.coeff[2]))
    
    geneticCovMatrix <- genPositiveDefMat("onion", dim = p, alphad=1, 
                                          eta=eta, rangeVar=rangeGenVar)$Sigma
    rownames(geneticCovMatrix) <- 
      colnames(geneticCovMatrix) <- paste0('V', 1:p)
    if (p.gen < 1) {
      if (fix.n.gen) {
        if (top.sort == TRUE) {
          traits.without.direct.gen.effect <- 
            rev(match(names(topological.sort(igraph.from.graphNEL(dag))),
                      names(edgeL(dag)))) [1:(round((1 - p.gen)*p))]
        } else { # top.sort == FALSE
          traits.without.direct.gen.effect <- 
            sample(1:p, size = round((1 - p.gen)*p))
        }
      } else { # fix.n.gen == F
        traits.without.direct.gen.effect <-
          which(rbinom(p=1-p.gen, size=1, n=p)==1)
      }
    } else { # p.gen >=1
      traits.without.direct.gen.effect <- integer()
    }
    geneticCovMatrix[traits.without.direct.gen.effect,] <- 0
    geneticCovMatrix[,traits.without.direct.gen.effect] <- 0
    errorVars <- runif(min=rangeEnvVar[1], max=rangeEnvVar[2], n=p)
    names(errorVars) <- paste0('V', 1:p)
  } else { # RunExample==T
    source(paste0('examples/',example.file))
    p <- length(errorVars)
  }
  
  gsem <- create.gsem(dag.adj=dag, genCovMatrix=geneticCovMatrix,
                      errorVars = errorVars)
  plot(gsem[['gdag']], main = "Simulated  genetic network")
  gsem.list[[sim]] <- gsem
  test <- rmvDAG_H2(ng=ng, r=r, gsem=gsem, simulate.QTLs = simulate.QTLs,
                    lambda.QTL = lambda.QTL, fix.number = fix.number,
                    alpha.QTL = alpha.QTL, p.QTL = p.QTL, 
                    K = K, use.marker.matrix = T, M = M)
  
  d <- test$fulldata[, 1:(p+1)]
  ####### ADD measurement error
  d_without_noise <- d
  for (j in (2:(p+1))) {
    d[,j] <- d[,j] + rnorm(n = nrow(d), 
                           sd = sqrt(extra.noise * gsem$Ve[j-1,j-1]))
  }
  tr.set <- rownames(K) #sample(rownames(K), size = nTr) # for randomness
  d <- d[(d$G %in% tr.set), ]      #c(pc.set, tr.set)
  d_without_noise <- d_without_noise[(d_without_noise$G %in% tr.set), ]   
  
  K <- K[d$G, d$G]
  d$G <- factor(as.character(d$G))
  d_without_noise$G <- factor(as.character(d$G))
  sim.data.list[[sim]] <- d
  sim.data.list.without.noise[[sim]] <- d_without_noise
  ####################################
  
  if (1 %in% method.numbers){ # pcgen
    pc.fit1 <- pcgen(d_without_noise, K = K, verbose = T, NAdelete = F)
    plot(pc.fit1, main = "Learnt pcgen network")
    sim.result.list[[1]][[sim]] <- pc.fit1
  }
  
  if (2 %in% method.numbers){ # pcgen_with_res
    cor <- cor(getResiduals(suffStat = d_without_noise, K = K))
    pc.fit2 <- pcgen(d_without_noise, use.res = T, K = K, res.cor = cor) 
    plot(pc.fit2, main = "Learnt pcgen and residuals network")
    sim.result.list[[2]][[sim]] <- pc.fit2
  }
  
  if (3 %in% method.numbers){ # pcRes
    pc.fit3 <- pcRes(d_without_noise, K = K,verbose = T)
    plot(pc.fit3, main = "Learnt pcRes network")
    sim.result.list[[3]][[sim]] <- pc.fit3
  }
  
} # end loop across simulation

#iplotPCgen(pc.fit2)
# simdata <- d; save(simdata, file = 'simdata.RData')
# save.image(file='F:/debugPCgen.RData')

(end.date <- date())

######################################################
# assess results

dag.results <- list()

gen.results <- list()

# the results for the edges among the traits
for (i in 1:n.method) {dag.results[[i]] <- matrix(0,n.sim,4)}

# the results for edges G --> traits
for (i in 1:n.method) {gen.results[[i]] <- matrix(0,n.sim,3)}

for (sim in 1:n.sim) { # sim=1
  for (mmm in intersect(c(1:2), method.numbers)) { # mmm=2
    fitted.subgraph <- subGraph(snodes = setdiff(
      nodes(sim.result.list[[mmm]][[sim]]@graph),
      c('G', names(sim.data.list[[sim]])[grep(pattern='qtl_', 
                                              x=names(sim.data.list[[sim]]))])),
      sim.result.list[[mmm]][[sim]]@graph)
    
    fitted.graph         <- sim.result.list[[mmm]][[sim]]@graph
    
    gen.results[[mmm]][sim,]    <- compareGeneticEffects(QTL_to_G(fitted.graph), gsem.list[[sim]]$gdag)
    
    dag.results[[mmm]][sim,1:3] <- compareGraphs(fitted.subgraph, gsem.list[[sim]]$dag)
    dag.results[[mmm]][sim,4]   <- shd(fitted.subgraph, gsem.list[[sim]]$dag)
  }
  # Comparing traits-only networks  
  for (mmm in intersect(c(3), method.numbers)) { # mmm = 3
    dag.results[[mmm]][sim,1:3] <- compareGraphs(sim.result.list[[mmm]][[sim]]@graph, gsem.list[[sim]]$dag)
    dag.results[[mmm]][sim,4]   <- shd(sim.result.list[[mmm]][[sim]]@graph, gsem.list[[sim]]$dag)
  }
}

for (mmm in 1:n.method) {
  colnames(dag.results[[mmm]]) <- c('tpr', 'fpr', 'tdr', 'shd')
  rownames(dag.results[[mmm]]) <- paste0('sim', 1:n.sim)
  
  colnames(gen.results[[mmm]]) <- c('tpr', 'fpr', 'tdr')
  rownames(gen.results[[mmm]]) <- paste0('sim', 1:n.sim)
}
###

dag.results.summary <- matrix(NA, n.method, 4)
rownames(dag.results.summary) <- method.names
colnames(dag.results.summary) <- colnames(dag.results[[1]])

gen.results.summary <- matrix(NA, n.method, 3)
rownames(gen.results.summary) <- method.names
colnames(gen.results.summary) <- colnames(gen.results[[1]])

for (i in 1:n.method) {
  dag.results.summary[i, ] <- apply(dag.results[[i]],2,mean)
  gen.results.summary[i, ] <- apply(gen.results[[i]],2,mean)
}

dag.results.summary
gen.results.summary
setwd(paste0(wd,'results'))
save.image(file = result.file)
