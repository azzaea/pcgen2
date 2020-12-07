# Loading dependencies
pacman::p_load("pcalg",
               "Matrix",
               "Hmisc",# needed for using the %nin% operator
               "gaston",# for type A tests
               "Rgraphviz", # for plotting resulting networks
               "sommer", # for calling pcRes
               "clusterGeneration")
source('~/github_repos/pcgen/R/pcgen.R')
source('~/github_repos/pcgen/R/skeleton2.R')
source('~/github_repos/pcgen/R/pcgenTest.R')
source('~/github_repos/pcgen/R/gen.var.test.R')
source('~/github_repos/pcgen/R/res.covar.test.R')
source('~/github_repos/pcgen/R/make.Z.matrix.R')
source('~/github_repos/pcgen/R/fitEM.R')
source('~/github_repos/pcgen/R/construct.block.R')
source('~/github_repos/pcgen/R/deviance2.R')
source('~/github_repos/pcgen/R/pc.cons.intern2.R')
source('~/github_repos/pcgen/R/checkTriple2.R')
source('~/github_repos/pcgen/R/udag2pdagRelaxed2.R')
source("~/github_repos/pcgen/R/pcRes.R")
source("~/github_repos/pcgen/R/getResiduals.R")
 
pc.fit <- pcRes(d, K = K,verbose = T)
plot(pc.fit, main = "Learnt pcRes network")

