
er_model, hub_model, scaleFree_model: er-model.R, hub-model.R, scale-free.R
  seed: 1279
  r: 1            # number of (genetically identical) replicates of each genotype
  useK: TRUE      # use a marker-based kinship matrix
  Kname: "regmap_100snps.RData" #'drops50k.RData' # 'regmap_kinship.RData'
  simulateQTLs: T # simulate QTLs, explaining part of the genetic variance ?
  lambdaQTL: 2    # Poisson parameter for the number of QTLs.
  fixnumber: T    # Fix number of QTLs instead of randomizing by Poisson dist 
  alphaQTL: 0.05  # The proportion of genetic variance explained by the QTLs
  pQTL: 0.5
  p: 3            # number of traits
  rangeGenVar: c(1,2) # Initial uniform limit of variances of the direct genetic effects
  pgen:           # Then, with probability (1-p.gen), they are set to zero
  fixngen: T      # If this is TRUE, a fixed number of round((1-p.gen) * p) is set to 0
  topsort: T      # If TRUE, nodes with highest topological order are set to zero
  rangeEnvVar: c(1,1) # Limits of Environmental variances uniform distribution
  dagcoeff: c(0.5, 1) # Lower and upper limites of DAG edge weights
  negative: FALSE # Create a DAG with negative edge weights?
  noise: T        # Add noise to the simulated graph data?
  extranoise: 0.2 # Additional noise amount
  @hub_model:
    N: 2            # Expected neighbours of a node
    islands: 2      # Islands in the model
  @scaleFree_model:
    N: 2            # Expected neighbours of a node
    alpha: 1        # Pereferntial attachmed parameter (power law)
  $true_dag: gsem 
  $data: d
  $K: K

pcgen, pcgen_res, pcres: pcgen.R, pcgen.res.R, pcres.R
  d: $data
  K: $K
  $learnt_dag: pc.fit

evaluate: eval.nw.R
  gsem: $true_dag
  testdag: $learnt_dag
  data: $data
  $dagfit: dag.results
  $genfit: gen.results

DSC:
  define:
    simulate: scaleFree_model, er_model #hub_model, 
    analyze: pcgen_res, pcres #, pcgen,  
  run: simulate * analyze * evaluate
  exec_path: bin
  output: dsc_benchmarks
  replicate: 10
