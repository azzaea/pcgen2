#!/usr/bin/env dsc
simulate: datamaker.R
  reps: 1                # number of (genetically identical) replicates of each genotype
  randomizeSamples: FALSE # reorder the kinship matrix rows and columns?
  negative: FALSE     # Create a DAG with negative edge weights?
  topsort: TRUE         # If TRUE, non-genetic traits are chosen from nodes with highest topological order
  fixngen: TRUE       # If TRUE, Number of traits with 0 genetic effects = ((1-pgen) * ntraits).
                      # If FALSE, the number is random with probability = 1 - pgen
  useK: TRUE          # use a marker-based kinship matrix
  target: "mouse"     # kinship of mice, from kgp or plant dataset
  simulateQTLs: TRUE  # simulate QTLs, explaining part of the genetic variance ?
  alphaQTL: 0         # The proportion of genetic variance explained by the QTLs
  fixnumber: TRUE     # If TRUE, Number of QTLs = lambdaQTL
                        # If FALSE, number is drawn from Poisson with lambdaQTL
  lambdaQTL: 2        # Number of QTLs or the Poisson parameter for their number
  pQTL: 0.5           # The probability a QTL has an effect on a trait
  noise: F            # Add noise to the simulated graph data?
  extranoise: 0.2     # Additional noise amount
  pgen: 10           # Probability a trait has non-zero genetic variance.
  ## Are these values standard? Or are they arbitrary..
  ntraits: 3         # number of traits
  dagcoeff:   .25, .5, .75, 1   # Lower and upper limits of DAG edge weights 
  rangeGenVar: .01, .1, .5, 1   # Range of variances of direct genetic effects  ---------------> make it low
  rangeEnvVar: .01, .1, .5, 1   # Range of error variances (environment effects) ------------------> make it low
  diffVar: FALSE 
  $true_dag: gsem 
  $data: nw.data 
  $K: K
  $cor_traits: cor.traits

pcres, pcgen, vanilla: pcresbuilder.R, pcgenbuilder.R, vanillapc.R
  data: $data
  K: $K
  gsem: $true_dag
  $learnt_dag: learnt.nw 
  $dagfit_tpr: dag.results['tpr']
  $dagfit_fpr: dag.results['fpr']
  $dagfit_tdr: dag.results['tdr']
  $dagfit_shd: dag.results['shd']
  $genfit_tpr: gen.results['tpr']
  $genfit_fpr: gen.results['fpr']
  $genfit_tdr: gen.results['tdr']
  $cor_residuals: cor.residuals

DSC:
    define:
      learn: pcres, pcgen, vanilla 
    run: simulate * learn 
    exec_path: bin
    output: mouse.simulation.2gtrait.identicalvar
    replicte: 4 
        
