## ------------------------------------------------ ##
## ----          SUMMARISE RESULTS           ------ ##
## ------------------------------------------------ ##
## This script is based on Kostas Angelis' code.    ##
## It was first used to summarise the results       ##
## obtained in Angelis et al., 2017,                ## 
## https://doi.org/10.1093/sysbio/syx061            ##
## ------------------------------------------------ ##

## -------------------------------------- ##
## ----    INITIALIZE ENVIRONMENT   ----- ##
## -------------------------------------- ##
# Clean environment 
rm( list = ls() )

# Load working directory
# NOTE: Change the path according to the directory where you have the analysis running!
wd <- "simulations/effect_x/"


## -------------------------------- ##
## ----    CREATE VARIABLES    ---- ##
## -------------------------------- ##
# Num. replicates
reps <- 1000
# Num. species
sp <- 8


## ------------------------------------ ##
## ----    SUMMARISE REPLICATES    ---- ##
## ------------------------------------ ##

# We carried out the Bayesian inference of divergence times for p = 100, p = 1,000; and p = 10,000 characters.
# A summary file called "results.tsv" is generated from the "out.txt" output by MCMCtree, which just greps 
# the summary statistics about the estimated posterior parameters (divergence times, mean rate, diffusion rate).
# Its format is the following:
#
# Time	     Posterior mean	 95% equal tail CI--low	 95% equal tail CI--up	  95% HPD CI--low	 95% HPD CI--up	  HPD-CI-width
# t_nXX	     A               B	        	            C                       D	               E	              F	
# .
# .
# muYY	     A'	             B'	                      C'	                    D'	             E'	              F'
# .
# .
# sigma2_ZZ  A'	             B'	                      C'	                    D'	             E'	              F'
# .
# .
#
# There are many "t_nXX" as internal nodes, where "XX" is the node number; and as many "muYY" and 
# "sigma2_ZZ" as data partitions (one estimated mean rate and mean sigma^2 per partition)
#
# This output file is then used to obtain the node, the posterior mean, and the lower and upper 95% equal tail CI 
# intervals for the summary statistics.

# p = 100 characters
# Remember to change the path according to the directory where you have the results 
# generated by MCMCtree
rt.p100 <- as.list( numeric( reps ) )
for( i in 1:reps ){
  # Read only entries with divergence times estimates
  # Ommit columns with HPD-CI
  rt.p100[[ i ]] <- read.csv( paste( wd, "p100/replicates/rep_", i,
                                    "/results.tsv", sep = "" ), skip=1, header=F, sep="\t" )[1:(sp-1),-(c(5:8))]
  colnames( rt.p100[[ i ]] ) <- c( "node", "mean", "L", "U" )
  rownames( rt.p100[[ i ]] ) <- rt.p100[[ i ]][,1]
  # Remove first column with the node names as now they are the rownames
  rt.p100[[ i ]] <- rt.p100[[ i ]][,-1]
}

# p = 1,000 characters
# Remember to change the path according to the directory where you have the results 
# generated by MCMCtree
rt.p1000 <- as.list( numeric( reps ) )
for( i in 1:reps ){
  # Read only entries with divergence times estimates
  # Ommit columns with HPD-CI
  rt.p1000[[ i ]] <- read.csv( paste( wd, "p1000/replicates/rep_", i,
                                      "/results.tsv", sep = "" ), skip=1, header=F, sep="\t" )[1:(sp-1),-(c(5:8))]
  colnames( rt.p1000[[ i ]] ) <- c( "node", "mean", "L", "U" )
  rownames( rt.p1000[[ i ]] ) <- rt.p1000[[ i ]][,1]
  # Remove first column with the node names as now they are the rownames
  rt.p1000[[ i ]] <- rt.p1000[[ i ]][,-1]
}

# p = 10,000 characters
# Remember to change the path according to the directory where you have the results 
# generated by MCMCtree
rt.p10000 <- as.list( numeric( reps ) )
for( i in 1:reps ){
   # Read only entries with divergence times estimates
   # Ommit columns with HPD-CI
   rt.p10000[[ i ]] <- read.csv( paste( wd, "p10000/replicates/rep_", i,
                                        "/results.tsv", sep = "" ), skip=1, header=F, sep="\t" )[1:(sp-1),-(c(5:8))]
   colnames( rt.p10000[[ i ]] ) <- c( "node", "mean", "L", "U" )
   rownames( rt.p10000[[ i ]] ) <- rt.p10000[[ i ]][,1]
   # Remove first column with the node names as now they are the rownames
   rt.p10000[[ i ]] <- rt.p10000[[ i ]][,-1]
}


## ------------------------------------ ##
## ************************************ ##