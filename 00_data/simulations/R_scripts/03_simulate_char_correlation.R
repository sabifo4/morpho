#------------------------------------------------#
# EFFECT OF WITHIN-LINEAGE CHARACTER CORRELATION #
#------------------------------------------------#

# Before running this code, you should have in your
# working directory the following architecture or something similar.
#
# The data are simulated to have a fixed population noise (c=0.25)
# under the constant correlation model with a specific correlation value "rho",
# that is, with equal correlations "rho". We test different values
# of rho=0,0.25,0.35,0.50,0.70,0.80,0.90.
#
# Note that scaling was performed using the estimate of the population noise 
# calculated on the simulated sampled population so we could simulate 
# what we need to do when real data are used.
# We only scale by the true population noise the data under
# the directory `MstrueR.Rtrue`.
# 
# - your_wd
#     |-alignments
#         |- rho<X>                                                                           # correct c?           correct R?
#              |- MnR            # not scaled | not transformed                           -->   N | c=0 (ignored)      N | R=0  (ignored)
#              |   |-p100
#              |   |-p1000
#              |   |-p10000
#              |   
#              |- MsR            # scaled     | not transformed                           -->   Y | c=1 (est.c)        N | R=0  (ignored)
#              |   |-p100
#              |   |-p1000
#              |   |-p10000
#              |   
#              |- MsR.R0.01      # scaled     | transformed, shrunk R with fixed d = 0.01 -->   Y | c=1 (est.c)        Y | R!=0 (corrected)
#              |   |-p100
#              |   |-p1000
#              |   |-p10000
#              |   
#              |- MsR.Rdef       # scaled     | transformed, shrunk R with est. d = def   -->   Y | c=1 (est.c)        Y | R!=0 (corrected)
#              |   |-p100
#              |   |-p1000
#              |   |-p10000
#              |   
#              |- MstrueR.Rtrue  # scaled     | transformed, true R                       -->   Y | c=1 (true.c)       Y | R!=0 (corrected)
#                  |-p100
#                  |-p1000
#                  |-p10000
#                  
#            
# If you want to change your file architecture, change
# the path to the output file in the "mmcmc3::write.morpho"
# function in the "filename" option. This is used to save
# the alignments following the architecture
# described above. 
#
# This script should be run as it follows: 
#   $ RScript 03_simulate_char_correlation.R <path_to_wd> <num_chars> <rho_val>
#
# E.g.: Simulate alignments with p = 100 characters and rho = 0.50
#   $ RScript 03_simulate_char_correlation.R /home/morpho/p100/ 100 0.50

# Clean environment
rm( list = ls( ) )

# Load libraries
library( mcmc3r )
library( ape )

# Collect args 
args  <- commandArgs(TRUE)

# Set working directory to this file location
# It uses the first arg
wd <- args[1]
setwd( wd )
cat( "===========================\n")
cat( " SETTING WORKING DIRECTORY\n")
cat( "===========================\n")
cat( "Your working directory is ", wd, "\n\n" )

# Set variables for number of characters. Takes argument 2
pvar   <- paste( "p", args[2], sep = "" )
p      <- pvar2 <- as.numeric( args[2] )
cat( "=============================\n")
cat( " SETTING GLOBAL VARIABLES...\n")
cat( "=============================\n")
cat( "A. You are simulating matrices with p = ", p, " characters\n\n" )

# Set number of replicates 
reps <- 1000

# Set number of species 
sp <- 8

# Set variance
c <- 0.25

# Set rho value to evaluate low and high
# correlations among characters
# This is argument 3
cat( "B. The size of your correlation matrix is ", p, "x", p, "\n\n" )
rho     <- as.numeric( args[3] )
rho.wd  <- paste( "rho", args[3], "/", sep = "" )
R       <- matrix( rho, p, p )
diag(R) <- 1

# Set number of specimens to sample for 
# the simulated population
psample = 20

# Get tree 
tree <- ape::read.tree(text="((((A^1:0.1,B^1:0.1):0.2,(F^0.9:0.1,C^1:0.2):0.1):0.5,H^0.3:0.1):0.2,(D^1:0.7,(G^0.7:0.2,E^1:0.5):0.2):0.3);")


# Write each alignment in one file
for ( i in seq( 1:reps ) ){
  
  cat( "=================================================\n")
  cat( " GENERATING ALIGNMENTS FOR REPLICATE ", i, "...\n")
  cat( "=================================================\n\n")
  
  # The seed number changes for every alignment using values from 1 to "R"
  # E.g.:
  # For alignment "1", seed is 1231 (123[1]),
  # for alignment "2", it is 1232 (123[2]), and so on
  seed  <- paste( 123, i, sep = "" )
  cat( "The seed used for replicate ", i, "is ", seed, "\n\n" )
  set.seed( seed )

  # Simulate  n = 100 continuous characters with within population variance c = 0.25
  # and correlation rho = 0.50.
  sim.cont.chars    <- mcmc3r::sim.morpho( tree = tree, n = p, c = c, R = R )
  
  # Simulate a population sample with psample = 20 specimens, within population 
  # variance c = 0.25 and correlation rho = 0.50, and n = 100 characters
  cat( "\nSimulating populatin matrix...\nEstimating R.sh...\n" )
  simPop            <- mcmc3r::sim.pop( psample = psample, n = p, c = c, R = R ) 
  
  # Write the alignment that has not been scaled nor transformed, i.e.,
  # ignores population noise and character correlation.
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'alignments/', rho.wd, 'MnR/', pvar, '/rtraitcont', pvar2, '_MnR_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
  
  # Write alignment after scaling by estimated population noise and matrix transformation
  # is carried out with the shrinkage correlation matrix estimated with delta = default mode
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'alignments/', rho.wd, 'MsR.Rdef/', pvar, '/rtraitcont', pvar2, '_MsR.Rdef_',
                                          i,'.txt', sep='' ),
                        c = simPop$var, R = cbind(simPop$Rsh), method = "chol", names = rownames( sim.cont.chars ) )
  
  # Write alignment after scaling by estimated population noise and matrix transformation
  # is carried out with the shrinkage correlation matrix estimated with delta = 0.01 (fixed value)
  cat( "\nEstimating R.sh.0.01...\n" )
  Rsh.0.01 <- as.matrix( corpcor::cor.shrink( simPop$P, lambda = 0.01 ) )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'alignments/', rho.wd, 'MsR.R0.01/', pvar, '/rtraitcont', pvar2, '_MsR.R0.01_',
                                          i,'.txt', sep='' ),
                        c = simPop$var, R = cbind(Rsh.0.01), method = "chol", names = rownames( sim.cont.chars ) )
  
  # Write alignment after scaling by true population noise and matrix transformation is 
  # carried out using the true correlation matrix
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'alignments/', rho.wd, 'MstrueR.Rtrue/', pvar, '/rtraitcont', pvar2, '_MstrueR.Rtrue_',
                                          i,'.txt', sep='' ),
                        c = c, R = R, method = "chol", names = rownames( sim.cont.chars ) )
  
  cat( "All alignments for replicate ", i, "have been saved!\n\n*\n*\n*\n*\n\n")
}
