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
# This script is thought to be run as a job array in an HPC.
# You should write the following within the bash script you will 
# submit to your cluster: 
#   $ RScript 03_simulate_char_correlation_cluster.R <path_to_wd> <task_ID> \
#     <num_chars> <rho_val>
#
# E.g.: Simulate alignment 1 with p = 100 characters and rho = 0.50. 
#       The variable we use as a tag ID in the Apocrita HPC is $SGE_TASK_ID, 
#       but please modify according to the parameters set in your cluster 
#       called $SGE_TASK_ID
#   $ RScript 03_simulate_char_correlation_cluster.R /home/scratch/usr12345/morpho/ \
#     $SGE_TASK_ID 100 0.50

# Clean environment
rm( list = ls( ) )

# Load libraries
library( mcmc3r )
library( ape )

# Collect args 
args  <- commandArgs(TRUE)

cat( "====================================\n")
cat( " START JOB FOR REPLICATE NUMBER", args[2], "\n")
cat( "====================================\n\n")


# Set working directory to this file location
# It uses the first arg
wd <- args[1]
setwd( wd )
cat( "===========================\n")
cat( " SETTING WORKING DIRECTORY\n")
cat( "===========================\n")
cat( "Your working directory is ", wd, "\n\n" )

# This script is thought to run as a job array in an HPC.
# The seed number changes for every alignment as values from 1 to "R"
# are passed as an argument.
# E.g.:
# For alignment "1", seed is 1231 (123[1]),
# for alignment "2", it is 1232 (123[2]), and so on
seed  <- paste( 123, args[2], sep = "" ) 
cat( "=================\n")
cat( " SETTING SEED...\n")
cat( "=================\n")
cat( "The seed used for replicate ", args[2], "is ", seed, "\n\n" )
set.seed( seed )
setwd( wd )

# Set variables for number of characters. Takes argument 3
pvar   <- paste( "p", args[3], sep = "" )
p      <- pvar2 <- as.numeric( args[3] )
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
# This is argument 4
cat( "B. The size of your correlation matrix is ", p, "x", p, "\n\n" )
rho     <- args[4]
rho.wd  <- paste( "rho", args[4], "/", sep = "" )
R       <- matrix(rho, p, p)
diag(R) <- 1

# Set number of specimens to sample for 
# the simulated population
psample = 20

# Get tree 
tree <- read.tree(text="((((A^1:0.1,B^1:0.1):0.2,(F^0.9:0.1,C^1:0.2):0.1):0.5,H^0.3:0.1):0.2,(D^1:0.7,(G^0.7:0.2,E^1:0.5):0.2):0.3);")

# Write each alignment in one file
# We are not running a for loop in the cluster 
# because the array job is like a for loop in which each job
# has an ID. This ID would be the "i" in this for loop,
# so num.i <- args[2]
num.i <- args[2]

cat( "==========================\n")
cat( " GENERATING ALIGNMENTS...\n")
cat( "==========================\n\n")

# Simulate  n = 100 continuous characters with within population variance c = 0.25
# and correlation rho = 0.50.
sim.cont.chars    <- sim.morpho( tree = tree, n = p, c = c, R = R )
  
# Simulate a population sample with psample = 20 specimens, within population 
# variance c = 0.25 and correlation rho = 0.50, and n = 100 characters
cat( "\nSimulating populatin matrix...\nEstimating R.sh...\n" )
simPop            <- sim.pop( psample = psample, n = p, c = c, R = R ) 
  
# Write the alignment that has not been scaled nor transformed, i.e.,
# ignores population noise and character correlation.
write.morpho( M = sim.cont.chars,
              filename = paste( wd, 'alignments/', rho.wd, 'MnR/', pvar, '/rtraitcont', pvar2, '_MnR_',
                                num.i,'.txt', sep='' ),
              names = rownames( sim.cont.chars ) )

# Write alignment after scaling by estimated population noise and matrix transformation
# is carried out with the shrinkage correlation matrix estimated with delta = default mode
write.morpho( M = sim.cont.chars,
              filename = paste( wd, 'alignments/', rho.wd, 'MsR.Rdef/', pvar, '/rtraitcont', pvar2, '_MsR.Rdef_',
                                num.i,'.txt', sep='' ),
              c = simPop$var, R = cbind(simPop$Rsh), method = "chol", names = rownames( sim.cont.chars ) )

# Write alignment after scaling by estimated population noise and matrix transformation
# is carried out with the shrinkage correlation matrix estimated with delta = 0.01 (fixed value)
cat( "\nEstimating R.sh.0.01...\n" )
Rsh.0.01 <- as.matrix( corpcor::cor.shrink( simPop$P, lambda = 0.01 ) )
write.morpho( M = sim.cont.chars,
              filename = paste( wd, 'alignments/', rho.wd, 'MsR.R0.01/', pvar, '/rtraitcont', pvar2, '_MsR.R0.01_',
                                num.i,'.txt', sep='' ),
              c = simPop$var, R = cbind(Rsh.0.01), method = "chol", names = rownames( sim.cont.chars ) )

# Write alignment after scaling by true population noise and matrix transformation is 
# carried out using the true correlation matrix
write.morpho( M = sim.cont.chars,
              filename = paste( wd, 'alignments/', rho.wd, 'MstrueR.Rtrue/', pvar, '/rtraitcont', pvar2, '_MstrueR.Rtrue_',
                                num.i,'.txt', sep='' ),
              c = c, R = R, method = "chol", names = rownames( sim.cont.chars ) )

cat( "All alignments for replicate ", args[2], "have been saved!\n")
cat( "\n*\n*\n*\n*\n\n")