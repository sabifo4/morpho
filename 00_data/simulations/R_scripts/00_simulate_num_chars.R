#--------------------------------#
# EFFECT OF NUMBER OF CHARACTERS #
#--------------------------------#

# Before running this code, you should have in your
# working directory the following architecture or something similar.
#
# - your_wd
#     |-p100
#     |-p1000
#     |-p10000
#
# This script should be run as it follows: 
#   $ RScript 00_simulate_num_chars.R <path_to_wd>
#
# E.g.: RScript 00_simulate_num_chars.R /home/morpho/num_chars

# Collect args 
args  <- commandArgs(TRUE)

# Set working directory
wd <- args[1]
setwd( wd )

# Load libraries 
library( ape )
library( mcmc3r )

# Set seed 
set.seed( 12345 )

# Set variables for replicates
p1 <- 100   # 100 replicates
p2 <- 1000  # 1,000 replicates
p3 <- 10000 # 10,000 replicates

# Set number of replicates 
reps <- 1000

# Set number of species 
sp <- 8

# Define tree 
tree <- ape::read.tree(text="((((A^1:0.1,B^1:0.1):0.2,(F^0.9:0.1,C^1:0.2):0.1):0.5,H^0.3:0.1):0.2,(D^1:0.7,(G^0.7:0.2,E^1:0.5):0.2):0.3);")

# Write each alignment in one file - p100
cat( "====================================================\n")
cat( " GENERATING ALIGNMENTS WITH p = 100 characters...\n")
cat( "====================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree, n = p1 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'p100/rtraitcont100_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

# Write each alignment in one file - p1000
cat( "====================================================\n")
cat( " GENERATING ALIGNMENTS WITH p = 1,000 characters...\n")
cat( "====================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree, n = p2 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'p1000/rtraitcont1000_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

# Write each alignment in one file - p10000
cat( "====================================================\n")
cat( " GENERATING ALIGNMENTS WITH p = 10,000 characters...\n")
cat( "====================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree, n = p3 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'p10000/rtraitcont10000_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}
