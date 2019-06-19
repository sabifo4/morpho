#----------------------#
# EFFECT OF FOSSIL AGE #
#----------------------#

# Before running this code, you should have in your
# working directory the following architecture or something similar.
#
# - your_wd
#     |- spH.1_0.7      
#     |        |-p100
#     |        |-p1000
#     |        |-p10000
#     |   
#     |- spH.2_0.5      
#     |        |-p100
#     |        |-p1000
#     |        |-p10000
#     |   
#     |- spH.3_0.3      
#     |        |-p100
#     |        |-p1000
#     |        |-p10000
#     |   
#     |- spH.4_0.1      
#              |-p100
#              |-p1000
#              |-p10000
#
# If you want to change your file architecture, change
# the path to the output file in the "mmcmc3::write.morpho"
# function in the "filename" option. This is used to save
# the alignments following the architecture
# described above. 
#
# This script should be run as it follows: 
#   $ RScript 01_simulate_fossil_age.R <path_to_wd>
#
# E.g.: RScript 01_simulate_fossil_age.R /home/morpho/fossil_age

   
# Clean environment
rm( list = ls( ) )

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

# Set variables for number of characters
p1 <- 100   # 100 characters
p2 <- 1000  # 1,000 characters
p3 <- 10000 # 10,000 character

# Set number of replicates 
reps <- 1000

# Set number of species 
sp <- 8

# Define trees 
tree.spH1 <- ape::read.tree(text="((((A^1:0.1,B^1:0.1):0.2,(F^0.9:0.1,C^1:0.2):0.1):0.5,H^0.3:0.1):0.2,(D^1:0.7,(G^0.7:0.2,E^1:0.5):0.2):0.3);")
tree.spH2 <- ape::read.tree(text="((((A^1:0.1,B^1:0.1):0.2,(F^0.9:0.1,C^1:0.2):0.1):0.5,H^0.5:0.3):0.2,(D^1:0.7,(G^0.7:0.2,E^1:0.5):0.2):0.3);")
tree.spH3 <- ape::read.tree(text="((((A^1:0.1,B^1:0.1):0.2,(F^0.9:0.1,C^1:0.2):0.1):0.5,H^0.7:0.5):0.2,(D^1:0.7,(G^0.7:0.2,E^1:0.5):0.2):0.3);")
tree.spH4 <- ape::read.tree(text="((((A^1:0.1,B^1:0.1):0.2,(F^0.9:0.1,C^1:0.2):0.1):0.5,H^0.9:0.7):0.2,(D^1:0.7,(G^0.7:0.2,E^1:0.5):0.2):0.3);")

#--------------#
# -- TREE 1 -- #
#--------------#

# Write each alignment in one file - p100
cat( "==========================================================\n")
cat( " TREE 1: GENERATING ALIGNMENTS WITH p = 100 characters...\n")
cat( "==========================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH1, n = p1 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.1_0.7/p100/rtraitcont100_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

# Write each alignment in one file - p1000
cat( "=============================================================\n")
cat( " TREE 1: GENERATING ALIGNMENTS WITH p = 1,000 characters...\n")
cat( "=============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH1, n = p2 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.1_0.7/p1000/rtraitcont1000_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

# Write each alignment in one file - p10000
cat( "==============================================================\n")
cat( " TREE 1: GENERATING ALIGNMENTS WITH p = 10,000 characters...\n")
cat( "==============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH1, n = p3 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.1_0.7/p10000/rtraitcont10000_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

#--------------#
# -- TREE 2 -- #
#--------------#

# Write each alignment in one file - p100
cat( "=============================================================\n")
cat( " TREE 2: GENERATING ALIGNMENTS WITH p = 100 characters...\n")
cat( "=============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH2, n = p1 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.2_0.5/p100/rtraitcont100_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

# Write each alignment in one file - p1000
cat( "=============================================================\n")
cat( " TREE 2: GENERATING ALIGNMENTS WITH p = 1,000 characters...\n")
cat( "=============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH2, n = p2 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.2_0.5/p1000/rtraitcont1000_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

# Write each alignment in one file - p10000
cat( "==============================================================\n")
cat( " TREE 2: GENERATING ALIGNMENTS WITH p = 10,000 characters...\n")
cat( "==============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH2, n = p3 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.2_0.5/p10000/rtraitcont10000_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

#--------------#
# -- TREE 3 -- #
#--------------#

# Write each alignment in one file - p100
cat( "=============================================================\n")
cat( " TREE 3: GENERATING ALIGNMENTS WITH p = 100 characters...\n")
cat( "=============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH3, n = p1 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.3_0.3/p100/rtraitcont100_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

# Write each alignment in one file - p1000
cat( "=============================================================\n")
cat( " TREE 3: GENERATING ALIGNMENTS WITH p = 1,000 characters...\n")
cat( "=============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH3, n = p2 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.3_0.3/p1000/rtraitcont1000_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

# Write each alignment in one file - p10000
cat( "=============================================================\n")
cat( " TREE 3: GENERATING ALIGNMENTS WITH p = 10,000 characters...\n")
cat( "=============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH3, n = p3 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.3_0.3/p10000/rtraitcont10000_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

#--------------#
# -- TREE 4 -- #
#--------------#

# Write each alignment in one file - p100
cat( "=============================================================\n")
cat( " TREE 4: GENERATING ALIGNMENTS WITH p = 100 characters...\n")
cat( "=============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH4, n = p1 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.4_0.1/p100/rtraitcont100_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
  
}

# Write each alignment in one file - p1000
cat( "=============================================================\n")
cat( " TREE 4: GENERATING ALIGNMENTS WITH p = 1,000 characters...\n")
cat( "=============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH4, n = p2 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.4_0.1/p1000/rtraitcont1000_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}

# Write each alignment in one file - p10000
cat( "=============================================================\n")
cat( " TREE 4: GENERATING ALIGNMENTS WITH p = 10,000 characters...\n")
cat( "=============================================================\n\n")
for ( i in seq( 1:reps ) ){
  sim.cont.chars <- mcmc3r::sim.morpho( tree = tree.spH4, n = p3 )
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( 'spH.4_0.1/p10000/rtraitcont10000_',
                                          i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
}
