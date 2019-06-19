#----------------------------#
# EFFECT OF POPULATION NOISE #
#----------------------------#

# Before running this code, you should have in your
# working directory the following architecture or something similar.
#
# Note that scaling was performed using the estimate of the population noise 
# calculated on the simulated sampled population so we could simulate 
# what we need to do when real data are used.
# We only scale by the true population noise the data under
# the directory `c_0.25_true` and `c_0.50_true`.
# 
# - your_wd
#     |-alignments                                            # corrected pop.noise?
#              |- c_0.25_Mn       # not scaled            -->     N | c=0 (ignored) 
#              |   |-p100
#              |   |-p1000
#              |   |-p10000
#              |   
#              |- c_0.25          # scaled, using est.c   -->     Y | c=1 (corrected with est.c)
#              |   |-p100
#              |   |-p1000
#              |   |-p10000
#              |   
#              |- c_0.25_true     # scaled, using true.c  -->     Y | c=1 (corrected with true.c)  
#              |   |-p100
#              |   |-p1000
#              |   |-p10000
#              |   
#              |- c_0.50_Mn       # not scaled            -->     N | c=0 (ignored) 
#              |   |-p100
#              |   |-p1000
#              |   |-p10000
#              |   
#              |- c_0.50          # scaled, using est.c   -->     Y | c=1 (corrected with est.c)
#              |   |-p100
#              |   |-p1000
#              |   |-p10000
#              |   
#              |- c_0.50_true     # scaled, using true c  -->     Y | c=1 (corrected with true.c)    
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
#   $ RScript 02_Simulate_popnoise.R <path_to_wd> 
#
# E.g.: 
#   $ RScript 02_Simulate_popnoise.R /home/morpho/p100/ 

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

# Set number of specimens to sample
# from the population
psample <- 20 

# Set population noise
c.low  <- 0.25  # Small population noise
c.high <- 0.50 # Large population noise

# Define tree 
tree <- ape::read.tree(text="((((A^1:0.1,B^1:0.1):0.2,(F^0.9:0.1,C^1:0.2):0.1):0.5,H^0.3:0.1):0.2,(D^1:0.7,(G^0.7:0.2,E^1:0.5):0.2):0.3);")

#---------------------------#
# -- p = 100 characters --  #
#---------------------------#

# Write each alignment in one file - p100, c= 0.25
for ( i in seq( 1:reps ) ){

  # Simulate  n = 100 continuous characters with within population variance c = 0.25
  sim.cont.chars    <- mcmc3r::sim.morpho( tree = tree, n = p1, c = c.low )

  # Simulate a population sample with psample = 20 specimens, within population 
  # variance c = 0.25 and n = 100 characters
  simPop            <- mcmc3r::sim.pop( psample = psample, n = p1, c = c.low )

  # Write alignment that has not been scaled nor corrected for popnoise,
  # i.e., it ignores population noise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.25_Mn/p', p1, '/rtraitcont', p1, '_Mn_', i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )

  # Write alignment that has been scaled and corrected for estimated popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.25/p', p1, '/rtraitcont', p1, '_Ms_', i,'.txt', sep='' ),
                        c = simPop$var, names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for true popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.25_true/p', p1, '/rtraitcont', p1, '_Mstrue_', i,'.txt', sep='' ),
                        c = c.low, names = rownames( sim.cont.chars ) )

}

# Write each alignment in one file - p100, c= 0.50
for ( i in seq( 1:reps ) ){
  
  # Simulate  n = 100 continuous characters with within population variance c = 0.50
  sim.cont.chars    <- mcmc3r::sim.morpho( tree = tree, n = p1, c = c.high )
  
  # Simulate a population sample with psample = 20 specimens, within population 
  # variance c = 0.50 and n = 100 characters
  simPop            <- mcmc3r::sim.pop( psample = psample, n = p1, c = c.high )
  
  # Write alignment that has not been scaled nor corrected for popnoise
  # i.e., it ignores population noise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.50_Mn/p', p1, '/rtraitcont', p1, '_Mn_', i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for estimated popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.50/p', p1, '/rtraitcont', p1, '_Ms_', i,'.txt', sep='' ),
                        c = simPop$var, names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for true popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.50_true/p', p1, '/rtraitcont', p1, '_Mstrue_', i,'.txt', sep='' ),
                        c = c.high, names = rownames( sim.cont.chars ) )
  
}

#----------------------------#
# -- p = 1,000 characters -- #
#----------------------------#

# Write each alignment in one file - p1000, c= 0.25
for ( i in seq( 1:reps ) ){
  
  # Simulate  n = 1000 continuous characters with within population variance c = 0.25
  sim.cont.chars    <- mcmc3r::sim.morpho( tree = tree, n = p2, c = c.low )
  
  # Simulate a population sample with psample = 20 specimens, within population 
  # variance c = 0.25 and n = 1000 characters
  simPop            <- mcmc3r::sim.pop( psample = psample, n = p2, c = c.low )
  
  # Write alignment that has not been scaled nor corrected for popnoise
  # i.e., it ignores population noise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.25_Mn/p', p2, '/rtraitcont', p2, '_Mn_', i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for estimated popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.25/p', p2, '/rtraitcont', p2, '_Ms_', i,'.txt', sep='' ),
                        c = simPop$var, names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for true popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.25_true/p', p2, '/rtraitcont', p2, '_Mstrue_', i,'.txt', sep='' ),
                        c = c.low, names = rownames( sim.cont.chars ) )
  
}

# Write each alignment in one file - p1000, c= 0.50
for ( i in seq( 1:reps ) ){
  
  # Simulate  n = 1000 continuous characters with within population variance c = 0.50
  sim.cont.chars    <- mcmc3r::sim.morpho( tree = tree, n = p2, c = c.high )
  
  # Simulate a population sample with psample = 20 specimens, within population 
  # variance c = 0.50 and n = 1000 characters
  simPop            <- mcmc3r::sim.pop( psample = psample, n = p2, c = c.high )
  
  # Write alignment that has not been scaled nor corrected for popnoise
  # i.e., it ignores population noise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.50_Mn/p', p2, '/rtraitcont', p2, '_Mn_', i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for estimated popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.50/p', p2, '/rtraitcont', p2, '_Ms_', i,'.txt', sep='' ),
                        c = simPop$var, names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for true popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.50_true/p', p2, '/rtraitcont', p2, '_Mstrue_', i,'.txt', sep='' ),
                        c = c.high, names = rownames( sim.cont.chars ) )
  
}



#-----------------------------#
# -- p = 10,000 characters -- #
#-----------------------------#

# Write each alignment in one file - p10000, c= 0.25
for ( i in seq( 1:reps ) ){
  
  # Simulate  n = 10000 continuous characters with within population variance c = 0.25
  sim.cont.chars    <- mcmc3r::sim.morpho( tree = tree, n = p3, c = c.low )
  
  # Simulate a population sample with psample = 20 specimens, within population 
  # variance c = 0.25 and n = 10000 characters
  simPop            <- mcmc3r::sim.pop( psample = psample, n = p3, c = c.low )
  
  # Write alignment that has not been scaled nor corrected for popnoise
  # i.e., it ignores population noise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.25_Mn/p', p3, '/rtraitcont', p3, '_Mn_', i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for estimated popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.25/p', p3, '/rtraitcont', p3, '_Ms_', i,'.txt', sep='' ),
                        c = simPop$var, names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for true popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.25_true/p', p3, '/rtraitcont', p3, '_Mstrue_', i,'.txt', sep='' ),
                        c = c.low, names = rownames( sim.cont.chars ) )
  
}

# Write each alignment in one file - p10000, c= 0.50
for ( i in seq( 1:reps ) ){
  
  # Simulate  n = 10000 continuous characters with within population variance c = 0.50
  sim.cont.chars    <- mcmc3r::sim.morpho( tree = tree, n = p3, c = c.high )
  
  # Simulate a population sample with psample = 20 specimens, within population 
  # variance c = 0.50 and n = 10000 characters
  simPop            <- mcmc3r::sim.pop( psample = psample, n = p3, c = c.high )
  
  # Write alignment that has not been scaled nor corrected for popnoise
  # i.e., it ignores population noise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.50_Mn/p', p3, '/rtraitcont', p3, '_Mn_', i,'.txt', sep='' ),
                        names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for est.popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.50/p', p3, '/rtraitcont', p3, '_Ms_', i,'.txt', sep='' ),
                        c = simPop$var, names = rownames( sim.cont.chars ) )
  
  # Write alignment that has been scaled and corrected for true popnoise
  mcmc3r::write.morpho( M = sim.cont.chars,
                        filename = paste( wd, 'c_0.50_true/p', p3, '/rtraitcont', p3, '_Mstrue_', i,'.txt', sep='' ),
                        c = c.high, names = rownames( sim.cont.chars ) )
  
}
