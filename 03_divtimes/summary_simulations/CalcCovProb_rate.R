## ------------------------------------------------ ##
## ----    COMPARE COVERAGE PROBABILITIES    ------ ##
## ------------------------------------------------ ##
## This script is based on Kostas Angelis' code.    ##
## It was first used to summarise the results       ##
## obtained in Angelis et al., 2017,                ## 
## https://doi.org/10.1093/sysbio/syx061            ##
## ------------------------------------------------ ##
##   Calculate the percentage of replicates which   ##
##   credibility interval (CI) of node i contains   ##
##   the true rate                                  ##
## ------------------------------------------------ ##

## -------------------------------------- ##
## ----    INITIALIZE ENVIRONMENT   ----- ##
## -------------------------------------- ##

# Clean environment
rm( list = ls() )

# Read estimates
# NOTE: Change the path according to the directory where you have the analysis running!
cat("Load estimates...\n")
source( "SumResults_rate.R" )

# Read functions Konstas used
# NOTE: You can change the location to this "functions.R" script, but then remember to 
# change this path
source( "functions.R" )

# Get wd
# NOTE: Remember to change this to the path where you want to save the summary results!
wd <- "my_wd"


## -------------------------------- ##
## ----    CREATE VARIABLES    ---- ##
## -------------------------------- ##

# Num. replicates
reps <- 1000
# Num. species
sp <- 8
# True rate
trueRate <- c( 1 )
names( trueRate ) <- paste( "mu" )
# Node names
nodeNames <- paste( "mu" )


## --------------------------------------------------------- ##
## ----    FIND WHETHER A CI CONTAINS THE TRUE VALUE    ---- ##
## --------------------------------------------------------- ##

cat("Start analysis!\n")

# Create matrices to store the coverage
cover_rt.p100      <- cover_rt.p1000      <- cover_rt.p10000 <-
matrix( numeric( reps * (1) ), reps, 1 )

# Look if the true age is in the 95%CI interval estimated by MCMCTree
# with the simulated data by using the function xInRange
# This function is inside the "functions.R" script
# used in Konstas et. al, 2017, https://doi.org/10.1093/sysbio/syx061
for( j in 1 : reps ){

   cover_rt.p100[ j , ]        <- xInRange( trueRate, rt.p100[[ j ]][ , 2:3 ] )
   cover_rt.p1000[ j , ]       <- xInRange( trueRate, rt.p1000[[ j ]][ , 2:3 ] )
   cover_rt.p10000[ j , ]      <- xInRange( trueRate, rt.p10000[[ j ]][ , 2:3 ] )

}

# Create matrices to later fill in with probabilities
# of node i containing the true age ti in each simulation
prob_rt.p100_node      <- prob_rt.p1000_node      <- prob_rt.p10000_node <-
matrix( numeric( 1 * (1) ), 1, 1 )

# Use nodeNames for the colnames in probability matrices
colnames( prob_rt.p100_node )      <- colnames( prob_rt.p1000_node )      <- colnames( prob_rt.p10000_node ) <-
nodeNames

# Calculate the coverage (in %)
# Mean of '1' (true, CI contains true value) and '0' values
# (false, CI does not contain true value) in cover_* tables
# which were calculated by xInRange

# Replicates (100, 1000, 10000)
prob_rt.p100_node[ 1, ]   <- round( colMeans( cover_rt.p100 ), 3 )
prob_rt.p100_node         <- t( prob_rt.p100_node )
prob_rt.p100_node         <- as.data.frame( 100 * prob_rt.p100_node )

prob_rt.p1000_node[ 1, ]  <- round( colMeans( cover_rt.p1000 ), 3 )
prob_rt.p1000_node        <- t( prob_rt.p1000_node )
prob_rt.p1000_node        <- as.data.frame( 100 * prob_rt.p1000_node )

prob_rt.p10000_node[ 1, ] <- round( colMeans( cover_rt.p10000 ), 3 )
prob_rt.p10000_node       <- t( prob_rt.p10000_node )
prob_rt.p10000_node       <- as.data.frame( 100 * prob_rt.p10000_node )


# Convert prob. tables in nice format
summ <- matrix( 0, nrow = 3, ncol = 1 )
summ[1,1] <- paste( prob_rt.p100_node[1], sep = "" )
summ[2,1] <- paste( prob_rt.p1000_node[1], sep = "" )
summ[3,1] <- paste( prob_rt.p10000_node[1], sep = "" )

rownames( summ ) <- c( "r.100", "r.1000", "r.10000" )


## ------------------------- ##
## ----    RESULTS    ------ ##
## ------------------------- ##

prob_node             <- cbind( prob_rt.p100_node, prob_rt.p1000_node, prob_rt.p10000_node )
colnames( prob_node ) <- c( "rt100", "rt1000", "rt10000" )

write.csv( summ, paste( wd, "coverage_node_rate.csv", sep = '' ), quote = F )

## -------------------------- ##
## ************************** ##

