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
##   the true age ti                                ##
## ------------------------------------------------ ##

## -------------------------------------- ##
## ----    INITIALIZE ENVIRONMENT   ----- ##
## -------------------------------------- ##

# Clean environment
rm( list = ls() )

# Read estimates
# NOTE: Change the path according to the directory where you have the analysis running!
cat("Load estimates...\n")
source( "SumResults.R" )

# Read functions file
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
# True ages (in 100My)
#     Order:  t9, t10, t11, t12, t13, t14, t15
trueAges <- c( 1, 0.8, 0.3, 0.1, 0.2, 0.7, 0.5 )
names( trueAges ) <- paste( rep( "t_n", sp-1 ), (sp+1):(2*sp-1), sep = "" )
# Node names
nodeNames <- paste( rep( "n", sp-1 ), (sp+1):(2*sp-1), sep = "" )


## --------------------------------------------------------- ##
## ----    FIND WHETHER A CI CONTAINS THE TRUE VALUE    ---- ##
## --------------------------------------------------------- ##

cat("Start analysis!\n")

# Create matrices to store the coverage
cover_rt.p100      <- cover_rt.p1000      <- cover_rt.p10000 <-
matrix( numeric( reps * (sp-1) ), reps, sp-1 )

# Look if the true age is in the 95%CI interval estimated by MCMCTree
# with the simulated data by using the function xInRange
# This function is inside the "functions.R" script
# used in Konstas et. al, 2017, https://doi.org/10.1093/sysbio/syx061
for( j in 1 : reps ){

  cover_rt.p100[ j , ]        <- xInRange( trueAges, rt.p100[[ j ]][ , 2:3 ] )
  cover_rt.p1000[ j , ]       <- xInRange( trueAges, rt.p1000[[ j ]][ , 2:3 ] )
  cover_rt.p10000[ j , ]      <- xInRange( trueAges, rt.p10000[[ j ]][ , 2:3 ] )

}

# Create matrices to later fill in with probabilities
# of node i containing the true age ti in each simulation
prob_rt.p100_node <- prob_rt.p1000_node  <- prob_rt.p10000_node <-
matrix( numeric( 1 * (sp-1) ), 1, sp-1 )

# Use nodeNames for the colnames in probability matrices
colnames( prob_rt.p100_node ) <- colnames( prob_rt.p1000_node ) <- colnames( prob_rt.p10000_node ) <-
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
summ <- matrix( 0, nrow = 28, ncol = 1 )
j=0
for (i in 1:7){
   
  summ[1+j,1] <- paste( prob_rt.p100_node[i,1], sep="" )
  summ[2+j,1] <- paste( prob_rt.p1000_node[i,1], sep="" )
  summ[3+j,1] <- paste( prob_rt.p10000_node[i,1], sep="" )
  summ[4+j,1] <- ""
  
  # Keep adding 4 more
  j <- j + 4
  
}

rownames( summ ) <- c( rep( "t9", 3 ), "", rep( "t10", 3 ), "", rep( "t11", 3 ), "", 
                       rep( "t12", 3 ), "", rep( "t13", 3 ), "", rep( "t14", 3 ), "", 
                       rep( "t15", 3 ), "")


## ------------------------- ##
## ----    RESULTS    ------ ##
## ------------------------- ##

prob_node             <- cbind( prob_rt.p100_node, prob_rt.p1000_node, prob_rt.p10000_node )
colnames( prob_node ) <- c( "rt100", "rt1000", "rt10000" )

write.csv( summ, paste( wd, "coverage_node.csv", sep = '' ), quote = F )

## -------------------------- ##
## ************************** ##
