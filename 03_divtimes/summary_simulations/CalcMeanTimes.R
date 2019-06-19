## ------------------------------------------------------------ ##
## ----    AVERAGE OVER POSTERIOR MEAN TIME ESTIMATES    ------ ##
## ------------------------------------------------------------ ##
## This script is based on Kostas Angelis' code.                ##
## It was first used to summarise the results obtained in       ##
## Angelis et al., 2017, https://doi.org/10.1093/sysbio/syx061  ##
## ------------------------------------------------------------ ##
## Estimate averages of posterior time means and CIs over all   ##
## replicates for each node separately, for each combination of ##
## rate prior, calibration and rate-drift model.                ##
## ------------------------------------------------------------ ##


## -------------------------------------- ##
## ----    INITIALIZE ENVIRONMENT   ----- ##
## -------------------------------------- ##

# Clean environment
rm( list = ls() )

# Load libraries
library( ape )

# Read estimates
# NOTE: Change the path according to the directory where you have the analysis running!
cat("Load estimates...\n")
source( "SumResults.R" )

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
# Names for the rows and cols of the 3D arrays
nodes <- c( "t9", "t10", "t11", "t12", "t13", "t14", "t15" )
estimates <- c( "mean", "L", "U" )
rep.name  <- paste( rep( "Rep"), 1:reps, sep = "" )


## ----------------------------------- ##
## ----    GET MEAN POST.TIMES    ---- ##
## ----------------------------------- ##

# Create 3D arrays for the different simulations

# Replicates ( p = 100, 1000, 10000 characters)
rt.p100_ar    <- array( unlist( rt.p100 ), dim = c( (sp-1), 3, reps ),
                        dimnames = list( nodes, estimates, rep.name ) )
rt.p1000_ar   <- array( unlist( rt.p1000 ), dim = c( (sp-1), 3, reps ),
                        dimnames = list( nodes, estimates, rep.name ) )
rt.p10000_ar  <- array( unlist( rt.p10000 ), dim = c( (sp-1), 3, reps ),
                         dimnames = list( nodes, estimates, rep.name ) )


# Estimate mean across all 1000 replicates

rt.p100_mean    <- apply( rt.p100_ar, 1:2, mean )
rt.p1000_mean   <- apply( rt.p1000_ar, 1:2, mean )
rt.p10000_mean  <- apply( rt.p10000_ar, 1:2, mean )


# Calculate quantiles 
rt.p100_quantiles <- apply( rt.p100_ar, 1:2, function( x ) quantile( x,c( .025,.975 ) ) )[,,1]

rt.p1000_quantiles <- apply( rt.p1000_ar, 1:2, function( x ) quantile( x, c( .025,.975 ) ) )[,,1]

rt.p10000_quantiles <- apply( rt.p10000_ar, 1:2, function( x ) quantile( x, c( .025,.975 ) ) )[,,1]


# Format table to copy paste in the formt of the SS material

summ <- summ2 <- matrix( 0, nrow = 28, ncol = 1 )

# Use mean for L and U

j = 0
for ( i in 1:7 ){
  
  summ[1+j,1] <- paste( "\"", round( rt.p100_mean[i,1], 3 ), " (", 
                        round( rt.p100_mean[i,2], 3 ), ",",  
                        round( rt.p100_mean[i,3], 3 ), ")", "\"", sep="" )
  summ[2+j,1] <- paste( "\"", round( rt.p1000_mean[i,1], 3 ), " (", 
                        round( rt.p1000_mean[i,2], 3 ), ",",  
                        round( rt.p1000_mean[i,3], 3 ), ")", "\"", sep="" )
  summ[3+j,1] <- paste( "\"", round( rt.p10000_mean[i,1], 3 ), " (", 
                         round( rt.p10000_mean[i,2], 3 ), ",",  
                         round( rt.p10000_mean[i,3], 3 ), ")", "\"", sep="" )
  summ[4+j,1] <- ""
  
  # Keep adding 4 more
  j <- j + 4
  
}

# Now use quantiles
j = 0
for ( i in 1:7 ){
  
  summ2[1+j,1] <- paste( "\"", round( rt.p100_mean[i,1], 3 ), " (", 
                        round( rt.p100_quantiles[1,i], 3 ), ",",  
                        round( rt.p100_quantiles[2,i], 3 ), ")", "\"", sep="" )
  summ2[2+j,1] <- paste( "\"", round( rt.p1000_mean[i,1], 3 ), " (", 
                        round( rt.p1000_quantiles[1,i], 3 ), ",",  
                        round( rt.p1000_quantiles[2,i], 3 ), ")", "\"", sep="" )
  summ2[3+j,1] <- paste( "\"", round( rt.p10000_mean[i,1], 3 ), " (", 
                         round( rt.p10000_quantiles[1,i], 3 ), ",",  
                         round( rt.p10000_quantiles[2,i], 3 ), ")", "\"", sep="" )
  summ2[4+j,1] <- ""
  
  # Keep adding 4 more
  j <- j + 4
  
}

rownames( summ ) <- rownames( summ ) <- 
  c( rep( "t9", 3 ), "", rep( "t10", 3 ), "", rep( "t11", 3 ), "", 
     rep( "t12", 3 ), "", rep( "t13", 3 ), "", rep( "t14", 3 ), "", 
     rep( "t15", 3 ), "" )


## ------------------------------- ##
## ----  SUMMARISE RESULTS    ---- ##
## ------------------------------- ##

write.csv( summ,
           paste( wd, "mean_post_times.csv", sep = '' ), quote = F )
write.csv( summ2,
           paste( wd, "mean_post_times_quantiles.csv", sep = '' ), quote = F )

## ----------------------- ##
## *********************** ##
