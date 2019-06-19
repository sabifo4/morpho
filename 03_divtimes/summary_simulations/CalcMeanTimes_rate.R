## ------------------------------------------------------------ ##
## ----    AVERAGE OVER POSTERIOR MEAN RATE ESTIMATES    ------ ##
## ------------------------------------------------------------ ##
## This script is based on Kostas Angelis' code.                ##
## It was first used to summarise the results obtained in       ##
## Angelis et al., 2017, https://doi.org/10.1093/sysbio/syx061  ##
## ------------------------------------------------------------ ##
## Estimate averages of posterior mean rate estimates and CIs   ##
## over all replicates for each node separately, for each       ##
## combination of rate prior, calibration and rate-drift model. ##
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
source( "/SumResults_rate.R" )

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
rates <- c( "mu" )
estimates <- c( "mean", "L", "U" )
rep.name  <- paste( rep( "Rep"), 1:reps, sep = "" )


## ----------------------------------- ##
## ----    GET MEAN POST. RATE    ---- ##
## ----------------------------------- ##

# Create 3D arrays for the different simulations

# Replicates ( p = 100, 1000, 10000 characters)
rt.p100_ar    <- array( unlist( rt.p100 ), dim = c( 1, 3, reps ),
                              dimnames = list( rates, estimates, rep.name ) )
rt.p1000_ar   <- array( unlist( rt.p1000 ), dim = c( 1, 3, reps ),
                              dimnames = list( rates, estimates, rep.name ) )
rt.p10000_ar  <- array( unlist( rt.p10000 ), dim = c( 1, 3, reps ),
                               dimnames = list( rates, estimates, rep.name ) )

#  Estimate mean across all 1000 replicates

rt.p100_mean    <- apply( rt.p100_ar, 1:2, mean )
rt.p1000_mean   <- apply( rt.p1000_ar, 1:2, mean )
rt.p10000_mean  <- apply( rt.p10000_ar, 1:2, mean )

# Calculate quantiles 

rt.p100_quantiles <- apply( rt.p100_ar, 1:2, function( x ) quantile( x, c( .025,.975 ) ) )[,,1]

rt.p1000_quantiles <- apply( rt.p1000_ar, 1:2, function( x ) quantile( x, c( .025,.975 ) ) )[,,1]

rt.p10000_quantiles <- apply( rt.p10000_ar, 1:2, function( x ) quantile( x, c( .025,.975 ) ) )[,,1]


# Format table

summ <- summ2 <- matrix( 0, nrow = 3, ncol = 1 )

# Use mean for L and U

summ[1,1] <- paste( "\"", round( rt.p100_mean[1], 3 ), " (", 
                     round( rt.p100_mean[2], 3 ), ",",  
                     round( rt.p100_mean[3], 3 ), ")", "\"", sep="" )
summ[2,1] <- paste( "\"", round( rt.p1000_mean[1], 3 ), " (", 
                     round( rt.p1000_mean[2], 3 ), ",",  
                     round( rt.p1000_mean[3], 3 ), ")", "\"", sep="" )
summ[3,1] <- paste( "\"", round( rt.p10000_mean[1], 3 ), " (", 
                      round( rt.p10000_mean[2], 3 ), ",",  
                      round( rt.p10000_mean[3], 3 ), ")", "\"", sep="" )

# Use quantiles for L and U

summ2[1,1] <- paste( "\"", round( rt.p100_mean[1], 3 ), " (", 
                    round( rt.p100_quantiles[1], 3 ), ",",  
                    round( rt.p100_quantiles[2], 3 ), ")", "\"", sep="" )
summ2[2,1] <- paste( "\"", round( rt.p1000_mean[1], 3 ), " (", 
                    round( rt.p1000_quantiles[1], 3 ), ",",  
                    round( rt.p1000_quantiles[2], 3 ), ")", "\"", sep="" )
summ2[3,1] <- paste( "\"", round( rt.p10000_mean[1], 3 ), " (", 
                     round( rt.p10000_quantiles[1], 3 ), ",",  
                     round( rt.p10000_quantiles[2], 3 ), ")", "\"", sep="" )

rownames( summ ) <- rownames( summ2 ) <- c( "r.100", "r.1000", "r.10000" )


## -------------------------------- ##
## ----  SUMMARISE  RESULTS    ---- ##
## -------------------------------- ##

write.csv( summ,
           paste( wd, "mean_post_times_rates.csv", sep = '' ), quote = F )
write.csv( summ2,
           paste( wd, "mean_post_times_rates_quantiles.csv", sep = '' ), quote= F )

## ----------------------- ##
## *********************** ##
