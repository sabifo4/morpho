## ------------------------------------------------ ##
## ----    MSE, BIAS, and WIDTH CALCULATION    ---- ##
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

# Read estimates
# NOTE: Change the path according to the directory where you have the analysis running!
cat("Load estimates...\n")
source( "SumResults_rate.R" )
source( "CalcMeanTimes_rate.R" )

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
trueRate <- c( 1 )
names( trueRate ) <- paste( "mu" )
# Node names
nodeNames <- paste( "mu" )


## ------------------------------------------ ##
## ----    CALCULATE MSE, BIAS, WIDTH    ---- ##
## ------------------------------------------ ##

cat("Start analysis!\n")

# Create matrices to store the MSE-lyx (epsilon)
mse_rt.p100  <- mse_rt.p1000  <- mse_rt.p10000  <-
matrix( numeric( reps * (1) ), reps, 1 )

# Create matrices to store the Bias (b)
bias_rt.p100  <- bias_rt.p1000  <- bias_rt.p10000  <-
matrix( numeric( reps * (1) ), reps, 1 )

# Create matrices to store the Width (w)
width_rt.p100 <- width_rt.p1000 <- width_rt.p10000 <-
matrix( numeric( reps * (1) ), reps, 1 )

# Calculate the width, the bias, and the MSE per replicate
# and store the corresponding values in the matrices
# created before
for( j in 1 : reps ){

  # Go through replicates ( p = 100, 1000, 10000 chars)
  
  # p = 100 characters
  w    <- rt.p100[[ j ]]$U - rt.p100[[ j ]]$L
  var  <- ( w / ( 2 * qnorm( 0.975 ) ) ) ^ 2
  bias <- rt.p100[[ j ]]$mean - trueRate
  mse_rt.p100[ j , ]   <- bias ^ 2
  bias_rt.p100[ j , ]  <- bias
  width_rt.p100[ j , ] <- w
  
  # p = 1,000 characters
  w    <- rt.p1000[[ j ]]$U - rt.p1000[[ j ]]$L
  var  <- ( w / ( 2 * qnorm( 0.975 ) ) ) ^ 2
  bias <- rt.p1000[[ j ]]$mean - trueRate
  mse_rt.p1000[ j , ]   <- bias ^ 2
  bias_rt.p1000[ j , ]  <- bias
  width_rt.p1000[ j , ] <- w
  
  # p = 10,000 characters
  w    <- rt.p10000[[ j ]]$U - rt.p10000[[ j ]]$L
  var  <- ( w / ( 2 * qnorm( 0.975 ) ) ) ^ 2
  bias <- rt.p10000[[ j ]]$mean - trueRate
  mse_rt.p10000[ j , ]   <- bias ^ 2
  bias_rt.p10000[ j , ]  <- bias
  width_rt.p10000[ j , ] <- w

}


# Create matrices to store MSE, bias, and width later

# MSE matrices
mse_rt.p100_node    <- mse_rt.p1000_node    <- mse_rt.p10000_node    <-
matrix( numeric( 1 * (1) ), 1, 1 )
# Relative MSE matrices
relmse_rt.p100_node <- relmse_rt.p1000_node <- relmse_rt.p10000_node <-
matrix( numeric( 1 * (1) ), 1, 1 )

# Bias matrices
bias_rt.p100_node    <- bias_rt.p1000_node    <- bias_rt.p10000_node    <-
matrix( numeric( 1 * (1) ), 1, 1 )
# Relative Bias matrices
relbias_rt.p100_node <- relbias_rt.p1000_node <- relbias_rt.p10000_node <-
matrix( numeric( 1 * (1) ), 1, 1 )

# Width matrices
width_rt.p100_node    <- width_rt.p1000_node    <- width_rt.p10000_node    <-
matrix( numeric( 1 * (1) ), 1, 1 )
# Relative Width matrices
relwidth_rt.p100_node <- relwidth_rt.p1000_node <- relwidth_rt.p10000_node <-
matrix( numeric( 1 * (1) ), 1, 1 )


# Add colnames to the matrices created before

# MSE
colnames( mse_rt.p100_node )    <- colnames( mse_rt.p1000_node )    <- colnames( mse_rt.p10000_node )    <-
nodeNames
# Relative MSE
colnames( relmse_rt.p100_node ) <- colnames( relmse_rt.p1000_node ) <- colnames( relmse_rt.p10000_node ) <-
nodeNames

# Bias
colnames( bias_rt.p100_node )    <- colnames( bias_rt.p1000_node )    <- colnames( bias_rt.p10000_node )    <-
nodeNames
# Relative Bias
colnames( relbias_rt.p100_node ) <- colnames( relbias_rt.p1000_node ) <- colnames( relbias_rt.p10000_node ) <-
nodeNames

# Width
colnames( width_rt.p100_node )    <- colnames( width_rt.p1000_node )    <- colnames( width_rt.p10000_node )    <-
nodeNames
# Relative Width
colnames( relwidth_rt.p100_node ) <- colnames( relwidth_rt.p1000_node ) <- colnames( relwidth_rt.p10000_node ) <-
nodeNames


# Average MSE, bias, and width for each node across replicates

# -- MSE
mse_rt.p100_node[ 1, ]   <- round( colMeans( mse_rt.p100 ), 6 )
mse_rt.p100_node         <- t( mse_rt.p100_node )
mse_rt.p100_node         <- as.data.frame( mse_rt.p100_node )

mse_rt.p1000_node[ 1, ]  <- round( colMeans( mse_rt.p1000 ), 6 )
mse_rt.p1000_node        <- t( mse_rt.p1000_node )
mse_rt.p1000_node        <- as.data.frame( mse_rt.p1000_node )

mse_rt.p10000_node[ 1, ] <- round( colMeans( mse_rt.p10000 ), 6 )
mse_rt.p10000_node       <- t( mse_rt.p10000_node )
mse_rt.p10000_node       <- as.data.frame( mse_rt.p10000_node )
# -- Relative MSE
relmse_rt.p100_node[ 1, ]   <- t(mse_rt.p100_node)/trueRate
relmse_rt.p100_node         <- t(relmse_rt.p100_node)
relmse_rt.p100_node         <- as.data.frame(100 * relmse_rt.p100_node)

relmse_rt.p1000_node[ 1, ]  <- t(mse_rt.p1000_node)/trueRate
relmse_rt.p1000_node        <- t(relmse_rt.p1000_node)
relmse_rt.p1000_node        <- as.data.frame(100 * relmse_rt.p1000_node)

relmse_rt.p10000_node[ 1, ] <- t(mse_rt.p10000_node)/trueRate
relmse_rt.p10000_node       <- t(relmse_rt.p10000_node)
relmse_rt.p10000_node       <- as.data.frame(100 * relmse_rt.p10000_node)


# -- Bias
bias_rt.p100_node[ 1, ]   <- round( colMeans( bias_rt.p100 ), 6 )
bias_rt.p100_node         <- t( bias_rt.p100_node )
bias_rt.p100_node         <- as.data.frame( bias_rt.p100_node )

bias_rt.p1000_node[ 1, ]  <- round( colMeans( bias_rt.p1000 ), 6 )
bias_rt.p1000_node        <- t( bias_rt.p1000_node )
bias_rt.p1000_node        <- as.data.frame( bias_rt.p1000_node )

bias_rt.p10000_node[ 1, ] <- round( colMeans( bias_rt.p10000 ), 6 )
bias_rt.p10000_node       <- t( bias_rt.p10000_node )
bias_rt.p10000_node       <- as.data.frame( bias_rt.p10000_node )
# -- Relative Bias
relbias_rt.p100_node[ 1, ]   <- t(bias_rt.p100_node)/trueRate
relbias_rt.p100_node         <- t(relbias_rt.p100_node)
relbias_rt.p100_node         <- as.data.frame(100 * relbias_rt.p100_node)

relbias_rt.p1000_node[ 1, ]  <- t(bias_rt.p1000_node)/trueRate
relbias_rt.p1000_node        <- t(relbias_rt.p1000_node)
relbias_rt.p1000_node        <- as.data.frame(100 * relbias_rt.p1000_node)

relbias_rt.p10000_node[ 1, ] <- t(bias_rt.p10000_node)/trueRate
relbias_rt.p10000_node       <- t(relbias_rt.p10000_node)
relbias_rt.p10000_node       <- as.data.frame(100 * relbias_rt.p10000_node)


# -- Width
width_rt.p100_node[ 1, ]   <- round( colMeans( width_rt.p100 ), 6 )
width_rt.p100_node         <- t( width_rt.p100_node )
width_rt.p100_node         <- as.data.frame( width_rt.p100_node )

width_rt.p1000_node[ 1, ]  <- round( colMeans( width_rt.p1000 ), 6 )
width_rt.p1000_node        <- t( width_rt.p1000_node )
width_rt.p1000_node        <- as.data.frame( width_rt.p1000_node )

width_rt.p10000_node[ 1, ] <- round( colMeans( width_rt.p10000 ), 6 )
width_rt.p10000_node       <- t( width_rt.p10000_node )
width_rt.p10000_node       <- as.data.frame( width_rt.p10000_node )
# -- Relative Width
relwidth_rt.p100_node[ 1, ]   <- t(width_rt.p100_node)/trueRate
relwidth_rt.p100_node         <- t(relwidth_rt.p100_node)
relwidth_rt.p100_node         <- as.data.frame(100 * relwidth_rt.p100_node)

relwidth_rt.p1000_node[ 1, ]  <- t(width_rt.p1000_node)/trueRate
relwidth_rt.p1000_node        <- t(relwidth_rt.p1000_node)
relwidth_rt.p1000_node        <- as.data.frame(100 * relwidth_rt.p1000_node)

relwidth_rt.p10000_node[ 1, ] <- t(width_rt.p10000_node)/trueRate
relwidth_rt.p10000_node       <- t(relwidth_rt.p10000_node)
relwidth_rt.p10000_node       <- as.data.frame(100 * relwidth_rt.p10000_node)


## -------------------------------- ##
## ----  SUMMARISE RESULTS     ---- ##
## -------------------------------- ##

# MSE
mse_node <- cbind( mse_rt.p100_node,       relmse_rt.p100_node,
                   mse_rt.p1000_node,      relmse_rt.p1000_node,
                   mse_rt.p10000_node,     relmse_rt.p10000_node
                   )
colnames(mse_node) <- c( "rt100-lyx",   "rel-rt100-lyx",
                         "rt1000-lyx",  "rel-rt1000-lyx",
                         "rt10000-lyx", "rel-rt10000-lyx"
                         )

summ.mse <- matrix( 0, nrow = 3, ncol = 1 )
summ.mse[1,1] <- paste( "\"", round( mse_rt.p100_node, 3 ), " (", 
                        round( relmse_rt.p100_node, 2 ), ")",  
                        "\"", sep="" )
summ.mse[2,1] <- paste( "\"", round( mse_rt.p1000_node, 3 ), " (", 
                        round( relmse_rt.p1000_node, 2 ), ")",  
                        "\"", sep="" )
summ.mse[3,1] <- paste( "\"", round( mse_rt.p10000_node, 3 ), " (", 
                        round( relmse_rt.p10000_node, 2 ), ")",  
                        "\"", sep="" )
rownames( summ.mse ) <- c( "r.100", "r.1000", "r.10000" )

# Bias
bias_node <- cbind( bias_rt.p100_node,       relbias_rt.p100_node,
                    bias_rt.p1000_node,      relbias_rt.p1000_node,
                    bias_rt.p10000_node,     relbias_rt.p10000_node
                    )
colnames(bias_node) <- c( "rt100",   "rel-rt100",
                          "rt1000",  "rel-rt1000",
                          "rt10000", "rel-rt10000"
                          )

summ.bias <- matrix( 0, nrow = 3, ncol = 1 )
summ.bias[1,1] <- paste( "\"", round( bias_rt.p100_node, 3 ), " (", 
                          round( relbias_rt.p100_node, 2 ), ")",  
                          "\"", sep="" )
summ.bias[2,1] <- paste( "\"", round( bias_rt.p1000_node, 3 ), " (", 
                          round( relbias_rt.p1000_node, 2 ), ")",  
                          "\"", sep="" )
summ.bias[3,1] <- paste( "\"", round( bias_rt.p10000_node, 3 ), " (", 
                           round( relbias_rt.p10000_node, 2 ), ")",  
                           "\"", sep="" )
rownames( summ.bias ) <- c( "r.100", "r.1000", "r.10000" )

# Width
width_node <- cbind( width_rt.p100_node,       relwidth_rt.p100_node,
                     width_rt.p1000_node,      relwidth_rt.p1000_node,
                     width_rt.p10000_node,     relwidth_rt.p10000_node
                     )
colnames(width_node) <- c( "rt100",   "rel-rt100",
                           "rt1000",  "rel-rt1000",
                           "rt10000", "rel-rt10000"
                           )


summ.width <- matrix( 0, nrow = 3, ncol = 1 )
summ.width[1,1] <- paste( "\"", round( width_rt.p100_node, 3 ), " (",
                          round( relwidth_rt.p100_node, 2 ), ")",
                          "\"", sep="" )
summ.width[2,1] <- paste( "\"", round( width_rt.p1000_node, 3 ), " (",
                          round( relwidth_rt.p1000_node, 2 ), ")",
                          "\"", sep="" )
summ.width[3,1] <- paste( "\"", round( width_rt.p10000_node, 3 ), " (",
                          round( relwidth_rt.p10000_node, 2 ), ")",
                          "\"", sep="" )
rownames( summ.width ) <- c( "r.100", "r.1000", "r.10000" )

# Write in csv files the results
write.csv( summ.mse, paste( wd, "MSE_node_rate.csv", sep='' ), quote = F  )
write.csv( summ.bias, paste( wd, "Bias_node_rate.csv", sep='' ), quote = F  )
write.csv( summ.width, paste( wd, "Width_node_rate.csv", sep='' ), quote = F  )



