#==============================================================#
# Kostas' Code                                                 #
#==============================================================#
# This is Kostas Angelis' function and it is used to calculate
# the coverage.
# It was first used to summarise the results obtained 
# in Angelis et al., 2017, https://doi.org/10.1093/sysbio/syx061

xInRange <- function(x, range){
# Returns 1 if x is within the 
# range of values, else 0
  n <- length( x )
  r <- numeric( n )
  for( i in 1:n ){
	if( x[i] > range[i,1] && x[i] < range[i,2] ) r[i] <- TRUE
	else r[i] <- FALSE
  }
  return( r )
}

