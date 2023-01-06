##########################################################################################
#           ***** ROUTINES *****
##########################################################################################

#--------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
# Kullback Leibler divergence 
#-----------------------------------------------------------------------------------------------
# par: parameters of the PGW distribution
# initW: initial values for the optimisation step (finding closest Weibull distribution)

# Lower integration limit
UL <- 0.001

minKL <- function(par,initW, method = "optim"){
  if(par[3]==1) out <- list(parW = par[1:2], min.KL = 0)
  if(par[3] != 1){
  initW <- log(initW)
  # Distance function in terms of the parameters of the Weibull dist.
  tempfW <- function(parW){
    # W parameters
    sigma0 <- exp(parW[1])
    nu0 <- exp(parW[2])
    # Distance between W and PGW
    integrand <- Vectorize(function(t) dpgw(t, par[1], par[2], par[3])*( 
      dpgw(t, par[1], par[2], par[3], log = T) - dweibull(t, scale = sigma0, shape = nu0, log = T) )  )
    val <- integrate(integrand, UL, Inf, subdivisions = 10000)$value # Consider controlling the error
    return(val)
  }
  # Distance minimisation
  if(method == "optim"){
  OPT <- optim(initW, tempfW, control = list(maxit = 10000))
  out <- list(parW = exp(OPT$par), min.KL = OPT$value)
  }
  if(method == "nlminb"){
    OPT <- nlminb(initW, tempfW, control = list(maxit = 10000))
    out <- list(parW = exp(OPT$par), min.KL = OPT$objective)
  }
  return(out)
  }
}


########################################################################################
# Minimum Hellinger between the PGW distribution and the Weibull family
# Returns the object from optim 
########################################################################################
# par: parameters of the PGW distribution
# initW: initial values for the optimisation step (finding closest Weibull distribution)

#-----------------------------------------------------------------------------------------------
# Minimum Hellinger distance 
#-----------------------------------------------------------------------------------------------

minHell <- function(par, initW){
  if(par[3]==1) out <- list(parW = par[1:2], min.dHell = 0)
  if(par[3] != 1){
    initW <- log(initW)
    # Distance function in terms of the parameters of the Weibull dist.
    tempfW <- function(parW){
      # W parameters
      sigma0 <- exp(parW[1])
      nu0 <- exp(parW[2])
      # Distance between W and PGW
      integrand <- Vectorize(function(t) ( sqrt(dpgw(t, par[1], par[2], par[3]))  - 
                                             sqrt(dweibull(t, sigma0, nu0) )  )^2)
      val <- sqrt( 0.5*integrate(integrand, 0, Inf, subdivisions = 10000)$value ) # Consider controlling the error
      return(val)
    }
    # Distance minimisation
    OPT <- optim(initW, tempfW, control = list(maxit = 10000))
    out <- list(parW = exp(OPT$par), min.dHell = OPT$value)
    return(out)
  }
}
