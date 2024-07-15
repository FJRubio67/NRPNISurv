## ----------------------------------------------------------------------------------------------------------------------------
########################################################################
# Case study II: Leukemia data set from the spBayesSurv R package
########################################################################

# See README file for more information
rm(list = ls())

# Routines
# https://github.com/FJRubio67/NRPNISurv
source("routines.R")

# Required packages (from CRAN)
# install.packages('PackageName')
library(knitr)
library(survival)
library(spBayesSurv)
library(numDeriv)

# Required packages (from GitHub)
#library(devtools)
#install_github("FJRubio67/HazReg")
library(HazReg)


#----------------------------------------------------------------------------------------------
# Leukemia data
#---------------------------------------------------------------------------------------------

data(LeukSurv)
?LeukSurv
head(LeukSurv)
dim(LeukSurv)

#******************
# Data preparation
#******************

n <- dim(LeukSurv)[1]  # number of individuals
# Design matrices
x <- as.matrix(cbind(scale(LeukSurv$age), LeukSurv$sex, scale(LeukSurv$wbc), scale(LeukSurv$tpi) ))
colnames(x) <- cbind("std age", "sex", "wbc", "tpi")
xt <- as.matrix(cbind(scale(LeukSurv$age), scale(LeukSurv$wbc), scale(LeukSurv$tpi)))

# Required quantities
status <- as.vector(LeukSurv$cens)
times <- as.vector(LeukSurv$time)/365.24 # in years


## ----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
# Maximum likelihood estimation
# See : https://github.com/FJRubio67/HazReg
#----------------------------------------------------------------------------------------------

# PGW model with no covariates
OPTPGW0 <- GHMLE(init = c(0,0,0), times = times, status = status, 
                  hstr = "baseline", dist = "PGW", method = "nlminb", maxit = 10000)

# PGW-PH model
OPTPGWPH <- GHMLE(init = c(0,0,0, rep(0,ncol(x))), times = times, status = status, 
     hstr = "PH", dist = "PGW", des = x, method = "nlminb", maxit = 10000)

# PGW-AFT model
OPTPGWAFT <- GHMLE(init = c(0,0,0,rep(0,ncol(x))), times = times, status = status, 
        hstr = "AFT", dist = "PGW", des = x, method = "nlminb", maxit = 10000)

# PGW-GH model
OPTPGWGH <- GHMLE(init = c(0,0,0, rep(0, ncol(xt)), rep(0,ncol(x))), times = times, status = status, 
            hstr = "GH", dist = "PGW", des = x, des_t = xt, method = "nlminb", maxit = 10000)

# MLE PGW no covariates
MLE.PGW0 <- exp(OPTPGW0$OPT$par[1:3])
kable(MLE.PGW0, digits = 3)

# MLE PGW-PH structure
MLE.PGWPH <- c(exp(OPTPGWPH$OPT$par[1:3]), OPTPGWPH$OPT$par[-c(1:3)])
kable(MLE.PGWPH, digits = 3)

# MLE PGW-AFT structure
MLE.PGWAFT <- c(exp(OPTPGWAFT$OPT$par[1:3]), OPTPGWAFT$OPT$par[-c(1:3)])
kable(MLE.PGWAFT, digits = 3)

# MLE PGW-GH structure
MLE.PGWGH <- c(exp(OPTPGWGH$OPT$par[1:3]), OPTPGWGH$OPT$par[-c(1:3)])
kable(MLE.PGWGH, digits = 3)


## ----------------------------------------------------------------------------------------------------------------------------
# AIC
AIC.PGW0 <- 2*OPTPGW0$OPT$objective + 2*length(MLE.PGW0)
AIC.PGWGH <- 2*OPTPGWGH$OPT$objective + 2*length(MLE.PGWGH)
AIC.PGWAFT <- 2*OPTPGWAFT$OPT$objective + 2*length(MLE.PGWAFT)
AIC.PGWPH <- 2*OPTPGWPH$OPT$objective + 2*length(MLE.PGWPH)

# Best model: PGW-GH
c(AIC.PGW0, AIC.PGWPH, AIC.PGWAFT, AIC.PGWGH)


## ----------------------------------------------------------------------------------------------------------------------------
######################################################################
# KL Divergence criterion
######################################################################
# KL Divergence
dKL <- minKL(MLE.PGWGH[1:3], c(1,1))$min.KL
dKL

# Redundancy constant
k = length(MLE.PGWGH)
M = 0.05
ne <- n - 0.5*sum(!status)

UKL <- 0.5*k*M*length(MLE.PGWGH)*log(ne)/ne

# Test
(dKL < UKL)



######################################################################
# Hellinger distance criterion
######################################################################

# Hellinger Distance
DH <- minHell(MLE.PGWGH[1:3], c(1,1))$min.dHell
DH

# Hellinger criterion
kappa <- 0.05
UH <- log( 1 - (1-2*kappa)^2)/(2*log(1-DH^2))

# Test
(ne < UH)


######################################################################
# Hessian Method
######################################################################

# Hessian/Eigenvalues
HESS.PGW  <- -hessian(OPTPGWGH$log_lik, x = OPTPGWGH$OPT$par)
  
eigen.val <- eigen(HESS.PGW)$values

ref.val <- abs(as.vector(eigen.val))/max(abs(as.vector(eigen.val)))

sev <- sort(ref.val)

sev

# Upper bound
Un <- 0.001

# Test
(sev[1] < Un)


## ----------------------------------------------------------------------------------------------------------------------------
#######################################################################################
# Bootstrap
#######################################################################################

# Number of Boostrap samples 
B <- 1000

# Bootstrap MLEs
MLE.B <- matrix(0, ncol = length(MLE.PGWGH), nrow = B)
neB <- UKLB <- UHB <- DHB<- DKLB <- vector()
indHess <- vector()

for(i in 1:B){
  ind <- sample(1:n, replace = T)
  OPTB <- GHMLE(init = c(0,0,0, rep(0, ncol(xt)), rep(0,ncol(x))), times = times[ind], status = status[ind], 
         hstr = "GH", dist = "PGW", des = x[ind,], des_t = xt[ind,], method = "nlminb", maxit = 10000)
  MLE.B[i,] <- c(exp(OPTB$OPT$par[1:3]), OPTB$OPT$par[-c(1:3)])
  neB[i] <- n - 0.5*sum(!status[ind])
  
  UKLB[i] <- 0.5*k*M*length(MLE.B[i,])*log(neB[i])/neB[i]
  
  # Hessian/Eigenvalues
  HESSB  <- -hessian(OPTB$log_lik, x = OPTB$OPT$par)
  
  eigen.val <- eigen(HESSB)$values
  
  ref.val <- abs(as.vector(eigen.val))/max(abs(as.vector(eigen.val)))
  
  sev <- sort(ref.val)
  
  # Test
  indHess[i] <- as.numeric(sev[1] < Un)
}


######################################################################
# Bootstrap KL Divergence criterion
######################################################################

# Bootstrap KL Divergence

for(i in 1:B){
  DKLB[i] <- minKL(MLE.B[i,1:3], c(1,1))$min.KL
}

# Bootstrap probability of KL divergence criterion

indKLB <- (DKLB < UKLB) 

mean(indKLB)




######################################################################
# Bootstrap Hellinger distance criterion
######################################################################

# Hellinger Distance
for(i in 1:B){
DHB[i] <- minHell(MLE.B[i,1:3], c(1,1))$min.dHell

# Hellinger criterion
UHB[i] <- log( 1 - (1-2*kappa)^2)/(2*log(1-DHB[i]^2))
}

# Bootstrap probability of Hellinger distance criterion
indHB <- (neB < UHB)

mean(indHB)

######################################################################
# Bootstrap Hessian method
######################################################################


# Bootstrap probability of KL divergence criterion

mean(indHess) 


## ----------------------------------------------------------------------------------------------------------------------------
#######################################################################################
# Profile likelihood
#######################################################################################
p0 <- ncol(xt)
p1 <- ncol(x)

p <- length(MLE.PGWGH) # number of parameters
ML <- OPTPGWGH$OPT$objective # Maximum likelihood


# Profile likelihood function for parameter "ind"
prof.lik <- function(par1, ind){
  
  tempf <- function(par){
    tempv <- rep(0,p)
    tempv <- replace(x = tempv, c(1:p)[-ind] , par)
    tempv <- replace(x = tempv, ind , par1)
    out0 <- OPTPGWGH$log_lik(tempv)
    return(out0)
  } 
  
  out0 <-  -nlminb(OPTPGWGH$OPT$par[-ind],tempf, control = list(iter.max = 10000))$objective + ML
  out1 <-  -nlminb(OPTPGWGH$OPT$par[-ind],tempf, control = list(iter.max = 10000))$objective + ML
  out <- min(out0, out1)
  
  out2 <- ifelse(exp(out)<=1, exp(out), 0)
  return(out2)
}


# Profile likelihoods

# Profile likelihood of Parameter 1
prof1 <- Vectorize(function(par) prof.lik(log(par),1))
curve(prof1,0.05,0.2 , n = 200, lwd = 2, xlab = expression(sigma), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 2
prof2 <- Vectorize(function(par) prof.lik(log(par),2))
curve(prof2,0.8,1.3 , n = 200, lwd = 2, xlab = expression(nu), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 3
prof3 <- Vectorize(function(par) prof.lik(log(par),3))
curve(prof3,2,5 , n = 200, lwd = 2, xlab = expression(gamma), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 4
prof4 <- Vectorize(function(par) prof.lik(par,4))
curve(prof4,0.2,1.6, n = 200, lwd = 2, xlab = expression(alpha[1]), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 5
prof5 <- Vectorize(function(par) prof.lik(par,5))
curve(prof5,0.25, 1.75 , n = 200, lwd = 2, xlab = expression(alpha[2]), ylab = "Profile Likelihood")


# Profile likelihood of Parameter 6
prof6 <- Vectorize(function(par) prof.lik(par,6))
curve(prof6,-0.5, 1.25 , n = 200, lwd = 2, xlab = expression(alpha[3]), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 7
prof7 <- Vectorize(function(par) prof.lik(par,7))
curve(prof7,0.6,1.4, n = 200, lwd = 2, xlab = expression(beta[1]), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 8
prof8 <- Vectorize(function(par) prof.lik(par,8))
curve(prof8,-0.15, 0.35 , n = 200, lwd = 2, xlab = expression(beta[2]), ylab = "Profile Likelihood")


# Profile likelihood of Parameter 9
prof9 <- Vectorize(function(par) prof.lik(par,9))
curve(prof9,0.2, 1.2 , n = 200, lwd = 2, xlab = expression(beta[3]), ylab = "Profile Likelihood")


# Profile likelihood of Parameter 10
prof10 <- Vectorize(function(par) prof.lik(par,10))
curve(prof10,-0.2,0.8 , n = 200, lwd = 2, xlab = expression(beta[4]), ylab = "Profile Likelihood")


## ----setup, warning=FALSE----------------------------------------------------------------------------------------------------
##################################################################################
# Alternative models
##################################################################################


# Exponentiated Weibull

OPTEWGH <- GHMLE(init = c(0,0,0, rep(0, ncol(xt)), rep(0,ncol(x))), times = times, status = status, 
                 hstr = "GH", dist = "EW", des = x, des_t = xt, method = "nlminb", maxit = 10000)

# MLE EW-GH structure
MLE.EWGH <- c(exp(OPTEWGH$OPT$par[1:3]), OPTEWGH$OPT$par[-c(1:3)])
kable(MLE.EWGH, digits = 3)
AIC.EWGH <- 2*OPTEWGH$OPT$objective + 2*length(MLE.EWGH)
AIC.EWGH


#######################################################################################
# Profile likelihood
#######################################################################################


p2 <- length(MLE.EWGH) # number of parameters
ML2 <- OPTEWGH$OPT$objective # Maximum likelihood


# Profile likelihood function for parameter "ind"
prof.lik2 <- function(par1, ind){
  
  tempf <- function(par){
    tempv <- rep(0,p2)
    tempv <- replace(x = tempv, c(1:p2)[-ind] , par)
    tempv <- replace(x = tempv, ind , par1)
    out0 <- OPTEWGH$log_lik(tempv)
    return(out0)
  } 
  
  out <-  -nlminb(OPTEWGH$OPT$par[-ind],tempf, control = list(iter.max = 10000))$objective + ML2

  out2 <- ifelse(exp(out)<=1, exp(out), 0)
  return(out2)
}


# Profile likelihoods

# Profile likelihood of Parameter 1
prof21 <- Vectorize(function(par) prof.lik2(log(par),1))
curve(prof21,0.0001,0.075 , n = 200, lwd = 2, xlab = expression(sigma), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 2
prof22 <- Vectorize(function(par) prof.lik2(log(par),2))
curve(prof22,0.125,0.325 , n = 200, lwd = 2, xlab = expression(nu), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 3
prof23 <- Vectorize(function(par) prof.lik2(log(par),3))
curve(prof23,3.5,25 , n = 200, lwd = 2, xlab = expression(gamma), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 4
prof24 <- Vectorize(function(par) prof.lik2(par,4))
curve(prof24,0.2,1.6, n = 200, lwd = 2, xlab = expression(alpha[1]), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 5
prof25 <- Vectorize(function(par) prof.lik2(par,5))
curve(prof25,0.25, 1.5 , n = 200, lwd = 2, xlab = expression(alpha[2]), ylab = "Profile Likelihood")


# Profile likelihood of Parameter 6
prof26 <- Vectorize(function(par) prof.lik2(par,6))
curve(prof26,-0.25, 1 , n = 200, lwd = 2, xlab = expression(alpha[3]), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 7
prof27 <- Vectorize(function(par) prof.lik2(par,7))
curve(prof27,0.7,1.2, n = 200, lwd = 2, xlab = expression(beta[1]), ylab = "Profile Likelihood")

# Profile likelihood of Parameter 8
prof28 <- Vectorize(function(par) prof.lik2(par,8))
curve(prof28,-0.125, 0.3 , n = 200, lwd = 2, xlab = expression(beta[2]), ylab = "Profile Likelihood")


# Profile likelihood of Parameter 9
prof29 <- Vectorize(function(par) prof.lik2(par,9))
curve(prof29,0.4, 1.1 , n = 200, lwd = 2, xlab = expression(beta[3]), ylab = "Profile Likelihood")


# Profile likelihood of Parameter 10
prof210 <- Vectorize(function(par) prof.lik2(par,10))
curve(prof210,0,0.75 , n = 200, lwd = 2, xlab = expression(beta[4]), ylab = "Profile Likelihood")


