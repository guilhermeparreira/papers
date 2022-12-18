# Poisson Multivariate ------------------------------------------------------------------------
######################
rm(list = ls())
# # Packages
library(mvabund) #Ant data
data("antTraits")
library(TMB)
# # External Codes
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "Ant_Data/SmallerModels"))
source(paste0(mainpath, "unify_packages.R"))
# .cpp File
model <- "02_poisson_multivariate" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
# Data
y <- antTraits$abund # Selecting response variables
X <- antTraits$env # Selecting covariates
nr <- ncol(y)
data <- list(X = as.matrix(cbind(1, X)),
             Y = as.matrix(y))
n_params <- 6
# Parameters
U <- matrix(0, nrow = nrow(y), ncol = ncol(y))
parameters <- list(beta = matrix(0, nrow = ncol(data$X), ncol = nr),
                   U = U,
                   rho = rep(0, (nr*(nr-1))/2),
                   sigma = sapply(y, sd))
gc()
threads <- 6
openmp(n=threads)
