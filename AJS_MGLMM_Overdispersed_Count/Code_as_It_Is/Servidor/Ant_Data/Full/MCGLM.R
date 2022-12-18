# https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multivariate-linear-models.html
# Package --------------------------------------------------------------------
rm(list = ls())
library(mcglm)
library(TMB)
library(mvabund) #Ant data
setwd("/home/guilherme/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/Ant_Data")
source("/home/guilherme/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/unify_packages.R")
# Dataset --------------------------------------------------------------------
data("antTraits")
##################################
# MCGLM
##################################
y <- antTraits$abund # Selecting response variables
sapply(y, function(x) rbind(var(x), mean(x)))
nr <- ncol(y)
names(y) <- paste0("y",1:nr)
X <- antTraits$env # Selecting covariates
data <- data.frame(y, X)
# Linear predictor
lp <- paste(names(y), "~", paste0(colnames(X), collapse = " + "))
form <- lapply(lp, formula)
# Matrix linear predictor
Z0 <- mc_id(data)
# Fitting multivariate model
fit1 <- mcglm(linear_pred = c(form), matrix_pred = rep(list(Z0), nr),
              link = rep("log", nr), variance = rep("poisson_tweedie", nr),
              control_algorithm = list(max_iter = 100),
              data = data)
res_cov <- as.matrix(residuals(fit1))
gc()
save.image("Poisson_MCGLM.RData")
