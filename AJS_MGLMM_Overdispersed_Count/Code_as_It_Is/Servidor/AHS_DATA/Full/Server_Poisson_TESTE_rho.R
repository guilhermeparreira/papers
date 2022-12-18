# https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multivariate-linear-models.html
# Package --------------------------------------------------------------------
rm(list = ls())
library(mcglm)
library(TMB)
setwd("/home/guilherme/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor")
source("unify_packages.R")
# Dataset --------------------------------------------------------------------
# Sample
## Full
set.seed(2390)
ss <- sample(1:nrow(ahs), 1000, replace = F)
ahs <- ahs[ss[-1], ] # COMMENT LATER
ahs$id <- 1:nrow(ahs)
n <- nrow(ahs)
# General variables
vars_resp <- c("Ndoc", "Nndoc", "Nmed", "Nhosp", "Nadm")
n_resp <- length(vars_resp)
n_params <- 11
units <- 60*60 # Second to hours
threads <- 6
## Fixed Effects -------------------------------------------------------------
form_Ndoc <- Ndoc ~ sex + age + income + levyplus + freepoor + freerepa + illness + actdays + chcond
# All Data sets for TMB --------------------------------------------------------------------
## Poisson
data_P <- list(Y = as.matrix(ahs[, vars_resp]),
               X = model.matrix(form_Ndoc, ahs))
# Start Parameters ------------------------------------------------------------------------
U = matrix(0, ncol = n_resp, nrow = n)
rho_start = rep(0, (n_resp*(n_resp-1))/2)
######################
# Poisson Multivariate ------------------------------------------------------------------------
######################
model <- "02_poisson_multivariate" #1st try
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
#### 300 sample
parameters <- list(beta = matrix(0, nrow = n_params, ncol = n_resp),
                   U = U,
                   rho = rho_start,
                   sigma = rep(1, n_resp))
openmp(n=threads)
obj <- MakeADFun(data = data_P,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = F,
                 random = "U")
gc()
a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(eval.max = 7000, iter.max = 7000,
                                                  abs.tol = 1e-04, rel.tol = 1e-04)))
theta <- fit_NLB$par[names(fit_NLB$par)=="rho"]
sigma <- fit_NLB$par[names(fit_NLB$par)=="sigma"]
cria_cor_matrix_tmb_sd_from_vector(n_resp, theta)$matrix_corre
getSigma(theta, sigma, varcov = F)

rep <- sdreport(obj)
rep_fe <- summary(rep, "fixed")
rep_cor <- summary(rep, "report");rep_cor
estimates_poisson_300 <- TMB_summary(rep_fe, rep_cor, 5)
estimates_poisson_300$cor_mat