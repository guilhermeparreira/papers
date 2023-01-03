rm(list = ls())
library(TMB)
library(mcglm)
library(gmailr)
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
complement <- "NHANES_DATA/SmallerModels"
pathlocal <- paste0(mainpath, complement)
setwd(paste0(pathlocal))
source(paste0(mainpath, "unify_packages.R"))
# Dataset --------------------------------------------------------------------
nhanes <- read.table(paste0(mainpath, "NHANES_DATA/nhanes.txt"), header = T)
# Sample
set.seed(2390)
nam <- 350
ss <- sample(1:nrow(nhanes), nam, replace = F)
nhanesam <- nhanes[ss, ]
## Full
n <- nrow(nhanes)
# General variables
vars_resp <- c("Nmsp", "Nmosp", "Nspfy")
nr <- length(vars_resp)
units <- 60*60 # Second to hours
threads <- 6
n_params <- 5
## Fixed Effects -------------------------------------------------------------
form_Nmsp <- Nmsp ~ Race + Education + Marital + Age
form_Nmosp <- Nmosp ~ Race + Education + Marital + Age
form_Nspfy <- Nspfy ~ Race + Education + Marital + Age
# All Data sets for TMB --------------------------------------------------------------------
data_am <- list(Y = as.matrix(nhanesam[, vars_resp]),
                   X = model.matrix(form_Nmsp, nhanesam))
data <- list(Y = as.matrix(nhanes[, vars_resp]),
               X = model.matrix(form_Nmsp, nhanes))
# Start Parameters ------------------------------------------------------------------------
U_am = matrix(0, ncol = nr, nrow = nam)
U = matrix(0, ncol = nr, nrow = n)
rho_start = rep(0, (nr*(nr-1))/2)
######
## Start Model - MCGLM -------------------------------------------------------------
######
Z0 <- mc_id(nhanes)
system.time(fit_refe <- mcglm(linear_pred = c(form_Nmsp, form_Nmosp, form_Nspfy),
                              matrix_pred = list(Z0,Z0,Z0),
                              link = rep("log", 3), variance = rep("tweedie", 3),
                              power_fixed = rep(TRUE, 3), data = nhanes,
                              control_algorithm = list(correct = FALSE, max_iter = 100)))
gc()
save.image("nhanes_mcglm.RData")
