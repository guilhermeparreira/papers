# https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multivariate-linear-models.html
# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "AHS_DATA/Full"))
#### 5190 full
# I took the sample data from nlminb. And I wat to adjust the MCMC
load("Poisson_sample.RData")
model <- "02_poisson_multivariate_sem_paralelizar"
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
rho.z <- (1/2)*log((1+estimates_poisson_am$rho)/(1-estimates_poisson_am$rho))
(exp(2*rho.z)-1)/(exp(2*rho.z)+1)
parameters <- list(beta = estimates_poisson_am$beta,
                   U = U,
                   rho = rho.z,
                   sigma = log(estimates_poisson_am$sigma))
dyn.load(dynlib(paste0(mainpath, model)))
FreeADFun(obj)
# 17986.51 # COrrect value
# 17436.4 # Primeiro
obj <- MakeADFun(data = data_am,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
obj$fn()

# DEU CERTO!!!!!!! Então, fechou, está correto o modelo! Falta limpar o código

model <- "02_poisson_multivariate"
parameters <- list(beta = estimates_poisson_am$beta,
                   U = U,
                   rho = estimates_poisson_am$rho,
                   sigma = estimates_poisson_am$sigma)
dyn.load(dynlib(paste0(mainpath, model)))
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
openmp(n=6)
fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 1e7, iter.max = 1e7,
                                 abs.tol = 1e-04, rel.tol = 1e-04))

gc()
obj$fn()


# ?rstan::sampling
rstan_options(auto_write = TRUE)
fit <- tmbstan::tmbstan(obj, silent = F, chains = 5, iter = 8000, cores = 5)
save.image("Poisson_stan.RData")
# a <- system.time(fit_NLB <- nlminb(obj$par, obj$fn, obj$gr,
#                                    control = list(eval.max = 7000, iter.max = 7000,
#                                                   abs.tol = 1e-04, rel.tol = 1e-04)))
# rep <- sdreport(obj)
# rep_fe <- summary(rep, "fixed")
# rep_cor <- summary(rep, "report")
# estimates_poisson <- TMB_summary(rep_fe, rep_cor, 5)
# estimates_poisson <- cbind(estimates_poisson, a[1]/units, fit_NLB$convergence, fit_NLB$objective, threads)
# start_beta <- matrix(estimates_poisson[grep("b\\d+_Y", rownames(estimates_poisson)), 1],ncol = n_resp)
# colnames(start_beta) <- paste0("beta", 1:5)
# start_sigma <- estimates_poisson[grep("sig", rownames(estimates_poisson)), 1]
# gc()
# save.image("Poisson_full.RData")
