## Negative Binomial distribution ----------------------------------------------
## Authors: Multivariate Statistical Modelling Group ---------------------------
## Date: June, 04 2020 ---------------------------------------------------------
rm(list=ls())
require(MASS)
require(mvtnorm)
require(TMB)
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
source("/home/guilherme/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/unify_packages.R")
## Compilation
model <- "03_negbin_bivariate_wag"
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
## Data simulation ------------------------------------------------------------
beta1 <- c(2, .8)
beta2 <- c(.5, -1)
s2_1 <- .5
s2_2 <- .3
true_rho <- 0
phi0 <- 1
phi1 <- 1
n <- 200
n_pars <- length(c(beta1, beta2, s2_1, s2_2, true_rho, phi0, phi1))
True <- c(beta1,beta2,true_rho,log(s2_1),log(s2_2),log(phi0),log(phi1))

Escala_True <- c(beta1 = beta1,beta2 = beta2,rho = true_rho,s2_1=s2_1,s2_2=s2_2,phi0=phi0,phi1=phi1)
## Change --------------------------------------------------------------------
# Size 100
sample_size <- 100
repeticoes <- 50
lists100 <- list()
for (i in 1:length(sample_size)) {
  for (j in 1:repeticoes){
  n <- sample_size[i]
  data <- rnegbinre(beta1, beta2, n, s2_1, s2_2, true_rho, phi0, phi1, 1234+j)
  parameters <- list(beta1 = beta1, beta2 = beta2,
                     U = matrix(0, ncol = 2, nrow = n),
                     rho = true_rho, sigma = c(0.1,0.1), phi = c(log(phi0), log(phi1)))
  obj <- MakeADFun(data, parameters, DLL=model, random = "U", silent = T)
  fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
  rep <- sdreport(obj)
  fixef <- summary(rep, "fixed")
  cor <- summary(rep, "report")
  pars_used <- TMB_summary_wag(fixef, cor, 2)
  lists100[[j]] <- pars_used
  print(lists100)
  }
}
# Size 300
sample_size <- 300
repeticoes <- 50
lists300 <- list()
for (i in 1:length(sample_size)) {
  for (j in 1:repeticoes){
    n <- sample_size[i]
    data <- rnegbinre(beta1, beta2, n, s2_1, s2_2, true_rho, phi0, phi1, 1234+j)
    parameters <- list(beta1 = beta1, beta2 = beta2,
                       U = matrix(0, ncol = 2, nrow = n),
                       rho = true_rho, sigma = c(0.1,0.1), phi = c(log(phi0), log(phi1)))
    obj <- MakeADFun(data, parameters, DLL=model, random = "U", silent = T)
    fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
    rep <- sdreport(obj)
    fixef <- summary(rep, "fixed")
    cor <- summary(rep, "report")
    # Joining all together
    pars_used <- TMB_summary_wag(fixef, cor, 2)
    lists300[[j]] <- pars_used
    print(lists300)
  }
}
# Size 500
sample_size <- 500
repeticoes <- 50
lists500 <- list()
for (i in 1:length(sample_size)) {
  for (j in 1:repeticoes){
    n <- sample_size[i]
    data <- rnegbinre(beta1, beta2, n, s2_1, s2_2, true_rho, phi0, phi1, 1234+j)
    parameters <- list(beta1 = beta1, beta2 = beta2,
                       U = matrix(0, ncol = 2, nrow = n),
                       rho = true_rho, sigma = c(0.1,0.1), phi = c(log(phi0), log(phi1)))
    obj <- MakeADFun(data, parameters, DLL=model, random = "U", silent = T)
    fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
    rep <- sdreport(obj)
    fixef <- summary(rep, "fixed")
    cor <- summary(rep, "report")
    # Joining all together
    pars_used <- TMB_summary_wag(fixef, cor, 2)
    lists500[[j]] <- pars_used
    print(lists500)
  }
}
# Size 1000
sample_size <- 1000
repeticoes <- 50
lists1000 <- list()
for (i in 1:length(sample_size)) {
  for (j in 1:repeticoes){
    n <- sample_size[i]
    data <- rnegbinre(beta1, beta2, n, s2_1, s2_2, true_rho, phi0, phi1, 1234+j)
    parameters <- list(beta1 = beta1, beta2 = beta2,
                       U = matrix(0, ncol = 2, nrow = n),
                       rho = true_rho, sigma = c(0.1,0.1), phi = c(log(phi0), log(phi1)))
    obj <- MakeADFun(data, parameters, DLL=model, random = "U", silent = T)
    fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
    rep <- sdreport(obj)
    fixef <- summary(rep, "fixed")
    cor <- summary(rep, "report")
    # Joining all together
    pars_used <- TMB_summary_wag(fixef, cor, 2)
    lists1000[[j]] <- pars_used
    # print(lists1000)
  }
}
# Size 5000
sample_size <- 5000
repeticoes <- 50
lists5000 <- list()
for (i in 1:length(sample_size)) {
  for (j in 1:repeticoes){
    n <- sample_size[i]
    data <- rnegbinre(beta1, beta2, n, s2_1, s2_2, true_rho, phi0, phi1, 1234+j)
    parameters <- list(beta1 = beta1, beta2 = beta2,
                       U = matrix(0, ncol = 2, nrow = n),
                       rho = true_rho, sigma = c(0.1,0.1), phi = c(log(phi0), log(phi1)))
    obj <- MakeADFun(data, parameters, DLL=model, random = "U", silent = T)
    fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
    rep <- sdreport(obj)
    fixef <- summary(rep, "fixed")
    cor <- summary(rep, "report")
    # Joining all together
    pars_used <- TMB_summary_wag(fixef, cor, 2)
    lists5000[[j]] <- pars_used
    # print(lists5000)
  }
}
sample_sizes <- c(100, 300, 500, 1000,5000)
saveRDS(list(lists100, lists300, lists500, lists1000, lists5000), "estudo_simulacao.rds")
all_lists <- list(lists100, lists300, lists500, lists1000, lists5000)
# 100
df100 <- plyr::ldply(lists100, rbind.data.frame)
df100$Par <- c("beta1_0", "beta1_1", "beta2_0", "beta2_1", "s2_1", "s2_2", "phi1", "phi2", "cor")
df100$Int_Inf <- df100$Estimate-1.96*df100$`Std. Error`
df100$Int_Sup <- df100$Estimate+1.96*df100$`Std. Error`
df100$n <- "100"
df100$ordem <- rep(1:50, each = 9)

ggplot(df100, aes(x = ordem, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Int_Inf, ymax = Int_Sup)) +
  facet_wrap(~Par, scales = "free") +
  coord_flip()


#300
df300 <- plyr::ldply(lists300, rbind.data.frame)
df300$Par <- c("beta1", "beta1", "beta2", "beta2", "s2_1", "s2_2", "phi1", "phi2", "cor")
df300$Int_Inf <- df300$Estimate-1.96*df300$`Std. Error`
df300$Int_Sup <- df300$Estimate+1.96*df300$`Std. Error`
df300$n <- "300"
#500
df500 <- plyr::ldply(lists500, rbind.data.frame)
df500$Par <- c("beta1", "beta1", "beta2", "beta2", "s2_1", "s2_2", "phi1", "phi2", "cor")
df500$Int_Inf <- df500$Estimate-1.96*df500$`Std. Error`
df500$Int_Sup <- df500$Estimate+1.96*df500$`Std. Error`
df500$n <- "500"
#1000
df1000 <- plyr::ldply(lists1000, rbind.data.frame)
df1000$Par <- c("beta1", "beta1", "beta2", "beta2", "s2_1", "s2_2", "phi1", "phi2", "cor")
df1000$Int_Inf <- df1000$Estimate-1.96*df1000$`Std. Error`
df1000$Int_Sup <- df1000$Estimate+1.96*df1000$`Std. Error`
df1000$n <- "1000"
#5000
df5000 <- plyr::ldply(lists5000, rbind.data.frame)
df5000$Par <- c("beta1", "beta1", "beta2", "beta2", "s2_1", "s2_2", "phi1", "phi2", "cor")
df5000$Int_Inf <- df5000$Estimate-1.96*df5000$`Std. Error`
df5000$Int_Sup <- df5000$Estimate+1.96*df5000$`Std. Error`
df5000$n <- "5000"
dfs <- rbind(df100, df300, df500, df1000, df5000)
dfs$ordem <- c(rep(rep(1:50, each = 9), 4), rep(1:46, each = 9))
saveRDS(dfs, "estudo_simulacao_em_df.rds")
library(ggplot2)
readRDS("estudo_simulacao_em_df.rds")
ggplot(dfs[dfs$Par=="cor", ], aes(x = ordem, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Int_Inf, ymax = Int_Sup)) +
  facet_wrap(~n) +
  scale_y_continuous(limits = c(-1,1)) +
  labs(title = "Valores dos coeficientes de correlação estimados nas 50 repetições de uma Binomial Negative",
       subtitle = "Verdadeiro valor da correlação = 0")
  
