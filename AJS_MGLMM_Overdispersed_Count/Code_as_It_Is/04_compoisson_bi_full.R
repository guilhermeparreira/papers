## Compoisson distribution ----------------------------------------------
## Authors: Multivariate Statistical Modelling Group ---------------------------
## Date: June, 18 2020 ---------------------------------------------------------
rm(list=ls())
require(MASS)
require(mvtnorm)
require(TMB)
require(compoisson)
require(mcglm)
# Real Dataset --------------------------------------------------------------------
set.seed(2390)
n <- 200
ss <- sample(1:nrow(ahs), n, replace = F)
ahs <- ahs[ss, ]
ahs$id <- 1:nrow(ahs)
# Theta is the dispersion parameters
m1 <- MASS::glm.nb(Nmed ~ actdays, data = ahs)
m2 <- MASS::glm.nb(Ndoc ~ actdays, data = ahs)
phi1 <- summary(m1)$theta
phi2 <- summary(m2)$theta
l_phi1 <- log(phi1)
l_phi2 <- log(phi2)
data <- list(Y1 = ahs$Nadm,
             Y2 = ahs$Ndoc,
             X = model.matrix(~actdays, ahs))
n <- nrow(ahs)
parameters <- list(beta1 = c(log(mean(ahs$Nadm)), 0),
                   beta2 = c(log(mean(ahs$Ndoc)), 0),
                   U = matrix(0, ncol = 2, nrow = n),
                   rho = cor(ahs$Nadm, ahs$Ndoc),
                   sigma = c(.2,
                             .4),
                   nu = c(l_phi1, l_phi2))

## Data Import TMB ------------------------------------------------------------
n <- 200
## Data set
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
data <- readRDS("com_poisson_data_TMB.rds")
with(data, cor(Y1, Y2))
set.seed(45)
vecs <- sort(sample(1:n, n))
data$Y1 <- data$Y1[vecs]
data$Y2 <- data$Y2[vecs]
data$X <- data$X[vecs, ]

## Parameters to generate data
beta1 <- c(log(2), .7)
beta2 <- c(log(2), -1)
s1 <- log(.3) #.3 é desvio padrão
s2 <- log(.15) #.15 é desvio padrão
true_rho <- .5
l_nu1 <- log(1); l_nu2 <- log(1);
parameters <- list(beta1 = beta1,
                   beta2 = beta2,
                   U = matrix(0, ncol = 2, nrow = n),
                   rho = 0.5,
                   sigma = c(exp(s1), exp(s2)),
                   nu = c(l_nu1, l_nu2))
## Data Import R -------------------------------------------------------------
data <- readRDS("com_poisson_data_R.rds")
n <- 100
## Data set
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
with(data, cor(Y1, Y2))
set.seed(45)
vecs <- sort(sample(1:n, n))
data$Y1 <- data$Y1[vecs]
data$Y2 <- data$Y2[vecs]
data$X <- data$X[vecs, ]
beta1 <- c(log(.5), .7)
beta2 <- c(log(.5), -1)
s1 <- log(sqrt(.3025)) #.55 é desvio padrão
s2 <- log(sqrt(.09)) #.3 é desvio padrão
true_rho <- .5
l_nu1 <- log(1)
l_nu2 <- log(1);
parameters <- list(beta1 = beta1,
                   beta2 = beta2,
                   U = matrix(0, ncol = 2, nrow = n),
                   rho = 0.5,
                   sigma = c(exp(s1), exp(s2)),
                   nu = c(l_nu1, l_nu2))
## Data Import - Negative Binomial R ------------------------------------------
data <- readRDS("neg_binomial_data.rds")
parameters <- readRDS("neg_binomial_parameters.rds")
names(parameters)[6] <- "nu"
## Compilation ----------------------------------------------------------------
model <- "04_compoisson_bi_full"
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))

obj <- MakeADFun(data, parameters, DLL=model, random = "U")
obj$fn()     #1073.514  valor original (1098,039 agora, último)
system.time(fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr))
# 1h10 com 1000 TMB
# 11min com 200 TMB
# 7.4min com 200 R
z_to_r(.5)
# obj$gr()
4184/60
fit_NLB$par

print(obj$report())
?nlminb
system.time(fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr))
system.time(fit_NM <- optim(par = obj$par, fn = obj$fn, gr = obj$gr))
system.time(fit_BFGS <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS"))
system.time(fit_CG <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "CG"))

# As diferenças nos valores de rho e sigma são em virtude da reparametrização
# Coeficientes ótimos "True" e "True_10000" são encontrados pela parametrização de Varcov
true <- c(beta1 = 2.072, beta1 = 0.811, beta2 = 0.602, beta2 = -1.185,
          rho = 0.852, sigma = -0.764, sigma = -1.358, phi = -0.04, phi = -0.391)
true_10000 <- c(beta1 = 2.034, beta1 = 0.777, beta2 = 0.548, beta2 = -1.006, 
                rho = 0.789, sigma = 0.68, sigma = 0.462, phi = -0.035, phi = -0.108)

cbind(Parametrizacao_Wag = true, 
      Parametrizacao_TMB= round(fit_NLB$par,3))

## Maximized log-likelihood values (1062.906)
cbind("Nelder-Mead" = fit_NM$value, "BFGS" = fit_BFGS$value,
      "CG" = fit_CG$value, "NLB" = fit_NLB$objective)

## Estimates
cbind("True" =  True, "Nelder-Mead" = fit_NM$par, "BFGS" = fit_BFGS$par,
      "CG" = fit_CG$par, "NLB" = fit_NLB$par)

cbind("True" =  True, "NLB" = fit_NLB$par)

# Inverso da Transf Fisher Z
inv.transf.fisher <- function(z){ 
  r <- (exp(2*z)-1)/(exp(2*z)+1)
  return(r)
}

# Otimizado
pars_both <- c(beta1 = 2.072, beta1 = 0.811, beta2 = 0.602, beta2 = -1.185, 
               rho = 0.959, sigma = 0.682, sigma = 0.507, phi = -0.04, phi = -0.391)
nll_both <- 1062.906

# Perfilando
# Profile na mão
par(mfrow=c(3,2), mai = rep(.3,4))
par(mfrow=c(1,1), mai = rep(1,4))
# Beta1
p1 <- tmbprofile(obj,1, trace = F)
plot(p1, main = "beta1")
# Rho
dev.off()
p5 <- tmbprofile(obj,5, parm.range = c(-.5,10), trace =F)
plot(p5, main = "rho")
# Sigma1
p6 <- tmbprofile(obj,6, trace = F)
plot(p6, main = "Sigma1")
# Sigma2
p7 <- tmbprofile(obj,7, trace = F)
plot(p7, main = "Sigma2")
# Phi1
p8 <- tmbprofile(obj,8, trace = F)
plot(p8, main = "phi1")
# Phi2
p9 <- tmbprofile(obj,9, trace = F)
plot(p9, main = "phi2")

## Obtaining parameters in a logical way
ss <- sdreport(obj)
rep_fixed <- summary(ss, "fixed")
rep_report <- summary(ss, "report")
source("/home/guilherme/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels/unify_packages.R")

saveRDS(round(TMB_summary(rep_fixed, rep_report), 4), "results_COM_Poisson_from_NB.rds")

# Variance
library(glmmTMB)
m1 <- glmmTMB(Nmed ~ actdays, data = ahs, family = compois(link = "log"))
summary(m1)
str(m1)

summary(ahs)
######## Interpreting parameters under Multivariate Normal via UNSTRUCTURRED_CORR and VSCALE
### Interpreting Sigma
# Comparing values between VSCALE and Sigma não parametrizada
# na UNSTRUCTURED ele estima sigma (nao sigma2)
s2_1 <- -.764
(s_1 <- sqrt(exp(s2_1)))
s2_2 <- -1.358
(s_2 <- sqrt(exp(s2_2)))
# Obtaining the true value of sigma (log(sigma2))
exp(s2_1) #Close to .5
exp(s2_2) #Close to .3

### Interpreting rho
# Inv da Transf de Fisher Z para a cor
z <- .852
(exp(2*z)-1)/(exp(2*z)+1) # rho estimado Wagner: .6921129
# O TMB não solta o RHO, e nem a covariância, ele solta o Theta (https://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html):
z <- .958
(cov <- (exp(2*z)-1)/(exp(2*z)+1)) # Inverso da Transf Fisher Z
cov/(s_1*s_2)

# Getting rho from UNSTRUCTURRED TMB -------------------------------------------------------
# 

theta <- fit_NLB$par[5]
(sigma <- fit_NLB$par[6:7])   #s
# Função que retorna a matriz Sigma (VarCov ou VarCorr!)
getSigma <- function(theta, sigma, byrow = T, varcov = T) {
  ns <- length(sigma)
  L <- diag(1, nrow = ns, ncol = ns)
  if (byrow){ #Fill in theta by row/col
    L[upper.tri(L, diag=FALSE)] <- theta
    L <- t(L)
  } else{
    L[lower.tri(L)] <- theta
  }
  
  recheio <- L %*% t(L)    
  D <- diag(recheio)
  Dinvsqrt <- diag(1/sqrt(D))     
  if (varcov){ # VarCov
    W <- diag(sigma)     
    out <- W %*% Dinvsqrt %*% recheio %*% Dinvsqrt %*% t(W) #retorna a matriz de variância-covariancia
  } else{      # VarCorr
    out <- Dinvsqrt %*% recheio %*% Dinvsqrt
    diag(out) <- sigma^2
  }
  return(out) 
}
getSigma(theta, sigma, varcov = F) # Var e correlation
cov2cor(getSigma(theta, sigma))    # Checa a correlação
getSigma(theta, sigma)             # Valores das variâncias batem
theta1 <- getSigma(theta, sigma_e)[1,2]
inv.transf.fisher(theta1)
z <- .2395617 # Ao considerar o s2 na conta
(exp(2*z)-1)/(exp(2*z)+1)
z <- .6922681 # Sem considerar o W na conta
(exp(2*z)-1)/(exp(2*z)+1)
z <- 3.384536 # Sem considerar o W e D matrix na conta




### Rubbish ---------------------------------------

# # Gera o dado
# # set.seed(1234)
# Y1 <- compoisson::rcom(n, lambda = lambda1, nu = nu1)
# 
# plot(Y1~x1)
# summary(Y1)
# summary(Y2)
# plot(Y1~mu1)
# plot(Y1~mu2)
# 
# # set.seed(32)
# Y2 <- compoisson::rcom(n, lambda = lambda2, nu = nu2)
# range(Y1)
# range(Y2)
# #### Study
# # cbind(lambda1, lambda2)
# range(mu1)
# range(mu2)
# 
# range(lambda1)
# range(lambda2)
# plot(lambda1 ~ lambda2)


# Regression adjust via compoisson
# tab <- table(Y1)
# mat <- matrix(0, nrow = length(tab), ncol = 2)
# mat[, 1] <- as.numeric(rownames(tab))
# mat[,2] <- as.numeric(tab)
# fit1 <- com.fit(mat)
# summary(fit1)
# 
# # Regression adjust via mpcmp
# ds <- data.frame(Y1 = Y1, 
#                  X = X)
# mpcmp::glm.cmp(Y1 ~ X.x1, ds)
# head(ds)
# 
# data(insurance)
# com.fit(Lemaire);
# 
# head(mat,4)
# 
# compoisson::dcom(Y1, )
# 
# mpcmp::dcomp()
# 