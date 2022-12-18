## Negative Binomial distribution ----------------------------------------------
## Authors: Multivariate Statistical Modelling Group ---------------------------
## Date: June, 04 2020 ---------------------------------------------------------
rm(list=ls())
require(MASS)
require(mvtnorm)
require(TMB)
## Data simulation ------------------------------------------------------------
beta1 <- c(2, 0.8)
beta2 <- c(0.5, -1)
n <- 1000
# n <- 10000
set.seed(1234)
x1 <- rnorm(n)
X <- model.matrix(~x1)

## Random effects covariance matrix
s1 <- log(0.3) #.5 é s2
s2 <- log(0.3) #.3 é s2
# Transformação Z de Fisher
rho <- 0.5*log((1 + 0.5)/ (1-0.5)) # .5 é o verdadeira valor de rho
Sigma <- matrix(NA, 2, 2)
Sigma[1,1] <- exp(s1)
Sigma[2,2] <- exp(s2)
Sigma[1,2] <- rho*sqrt(exp(s1))*sqrt(exp(s2))
Sigma[2,1] <- rho*sqrt(exp(s1))*sqrt(exp(s2))
set.seed(1234)
U <- rmvnorm(n, mean = c(0,0), sigma = Sigma)
mu1 <- exp(X%*%beta1 + U[,1])
mu2 <- exp(X%*%beta2 + U[,2])

## Response variables
phi0 <- log(50)
phi1 <- log(50)
set.seed(1234)
Y1 <- rnbinom(n, size = exp(phi0), mu = mu1)
set.seed(1234)
Y2 <- rnbinom(n, size = exp(phi1), mu = mu2)

## Simulated data set ------------------------------------------------------------
data <- list(Y1 = Y1, Y2 = Y2, X = X)
## Parameters
True <- c(2,0.8,0.5,-1,rho,log(0.3),phi0,phi1)
parameters <- list(beta1 = c(2, 0.8), beta2 = c(0.5, -1),
                   U = matrix(0, ncol = 2, nrow = n),
                   rho = 0, sigma = c(s1), phi = c(phi0, phi1))

## Real data set ------------------------------------------------------------
## Dados Reais
library(mcglm)
# apply(ahs[, c("Nadm", "Ndoc", "actdays")], 2, var)
# apply(ahs[, c("Nadm", "Ndoc", "actdays")], 2, mean)
set.seed(2390)
ss <- sample(1:nrow(ahs), 1500, replace = F)
ahs <- ahs[ss, ]
data <- list(Y1 = ahs$Nadm,
             Y2 = ahs$Ndoc,
             X = model.matrix(~actdays, ahs))

var_to_phi <- function(mu, var){
  (mu^2)/(var-mu)
}
var(ahs$Nadm)
n <- nrow(ahs)
parameters <- list(beta1 = c(log(mean(ahs$Nadm)), 0), 
                   beta2 = c(log(mean(ahs$Ndoc)), 0),
                   U = matrix(0, ncol = 2, nrow = n),
                   rho = with(ahs, cor(Nadm, Ndoc)), 
                   sigma = .3, 
                   phi = c(var_to_phi(mean(ahs$Nadm), var(ahs$Nadm)),
                           var_to_phi(mean(ahs$Ndoc), var(ahs$Ndoc))))

## Compilation ------------------------------------------------------------
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
# model <- "06_negbin_bivariate_var_comum"
model <- "06_negbin_bivariate_wag_var_comum"
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))

# https://stackoverflow.com/questions/14719349/error-c-stack-usage-is-too-close-to-the-limit/14719448#14719448
fit_NLB$objective #5600
TMB::openmp(4)
# options(expressions=500000)
obj <- MakeADFun(data, parameters, DLL=model, random = "U")
obj$fn()     #1073.514  valor original (1098,039 agora, último)
obj$fn(True)
# obj$gr()

print(obj$report())
?nlminb
system.time(fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)) #44seg
system.time(fit_NM <- optim(par = obj$par, fn = obj$fn, gr = obj$gr)) #2.8min
system.time(fit_BFGS <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS")) #35seg
system.time(fit_CG <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "CG")) # 2.083seg

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
fit_NLB <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
                  control = list(iter.max = 1000, eval.max = 1000,
                                 step.min = 1, step.max = 1))
par(mfrow=c(3,2), mai = rep(.3,4))
# Beta1
p1 <- tmbprofile(obj,1, trace = F)
plot(p1, main = "beta1")
# Rho
dev.off()
p5 <- tmbprofile(obj,5, parm.range = c(-.5,7), trace =F)
plot(p5, main = "rho")
# Sigma1
p6 <- tmbprofile(obj,6, trace = F)
plot(p6, main = "Sigma1")
# Phi1
p8 <- tmbprofile(obj,7, trace = F)
plot(p8, main = "phi2")
par(mfrow=c(1,2), mai = rep(.9,4))
plot(p8, main = "phi1", parm.range = c(0,10))
par(mfrow=c(1,3))
with(p8[p8$phi<5,], plot(value ~ phi, type = "l", main = "entre 0 e 5"))
with(p8[p8$phi<15,], plot(value ~ phi, type = "l", main = "entre 0 e 15"))
with(p9[p9$phi<15,], plot(value ~ phi, type = "l", main = "entre 0 e 15"))
with(p8[p8$phi<20,], plot(value ~ phi, type = "l", main = "entre 0 e 20"))
# Phi2
p9 <- tmbprofile(obj,8, trace = F)
plot(p9, main = "phi2")
plot(p9, main = "phi2", parm.range = c(10,15))

###### Obtaining correlation parameters
ss <- sdreport(obj)
tab <- summary(ss, "fixed", p.value = T)
tab_optim <- tab
tab_fit_NLB <- tab
row.names(tab)[c(1:4,7,8)] <- c(paste0("Nadm_", "beta0"), paste0("Nadm_", "beta1"),
                                paste0("Ndoc_", "beta0"), paste0("Ndoc_", "beta1"),
                                paste0("Nadm_", "phi"), paste0("Ndoc_", "phi"))
tab2 <- summary(ss, "report", p.value = T)
rownames(tab2) <- c("sd_comum","rho")
tab2 <- tab2[c(2,1), ]
tab2
tab[rownames(tab)==c("rho","sigma"), ] <- tab2
tab <- tab[c(1,3,2,4,7,8,6,6,5),]
dput(tab)
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
getSigma(theta, sigma) #Valores das variâncias batem
getSigma(theta, sigma, varcov = F) #Valores das variâncias batem
theta1 <- getSigma(theta, sigma_e)[1,2]
inv.transf.fisher(theta1)
z <- .2395617 # Ao considerar o s2 na conta
(exp(2*z)-1)/(exp(2*z)+1)
z <- .6922681 # Sem considerar o W na conta
(exp(2*z)-1)/(exp(2*z)+1)
z <- 3.384536 # Sem considerar o W e D matrix na conta

