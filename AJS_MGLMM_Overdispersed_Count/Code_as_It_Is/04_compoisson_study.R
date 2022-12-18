# Eduardo way to generate random numbers -------------------------
x = 0:10
px = mpcmp::dcomp(x, mu = 7, nu = 1)
sample(x, size = length(x), prob = px, replace = T)
# Só que é ineficiente pq vc tem que fazer isso pra casa mu_i
# E ainda, vc está trancando a geração de números aleatórios entre 0 e 10

# Different way to generate a random number
set.seed(123)
nu <- .1
mode <- 10
domain <- 0:100
prob <- dpois(domain, lambda=mode)^nu; prob <- prob / sum(prob)
sum(prob * domain) ## mean
x <- sample(domain, size=1e4, replace=TRUE, prob = prob)

## Github compoisson ---------------------------------------------
library(TMB)
compile("compois.cpp")
dyn.load(dynlib("compois"))
## TMB data
data <- list( x = x )
parameters <- list( logmu = 0, lognu = 0 )

## Parameterization through the mode
data$parameterization <- "mode"
obj <- MakeADFun(data, parameters, DLL="compois")
system.time(fit.mode <- nlminb(obj$par, obj$fn, obj$gr, obj$he))
rep.mode <- sdreport(obj)

## Parameterization through the mean
data$parameterization <- "mean"
obj <- MakeADFun(data, parameters, DLL="compois")
system.time(fit.mean <- nlminb(obj$par, obj$fn, obj$gr, obj$he))
rep.mean <- sdreport(obj)

summary(rep.mode, "report")
summary(rep.mean, "report")


## compoisson package ---------------------------------------------------------------------
lambda = 3
nu = .5

dcom(1, 3, .5)
com.mean(lambda = 3, nu = .5)
com.log.density(1, 3, .5, log.z = NULL)

# Steps to calculate the density
# First: Compute Z
com.compute.log.z(lambda, nu, log.error = 0.0001)
# Extra: Calculate the mean
com.expectation(function(x) x, lambda, nu, log.error = 0.001)


The COM-Poisson no pacote compoisson esta parametrizada em Lambda e nu.
The COM-Poisson no pacote TMB esta parametrizada em Mean e nu.

O pacote compoisson tem a função com.mean(Lambda, nu) que retorna a média.
Aí, o que eu fiz: 

compoisson::dcom(x = 1, lambda = 3, nu = .5) # retorna 0.008664
mpcmp::dcomp(1, mu = com.mean(lambda = 3, nu = .5), nu = .5)
TMB::dcompois2(x = 1, mean = com.mean(lambda = 3, nu = .5), nu = .5) #retorna .0086725 

mpcmp::rcomp(1, mu = 500, nu = 1)
exp(mu2lambda(50, 1)$lambda)
compoisson::rcom(1, 5000, 1)
mpcmp::rcomp()
rpois(1,9000)

Agora, como a partir da média e nu, eu chego no lambda? (pq aí, a média é o mu1, nu eu estabeleço; com ambos, eu conseguiria calcular o lambda, e gerar os valores)

lambda2mu <- function(lambda, nu) {
  phi <- log(nu)
  mu <- lambda^(1/nu) - (nu - 1)/(2 * nu)
  list("mu" = mu, "phi" = phi)
}

# Usar a função mu2lambda
http://www.leg.ufpr.br/doku.php/pessoais:eduardojr:publications-rbras2017
mu2lambda <- function(mu, phi) {
  nu <- exp(phi)
  lambda <- (mu + (nu-1)/(2*nu))^nu
  list("lambda" = lambda, "nu" = nu)
}
mu2lambda(500, 1)
range(rcom(1500, mu2lambda(mu1, .5)$lambda, .5))

## mpcmp package ------------------------------------------------------------------------
mpcmp::rcomp(n, mu = mu1, nu = .5)
mpcmp::rcomp(10000, mu = 9, nu = .5)
lmu1 <- length(mu1)
y <- numeric(lmu1)
for(i in 1:lmu1){
  y[i] <- mpcmp::rcomp(1, mu1[i], nu = .5)
}

## TMB simulation ------------------------------------------------------------
rm(list=ls())
require(mvtnorm)
require(TMB)
require(mpcmp)
beta1 <- c(log(10), .7)
beta2 <- c(log(10), -1)
n <- 50
# n <- 10000
# set.seed(1234)
x1 <- rnorm(n)
X <- model.matrix(~x1)
## Random effects covariance matrix
s1 <- log(.3) #.3 é s2
s2 <- log(.15) #.15 é s2
# Transformação Z de Fisher
rho <- 0.5*log((1 + 0.5)/ (1-0.5)) # .5 é o verdadeira valor de rho
Sigma <- matrix(NA, 2, 2)
Sigma[1,1] <- exp(s1)
Sigma[2,2] <- exp(s2)
Sigma[1,2] <- rho*sqrt(exp(s1))*sqrt(exp(s2))
Sigma[2,1] <- rho*sqrt(exp(s1))*sqrt(exp(s2))
# set.seed(61)
U <- rmvnorm(n, mean = c(0,0), sigma = Sigma)
mu1 <- exp(X%*%beta1 + U[,1])
plot(mu1~x1)
mu2 <- exp(X%*%beta2 + U[,2])
plot(mu2~x1)
## Response variables
# Define o nu
nu1 <- 1
nu2 <- 1
# nu1 <- .95
# nu2 <- .95
l_nu1 <- log(nu1) # .5 é o nu
l_nu2 <- log(nu2) # .5 é o nu
Y1 <- mpcmp::rcomp(n, mu = mu1, nu = nu1)
Y2 <- mpcmp::rcomp(n, mu = mu2, nu = nu2)
plot(log(Y2)~mu2)
plot(log(Y1)~mu2)

data <- list(Y = Y1, mu = X)
## Parameters

## Compilation
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
model <- "04_compoisson"
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))

Rinterface(paste0(model, ".cpp"))
obj <- MakeADFun(data, 
                 parameters = list(), 
                 # type = "Fun",
                 DLL=model
                 # ,checkParameterOrder=FALSE
                 )




## general study -------------------------------------------------------------------------
# COM-Poisson study
# lambda <- c(.5, 1, 5, 10, 30, 50, 100, 150, 200)
# llamb <- length(lambda)
# nu <- seq(0, 1, .1)
# ln <- length(nu)
# # nu < 1 is underdispersion
# mat <- matrix(0, nrow = ln, ncol = llamb,
#        dimnames = list(nu, lambda))
# for (i in 1:ln){
#   for (j in 1:llamb) {
#     # j <- 1
#     # i <- 1
#     value <- 0
#     # k <- 1
#     for (k in 1:1000){
#       # value <- value + k*log(lambda[j])-nu[i]*lfactorial(k)
#       value <- value + (lambda[j]^k)/(factorial(k)^(nu[i]))
#     }
#     mat[i, j] <- value
#   }
# }

# Rubbish ------------------------------

# mu2lambda <- function(mu, phi) {
#   nu <- exp(phi)
#   lambda <- (mu + (nu-1)/(2*nu))^nu
#   list("lambda" = lambda, "nu" = nu)
# }
# lambda2mu <- function(lambda, nu) {
#   phi <- log(nu)
#   mu <- lambda^(1/nu) - (nu - 1)/(2 * nu)
#   list("mu" = mu, "phi" = phi)
# }

# mu2lambda()
# lambda = 3
# nu = .5
# com.mean(lambda, nu)
# mu2lambda(9.5212, .5)
# lambda2mu(3, .5)
# 
# # Work
# dcom(1, lambda, nu)
# dcomp(1, lambda = lambda, nu = nu)
# 
# # 
# lambda2mu(3, .5)
# dcomp(1, mu = 9.5, nu = .5)


# é só aplicar o mu2lambda mesmo. O nu fica como está. 

