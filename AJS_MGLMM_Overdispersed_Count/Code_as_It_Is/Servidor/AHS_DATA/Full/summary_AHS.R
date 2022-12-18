rm()
path <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor"
localp <- "/AHS_DATA/Full"
locals <- "/AHS_DATA/SmallerModels"
setwd(paste0(path, localp))
##########################
# Model 1 - Poisson
##########################
load("Poisson_final.RData")
# ls()
# estimates_poisson_am
# estimates_poisson_am_2
# estimates_poisson
# estimates_poisson_s
# estimates_poisson_sd
# estimates_poisson_s2 # BFGS MELHOR
estimates_poisson_sd$summary$Model <- "Poisson "
estimates_poisson_sd$summary$np <- nrow(estimates_poisson_sd$ests)
estimates_poisson_sd$summary$SE <- TRUE
m01 <- estimates_poisson_sd$summary
rm(list = ls()[!ls()%in%c("m01", "path", "localp", "locals")])
##########################
# Model 2 - Negbin 
##########################
load("negbin.RData")
estimates_negbin_sd$summary$Model <- "NB"
estimates_negbin_sd$summary$np <- nrow(estimates_negbin_sd$ests)
estimates_negbin_sd$summary$SE <- TRUE
m02 <- estimates_negbin_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "path", "localp", "locals")])
##########################
# Model 3 - Negbin - Dispersion Fixed
##########################
setwd(paste0(path, locals))
load("Negbin_phi_fixed.RData")
# estimates_negbin$summary
# estimates_negbin_sd$summary
estimates_negbin_sd$summary$Model <- "Fixed dispersion NB"
estimates_negbin_sd$summary$np <- nrow(estimates_negbin_sd$ests)
estimates_negbin_sd$summary$SE <- TRUE
m03 <- estimates_negbin_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "path", "localp", "locals")])
##########################
# Model 4 - Negbin - Comum Variance
##########################
load("Negbin_var_comum.RData")
estimates_negbin_sd$summary$Model <- "Commom variance Var NB"
estimates_negbin_sd$summary$np <- nrow(estimates_negbin_sd$ests)
estimates_negbin_sd$summary$SE <- TRUE
m04 <- estimates_negbin_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04" ,"path", "localp", "locals")])
##########################
# Model 5 - Negbin - Fixed Variance
##########################
load("Negbin_var_fixed.RData")
estimates_negbin_sd$summary$Model <- "Fixed variance NB"
estimates_negbin_sd$summary$np <- nrow(estimates_negbin_sd$ests)
estimates_negbin_sd$summary$SE <- TRUE
m05 <- estimates_negbin_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05" ,"path", "localp", "locals")])
##########################
# Model 6 - Compoisson FULL (Próximo)
##########################
# Port primeira vez, e jogando para o port denovo, deu na mesma
setwd(paste0(path, localp))
load("compoisson_1.RData")
ls()
# fit_compoisson
# fit_compoisson2
# fit_compoisson_s
# estimates
# estimates_comp
# estimates_comp_s
estimates_compoisson_sd$summary$Model <- "COM-Poisson"
estimates_compoisson_sd$summary$np <- nrow(estimates_compoisson_sd$ests)
estimates_compoisson_sd$summary$SE <- TRUE
m6 <- estimates_compoisson_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m6" ,"path", "localp", "locals")])
##########################
# Model 7 - Compoisson Dispersion Fixed 
##########################
setwd(paste0(path, locals))
load("compoisson_nu_fixed.RData")
ls()
# estimates_compoisson_s
estimates_compoisson_sd$summary$Model <- "Fixed dispersion COM-Poisson"
estimates_compoisson_sd$summary$np <- nrow(estimates_compoisson_sd$ests)
estimates_compoisson_sd$summary$SE <- TRUE
m7 <- estimates_compoisson_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m6", "m7" ,"path", "localp", "locals")])
##########################
# Model 8 - Compoisson Comum Var
##########################
load("compoisson_var_comum_1.RData")
ls()
estimates_compoisson$summary$Model <- "Commom variance COM-Poisson"
estimates_compoisson$summary$np <- prod(dim(estimates_compoisson$beta)) + length(estimates_compoisson$sigma) + length(estimates_compoisson$rho) + length(estimates_compoisson$disp)
estimates_compoisson$summary$SE <- FALSE
m8 <- estimates_compoisson$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m6", "m7", "m8","path", "localp", "locals")])
##########################
# Model 9 - Compoisson Fixed Var
##########################
load("compoisson_var_fixed_1.RData")
ls()
estimates_compoisson$summary$Model <- "Fixed variance COM-Poisson"
estimates_compoisson$summary$np <- prod(dim(estimates_compoisson$beta)) + length(estimates_compoisson$sigma) + length(estimates_compoisson$rho) + length(estimates_compoisson$disp)
estimates_compoisson$summary$SE <- FALSE
m9 <- estimates_compoisson$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m6", "m7", "m8", "m9" ,"path", "localp", "locals")])
##########################
# Model 10 - Rho Fixed Poisson
##########################
setwd(paste0(path, locals))
load("Poisson_rho_fixed.RData")
estimates_poisson_sd$summary$Model <- "Rho fixed Poisson"
estimates_poisson_sd$summary$np <- nrow(estimates_poisson_sd$ests)
estimates_poisson_sd$summary$SE <- TRUE
m10 <- estimates_poisson_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m6", "m7", "m8", "m9", "m10" ,"path", "localp", "locals")])
##########################
# Model 11 - Rho Fixed NB
##########################
setwd(paste0(path, locals))
load("negbin_rho_fixed.RData")
estimates_negbin_sd1$summary$Model <- "Rho fixed NB"
estimates_negbin_sd1$summary$np <- nrow(estimates_negbin_sd1$ests)
estimates_negbin_sd1$summary$SE <- TRUE
m11 <- estimates_negbin_sd1$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m6", "m7", "m8", "m9", "m10", "m11" ,"path", "localp", "locals")])
##########################
# Model 12 - Rho Fixed COM-POISSON
##########################
setwd(paste0(path, locals))
load("compoisson_rho_fixed.RData")
estimates_compoisson_s$summary$Model <- "Rho fixed COM-Poisson"
estimates_compoisson_s$summary$np <- prod(dim(estimates_compoisson_s$beta)) + length(estimates_compoisson_s$sigma) + length(estimates_compoisson_s$disp)
estimates_compoisson_s$summary$SE <- FALSE
m12 <- estimates_compoisson_s$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m6", "m7", "m8", "m9", "m10", "m11", "m12" ,"path", "localp", "locals")])
#####################################################3
####### Joining all results
#####################################################3
# sample size 1281
mm <- rbind(m01,m10,
            m02,m03,m04,m05,m11,
            m6, m7, m8, m9,m12)
# mm <- rbind(m01,m02,m03,m04,m05,m6, m7)
mm$Loglik <- mm$Loglik*-1
mm$AIC <- round(-2*mm$Loglik + 2*mm$np)
mm$BIC <- round(-2*mm$Loglik + log(5190)*mm$np)
mm$Loglik <- round(mm$Loglik)
mm <- mm[, c(6,7,9,10,1,8,5)]
row.names(mm) <- NULL
mm$Optimizer <- with(mm, ifelse(Optimizer=="bfgs", "BFGS", 
                                ifelse(Optimizer%in%c("nlminb","nlimnb"), "PORT", Optimizer)))
mm <- mm[, -c(7)]

pri.maiuscula <- function(x) {
   s <- strsplit(x, " ")[[1]]
   paste(toupper(substring(s, 1,1)), substring(s, 2),
         sep="", collapse=" ")
}
pri.maiuscula <- Vectorize(pri.maiuscula, "x")
mm$Model <- as.character(pri.maiuscula(gsub("Full ","", gsub("Disp ", "dispersion ", gsub("Var ", "variance ", mm$Model)))))
mm$SE <- ifelse(mm$SE, "\\checkmark", "\\xmark")
library(knitr)
kable(mm, 
      caption = "GOODNESS-OF-FIT MEASURES FOR AHS DATA FROM DIFFERENT DISTRIBUTIONS AND SPECIFICATIONS",
      label = "ahsfit",
      align = c("l", rep("c", ncol(mm)-1)),
      booktabs = T,
      digits = 1,
      row.names = F,
      escape = F,
      format = "latex")

# minha conclusão: os algoritmos de maximização externa ainda não são espertos o suficiente para estimar esse modelo.
# veja o caso da Binomial Negative: com dispersão fixa, maior loglik que o modelo completo. claramente não teve uma estimação satisfatória.
# ainda, um resultado que algumas vezes era convergencia 1 para o PORT, era 0 para o BFGS
# valor da loglik para compoisson completa muito menor e distante ao fixar a dispersão.











