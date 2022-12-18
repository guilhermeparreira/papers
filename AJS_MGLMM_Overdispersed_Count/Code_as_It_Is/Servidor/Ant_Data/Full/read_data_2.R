rm()
path <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor"
localp <- "/Ant_Data/Full"
locals <- "/Ant_Data/SmallerModels"
setwd(paste0(path, localp))
##########################
# Model 1 - Poisson
##########################
load("Poisson.RData")
# ls()
# estimates_poisson_am
# estimates_poisson_am_2
# estimates_poisson
# estimates_poisson_s
# estimates_poisson_sd
# estimates_poisson_s2 # BFGS MELHOR
estimates_poisson_sd$summary$Model <- "Full Poisson "
estimates_poisson_sd$summary$np <- nrow(estimates_poisson_sd$ests)
estimates_poisson_sd$summary$SE <- TRUE
m01 <- estimates_poisson_sd$summary
rm(list = ls()[!ls()%in%c("m01", "path", "localp", "locals")])
##########################
# Model 2 - Negbin 
##########################
load("Negbin_full.RData")
# estimates_negbin$summary
# estimates_negbin_sd$summary
# estimates_negbin_am$summary
estimates_negbin_2$summary$Model <- "Full NB"
estimates_negbin_2$summary$np <- prod(dim(estimates_negbin_2$beta))+length(estimates_negbin_2$sigma)+length(estimates_negbin_2$rho)+length(estimates_negbin_2$disp)
estimates_negbin_2$summary$SE <- FALSE
m02 <- estimates_negbin_2$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "path", "localp", "locals")])
##########################
# Model 3 - Negbin - Dispersion Fixed 
##########################
setwd(paste0(path, locals))
load("Negbin_phi_fixed.RData")
estimates_negbin$summary$Model <- "Fixed Disp NB"
estimates_negbin$summary$np <- prod(dim(estimates_negbin$beta))+length(estimates_negbin$sigma)+length(estimates_negbin$rho)
estimates_negbin$summary$SE <- FALSE
m03 <- estimates_negbin$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "path", "localp", "locals")])
##########################
# Model 4 - Negbin - Comum Variance
##########################
load("Negbin_var_comum.RData")
estimates_negbin$summary$Model <- "Comum Var NB"
estimates_negbin$summary$np <- prod(dim(estimates_negbin$beta))+length(estimates_negbin$sigma)+length(estimates_negbin$disp)+length(estimates_negbin$rho)
estimates_negbin$summary$SE <- FALSE
m04 <- estimates_negbin$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04" ,"path", "localp", "locals")])
##########################
# Model 5 - Negbin - Fixed Variance
##########################
load("Negbin_var_fixed.RData")
estimates_negbin$summary$Model <- "Fixed Var NB"
estimates_negbin$summary$np <- prod(dim(estimates_negbin$beta))+length(estimates_negbin$sigma)+length(estimates_negbin$disp)+length(estimates_negbin$rho)
estimates_negbin$summary$SE <- FALSE
m05 <- estimates_negbin$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05" ,"path", "localp", "locals")])
##########################
# Model 6 - Compoisson FULL
##########################
setwd(paste0(path, localp))
load("Compoisson_full.RData")
ls()
# fit_compoisson
# fit_compoisson2
# fit_compoisson_s
# estimates
# estimates_comp
# estimates_comp_s
ss_sd$summary$Model <- "Full COM-Poisson"
ss_sd$summary$np <- nrow(ss_sd$ests)
ss_sd$summary$SE <- FALSE
m06 <- ss_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06" ,"path", "localp", "locals")])
##########################
# Model 7 - Compoisson Dispersion Fixed 
##########################
setwd(paste0(path, locals))
load("compoisson_nu_fixed.RData")
ls()
# estimates_compoisson_s
estimates_compoisson_sd$summary$Model <- "Fixed Disp COM-Poisson"
estimates_compoisson_sd$summary$np <- nrow(estimates_compoisson_sd$ests)
estimates_compoisson_sd$summary$SE <- TRUE
m7 <- estimates_compoisson_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m7" ,"path", "localp", "locals")])
##########################
# Model 8 - Compoisson Comum Var
##########################
load("compoisson_var_comum.RData")
ls()
estimates_compoisson_sd$summary$Model <- "Comum Var COM-Poisson"
estimates_compoisson_sd$summary$np <- nrow(estimates_compoisson_sd$ests)
estimates_compoisson_sd$summary$SE <- TRUE
m8 <- estimates_compoisson_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m7", "m8","path", "localp", "locals")])
##########################
# Model 9 - Compoisson Fixed Var
##########################
load("compoisson_var_fixed.RData")
ls()
estimates_compoisson_sd$summary$Model <- "Fixed Var COM-Poisson"
estimates_compoisson_sd$summary$np <- nrow(estimates_compoisson_sd$ests)
estimates_compoisson_sd$summary$SE <- TRUE
m9 <- estimates_compoisson_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m7", "m8", "m9" ,"path", "localp", "locals")])
##########################
# Model 10 - Poisson - Rho fixed
##########################
setwd(paste0(path, locals))
load("Poisson_rho_fixed.RData")
estimates_poisson_s1$summary$Model <- "Rho fixed Poisson"
estimates_poisson_s1$summary$np <- prod(dim(estimates_poisson_s1$beta))+length(estimates_poisson_s1$sigma)
estimates_poisson_s1$summary$SE <- FALSE
m10 <- estimates_poisson_s1$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m7", "m8", "m9", "m10" ,"path", "localp", "locals")])
##########################
# Model 11 - NB - Rho fixed
##########################
setwd(paste0(path, locals))
load("Negbin_rho_fixed.RData")
estimates_negbin_1$summary$Model <- "Rho fixed NB"
estimates_negbin_1$summary$np <- prod(dim(estimates_negbin_1$beta))+length(estimates_negbin_1$sigma)+length(estimates_negbin_1$disp)
estimates_negbin_1$summary$SE <- FALSE
m11 <- estimates_negbin_1$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m7", "m8", "m9", "m10", "m11" ,"path", "localp", "locals")])
##########################
# Model 12 - COM-Poisson - Rho fixed
##########################
setwd(paste0(path, locals))
load("Compoisson_rho_fixed.RData")
ls()
ss$summary$Model <- "Rho fixed COM-Poisson"
ss$summary$np <- prod(dim(ss$beta))+length(ss$sigma)+length(ss$disp)
ss$summary$SE <- FALSE
m12 <- ss$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m7", "m8", "m9", "m10", "m11", "m12","path", "localp", "locals")])
#####################################################3
####### Joining all results
#####################################################3
# sample size 1281
mm <- rbind(m01,m10,
            m02,m03,m04,m05,m11,
            m06, m7, m8, m9, m12)
mm$Loglik <- mm$Loglik*-1
mm$AIC <- -2*mm$Loglik + 2*mm$np
mm$BIC <- -2*mm$Loglik + log(1281)*mm$np

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
mm$SE <- ifelse(mm$SE=="TRUE", "\\checkmark", "\\xmark")
library(knitr)
mm
kable(mm, 
      caption = "Goodness-of-fit measures for ANT data from different distributions and specifications",
      label = "antfit",
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











