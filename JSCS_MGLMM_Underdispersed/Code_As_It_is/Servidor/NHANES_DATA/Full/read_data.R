rm()
path <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor"
localp <- "/NHANES_DATA/Full"
locals <- "/NHANES_DATA/SmallerModels"
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
load("negbin.RData")
# estimates_negbin$summary
# estimates_negbin_sd$summary
# estimates_negbin_am$summary
estimates_negbin_sd$summary$Model <- "Full NB"
estimates_negbin_sd$summary$np <- nrow(estimates_negbin_sd$ests)
estimates_negbin_sd$summary$SE <- TRUE
m02 <- estimates_negbin_sd$summary
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
m4 <- estimates_negbin$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4" ,"path", "localp", "locals")])
##########################
# Model 5 - Negbin - Fixed Variance
##########################
load("Negbin_var_fixed.RData")
estimates_negbin_sd$summary$Model <- "Fixed Var NB"
estimates_negbin_sd$summary$np <- nrow(estimates_negbin_sd$ests)
estimates_negbin_sd$summary$SE <- TRUE
m5 <- estimates_negbin_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5" ,"path", "localp", "locals")])
##########################
# Model 6 - Compoisson FULL
##########################
setwd(paste0(path, localp))
load("comp.RData")
estimates_comp$summary$Model <- "Full COM-Poisson"
estimates_comp$summary$np <- nrow(estimates_comp$ests)
estimates_comp$summary$SE <- TRUE
m6 <- estimates_comp$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6" ,"path", "localp", "locals")])
##########################
# Model 7 - Compoisson Dispersion Fixed 
##########################
setwd(paste0(path, locals))
load("compoisson_nu_fixed_1.5.RData")
# estimates_compoisson_s
estimates_compoisson$summary$Model <- "Fixed Disp COM-Poisson"
estimates_compoisson$summary$np <- prod(dim(estimates_compoisson$beta))+length(estimates_compoisson$sigma)+length(estimates_compoisson$disp)+length(estimates_compoisson$rho)
estimates_compoisson$summary$SE <- FALSE
m7 <- estimates_compoisson$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7" ,"path", "localp", "locals")])
##########################
# Model 8 - Compoisson Comum Var
##########################
load("compoisson_var_comum.RData")
estimates_compoisson_sd$summary$Model <- "Comum Var COM-Poisson"
estimates_compoisson_sd$summary$np <- nrow(estimates_compoisson_sd$ests)
estimates_compoisson_sd$summary$SE <- TRUE
m8 <- estimates_compoisson_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7", "m8","path", "localp", "locals")])
##########################
# Model 9 - Compoisson Fixed Var
##########################
# load("compoisson_var_fixed.RData")
m9 <- data.frame(Loglik = 33039.9, Convergence = 1, Threads = 12,
                 Time = NA, Optimizer = "nlminb", Model = "Fixed Var COM-Poisson",
                 np = 21, SE = TRUE)
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7", "m8", "m9" ,"path", "localp", "locals")])
##########################
# Model 10 - Poisson rho 0
##########################
setwd(paste0(path, locals))
load("Poisson_rho_fixed.RData")
# ls()
estimates_poisson$summary$Model <- "Rho zero Poisson"
estimates_poisson$summary$np <- prod(dim(estimates_poisson$beta))+length(estimates_poisson$sigma)
estimates_poisson$summary$SE <- FALSE
m10 <- estimates_poisson$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7", "m8", "m9", "m10","path", "localp", "locals")])
##########################
# Model 11 - NB rho 0
##########################
setwd(paste0(path, locals))
load("negbin_rho_fixed.RData")
estimates_negbin0$summary$Model <- "Rho zero NB"
estimates_negbin0$summary$np <- prod(dim(estimates_negbin0$beta))+length(estimates_negbin0$sigma)+length(estimates_negbin0$disp)
estimates_negbin0$summary$SE <- FALSE
m11 <- estimates_negbin0$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11","path", "localp", "locals")])
##########################
# Model 12 - CMP rho 0
##########################
setwd(paste0(path, locals))
load("comp_rho_fixed.RData")
estimates_comp$summary$Model <- "Rho zero COM-Poisson"
estimates_comp$summary$np <- nrow(estimates_comp$ests)-3
estimates_comp$summary$SE <- TRUE
m12 <- estimates_comp$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12","path", "localp", "locals")])
##########################
# Model 13 - DoublePoisson FULL
##########################
setwd(paste0(path, localp))
load("double_v2.RData")
estimates_double_5$summary$Model <- "Full double Poisson"
estimates_double_5$summary$np <- nrow(estimates_double_5_sd$ests)
estimates_double_5$summary$SE <- TRUE
m13 <- estimates_double_5$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "path", "localp", "locals")])
##########################
# Model 14 - DoublePoisson Dispersion Fixed 
##########################
setwd(paste0(path, locals))
load("doublepoisson_nu_fixed.RData")
estimates_double_nu_fixed_v3_sv$summary$Model <- "Fixed Disp double Poisson"
estimates_double_nu_fixed_v3_sv$summary$np <- prod(dim(estimates_double_nu_fixed_v3_sv$beta))+length(estimates_double_nu_fixed_v3_sv$sigma)+length(estimates_double_nu_fixed_v3_sv$disp)+length(estimates_double_nu_fixed_v3_sv$rho)
estimates_double_nu_fixed_v3_sv$summary$SE <- TRUE
m14 <- estimates_double_nu_fixed_v3_sv$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "m14", "path", "localp", "locals")])
##########################
# Model 15 - DoublePoisson Comum Var
##########################
load("doublepoisson_var_comum.RData")
estimates_double_3_sd$summary$Model <- "Comum Var double Poisson"
estimates_double_3_sd$summary$np <- nrow(estimates_double_3_sd$ests)
estimates_double_3_sd$summary$SE <- TRUE
m15 <- estimates_double_3_sd$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "m14", "m15", "path", "localp", "locals")])
##########################
# Model 16 - DoublePoisson Fixed Var
##########################
load("doublepoisson_var_fixed.RData")
estimates_doublepoisson_sd_2$summary$Model <- "Fixed Var double Poisson"
estimates_doublepoisson_sd_2$summary$np <- nrow(estimates_doublepoisson_sd_2$ests)
estimates_doublepoisson_sd_2$summary$SE <- TRUE
m16 <- estimates_doublepoisson_sd_2$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "m14", "m15", "m16", "path", "localp", "locals")])
##########################
# Model 17 - Double Poisson rho 0
##########################
load("doublepoisson_rho_fixed.RData")
estimates_double_sd_2$summary$Model <- "Rho zero double Poisson"
estimates_double_sd_2$summary$np <- nrow(estimates_double_sd_2$ests)-3
estimates_double_sd_2$summary$SE <- TRUE
m17 <- estimates_double_sd_2$summary
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "m14", "m15", "m16", "m17", "path", "localp", "locals")])

#####################################################3
####### Joining all results
#####################################################3
# sample size 1281
mm <- rbind(m01,m10,
            m02,m03,m4,m5,m11,
            m6, m7, m8, m9, m12,
            m13, m14, m15, m16, m17)
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
mm$Model <- as.character(pri.maiuscula(gsub("Full ","", gsub("Disp ", "dispersion ", gsub("Comum", "Common", gsub("Var ", "variance ", mm$Model))))))
mm$SE <- ifelse(mm$SE=="TRUE", "\\checkmark", "\\xmark")
library(knitr)
kable(mm, 
      caption = "Goodness-of-fit measures for NHANES data from different distributions and specifications",
      label = "nhanesfit",
      align = c("l", rep("r", ncol(mm)-1)),
      booktabs = T,
      digits = 1,
      row.names = F,
      escape = F,
      format = "latex")




# minha conclusão: os algoritmos de maximização externa ainda não são espertos o suficiente para estimar esse modelo.
# veja o caso da Binomial Negative: com dispersão fixa, maior loglik que o modelo completo. claramente não teve uma estimação satisfatória.
# ainda, um resultado que algumas vezes era convergencia 1 para o PORT, era 0 para o BFGS
# valor da loglik para compoisson completa muito menor e distante ao fixar a dispersão.