rm()
path <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor"
localp <- "/Ant_Data/Full"
locals <- "/Ant_Data/SmallerModels"
setwd(paste0(path, localp))
library(kableExtra)
library(reshape2)
library(latex2exp)
library(ggplot2)
library(tidyr)
library(dplyr)
library(knitr)
options(knitr.kable.NA = '')
#################################################################################
######################### Preciso fazer exp(disp) quando não tiver SE - tanto no NB como no COMPOISSON
#################################################################################
##########################
# Model 1 - Poisson
##########################
load("Poisson.RData")
# ls()
m01 <- as.data.frame(estimates_poisson_sd$ests)
m01$Model <- "Poisson"
m01$MainModel <- "Poisson"
rm(list = ls()[!ls()%in%c("m01", "path", "localp", "locals")])
##########################
# Model 2 - Negbin
##########################
load("Negbin_full.RData")
b1p <- melt(estimates_negbin_2$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_negbin_2$sigma)
r1p <- melt(estimates_negbin_2$rho)
d1p <- exp(melt(estimates_negbin_2$disp))
m02 <- data.frame(Estimate = unname(rbind(b1p, s1p, r1p, d1p)),
                  `Std. Error` = NA,
                  Z = NA,
                  `p value` = NA,
                  check.names = F,
                  row.names = c(row.names(m01), paste0('exp(phi)',1:41)))
m02$Model <- "NB"
m02$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "path", "localp", "locals")])
##########################
# Model 3 - Negbin - Disp Fixed
##########################
setwd(paste0(path, locals))
load("Negbin_phi_fixed.RData")
b1p <- melt(estimates_negbin$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_negbin$sigma)
r1p <- melt(estimates_negbin$rho)
m03 <- data.frame(Estimate = unname(rbind(b1p, s1p, r1p)),
           `Std. Error` = NA,
           Z = NA,
           `p value` = NA,
           check.names = F,
           row.names = row.names(m01))
m03$Model <- 'Fixed Disp NB'
m03$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "path", "localp", "locals")])
##########################
# Model 4 - Negbin - Comum Variance
##########################
load("Negbin_var_comum.RData")
b1p <- melt(estimates_negbin$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_negbin$sigma)
r1p <- melt(estimates_negbin$rho)
d1p <- exp(melt(estimates_negbin$disp)) # Transformação aqui.
nbetastotal = (length(colnames(data$X))*length(colnames(data$Y)))
ncorre = (nr*nr-nr)/2
m04 <- data.frame(Estimate = unname(rbind(b1p, s1p, r1p, d1p)),
                 `Std. Error` = NA,
                 Z = NA,
                 `p value` = NA,
                 check.names = F,
                 row.names = c(row.names(m01)[1:nbetastotal], 
                               paste0('sig',1), 
                               row.names(m01)[(nbetastotal+nr+1):(nbetastotal+nr+ncorre)], 
                               paste0('exp(phi)',1:41)))
m04$Model <- 'Common Var NB'
m04$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04" ,"path", "localp", "locals", "nbetastotal", "ncorre")])
##########################
# Model 5 - Negbin - Fixed Var
##########################
load("Negbin_var_fixed.RData")
b1p <- melt(estimates_negbin$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
r1p <- melt(estimates_negbin$rho)
d1p <- exp(melt(estimates_negbin$disp)) # Transformação aqui.
m05 <- data.frame(Estimate = unname(rbind(b1p, r1p, d1p)),
                  `Std. Error` = NA,
                  Z = NA,
                  `p value` = NA,
                  check.names = F,
                  row.names = c(row.names(m01)[1:nbetastotal], 
                                row.names(m01)[(nbetastotal+nr+1):(nbetastotal+nr+ncorre)], 
                                paste0('exp(phi)',1:41)))
m05$Model <- "Fixed Var NB"
m05$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05" ,"path", "localp", "locals", "nbetastotal", "ncorre")])
##########################
# Model 6 - Compoisson FULL
##########################
setwd(paste0(path, localp))
load("Compoisson_full.RData")
m06 <- as.data.frame(ss_sd$ests)
rownames(m06[(nbetastotal+nr+ncorre+1):nrow(m06),]) <- paste0('exp(nu)',1:41)
m06$Model <- "COM-Poisson"
m06$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06" ,"path", "localp", "locals", "nbetastotal", "ncorre")])
##########################
# Model 7 - Compoisson Disp Fixed 
##########################
setwd(paste0(path, locals))
load("compoisson_nu_fixed.RData")
# estimates_compoisson_s
m07 <- as.data.frame(estimates_compoisson_sd$ests)
m07$Model <- 'Fixed Disp COM-Poisson'
m07$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07" ,"path", "localp", "locals", "nbetastotal", "ncorre")])
##########################
# Model 8 - Compoisson Comum Var
##########################
load("compoisson_var_comum.RData")
m08 <- as.data.frame(estimates_compoisson_sd$ests)
m08$Model <- "Common Var COM-Poisson"
m08$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08","path", "localp", "locals", "nbetastotal", "ncorre")])
##########################
# Model 9 - Compoisson Fixed Var
##########################
load("compoisson_var_fixed.RData")
m09 <- as.data.frame(estimates_compoisson_sd$ests)
m09$Model <- "Fixed Var COM-Poisson"
m09$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09" ,"path", "localp", "locals", "nbetastotal", "ncorre")])
##########################
# Model 10 - Poisson Rho fixed
##########################
setwd(paste0(path, locals))
load("Poisson_rho_fixed.RData")
# estimates_compoisson_s
b1p <- melt(estimates_poisson_s1$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_poisson_s1$sigma)
m10 <- data.frame(Estimate = unname(rbind(b1p, s1p)),
                  `Std. Error` = NA,
                  Z = NA,
                  `p value` = NA,
                  check.names = F,
                  row.names = row.names(m01)[1:(nbetastotal+nr)])
m10$Model <- "Rho Zero Poisson"
m10$MainModel <- "Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09", "m10" ,"path", "localp", "locals", "nbetastotal", "ncorre")])
##########################
# Model 11 - NB Rho fixed
##########################
load("Negbin_rho_fixed.RData")
# estimates_compoisson_s
b1p <- melt(estimates_negbin_1$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_negbin_1$sigma)
d1p <- exp(melt(estimates_negbin_1$disp))
m11 <- data.frame(Estimate = unname(rbind(b1p, s1p, d1p)),
                  `Std. Error` = NA,
                  Z = NA,
                  `p value` = NA,
                  check.names = F,
                  row.names = c(row.names(m10), paste0('exp(phi)',1:nr)))
m11$Model <- "Rho Zero NB"
m11$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09", "m10", "m11" ,"path", "localp", "locals", "nbetastotal", "ncorre")])
##########################
# Model 12 - COM-Poisson rho fixed
##########################
load("Compoisson_rho_fixed.RData")
# estimates_compoisson_s
ls()
b1p <- melt(ss$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(ss$sigma)
d1p <- exp(melt(ss$disp))
m12 <- data.frame(Estimate = unname(rbind(b1p, s1p, d1p)),
                  `Std. Error` = NA,
                  Z = NA,
                  `p value` = NA,
                  check.names = F,
                  row.names = c(row.names(m10), paste0('exp(nu)',1:nr)))
m12$Model <- "Rho Zero COM-Poisson"
m12$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09", "m10", "m11", "m12" ,"path", "localp", "locals", "nbetastotal", "ncorre")])
#####################################################3
####### Joining all results
#####################################################3
# sample size 1281
m01$Parametros <- row.names(m01)
m02$Parametros <- row.names(m02)
m03$Parametros <- row.names(m03)
m04$Parametros <- row.names(m04)
m05$Parametros <- row.names(m05)
m06$Parametros <- row.names(m06)
m07$Parametros <- row.names(m07)
m08$Parametros <- row.names(m08)
m09$Parametros <- row.names(m09)
m10$Parametros <- row.names(m10)
m11$Parametros <- row.names(m11)
m12$Parametros <- row.names(m12)
mm <- rbind(m01,m02,m03,m04,m05,m06, m07, m08, m09, m10, m11, m12)
row.names(mm) <- NULL
# Creating variables
mm$Tipo <- ifelse(grepl("b[0-9]", mm$Parametros), "Média",
                  ifelse(grepl("sig[0-9]", mm$Parametros), "Variância",
                         ifelse(grepl("r_Y[0-9]", mm$Parametros), "Correlação",
                                ifelse(grepl("nu", mm$Parametros), "Dispersão",
                                       ifelse(grepl("phi", mm$Parametros), "Dispersão", "Deu Ruim")))))
mm <- mm[order(mm$MainModel),]
mm[(mm$Tipo=="Dispersão")&(mm$MainModel=="COM-Poisson"),"Parametros"] <- rep(paste0("exp(nu)", 1:41),4)


mm$Response <- ifelse(mm$Tipo=="Média", gsub(".*_Y", "", mm$Parametros),
                      ifelse(mm$Tipo=="Variância", gsub("sig", "", mm$Parametros),
                             ifelse(mm$Tipo=="Correlação", gsub("r_", "", mm$Parametros),
                                    ifelse(mm$Tipo=="Dispersão", gsub("[^0-9.-]", "", mm$Parametros), "BU"))))

varsnames <- c(" 1.Amblyopone", " 2.Aphaenogaster", " 3.Camponotus.Ci", " 4.Camponotus.Cl", 
               " 5.Camponotus.Co", " 6.Camponotus.Ni", " 7.Camponotus.Ns", " 8.Cardiocondyla", 
               " 9.Crematogaster", "10.Heteroponera", "11.Iridomyrmex.Bi", "12.Iridomyrmex.Dr", 
               "13.Iridomyrmex.Mj", "14.Iridomyrmex.Pu", "15.Iridomyrmex.Ru", 
               "16.Iridomyrmex.Su", "17.Iridomyrmex.Ss", "18.Melophorus.E", 
               "19.Melophorus.F", "20.Melophorus.H", "21.Meranoplus.A", "22.Monomorium.L", 
               "23.Monomorium.Ro", "24.Monomorium.Sy", "25.Myrmecia", "26.Notoncus.Ca", 
               "27.Notoncus.Ec", "28.Nylanderia", "29.Ochetellus", "30.Paraparatrechina", 
               "31.Pheidole.A", "32.Pheidole.B", "33.Pheidole.E", "34.Pheidole.J", 
               "35.Polyrhachis", "36.Rhytidoponera.A", "37.Rhytidoponera.B", 
               "38.Solenopsis", "39.Stigmacros", "40.Tapinoma", "41.Tetramorium")
mm$Response2 <- factor(mm$Response, 
                      levels = as.character(1:41),
                      labels = varsnames)



names(mm)[1:4] <- c("Estimativas", "ErroPadrao", "Z", "p value")
mm$NAN <- factor(as.character(ifelse(is.nan(mm$ErroPadrao), 1, 20)))
mm$ModelSpec <- ifelse(mm$Model%in%c("COM-Poisson", "NB", "Poisson"), "Full",
                       ifelse(mm$Model%in%c("Fixed Disp COM-Poisson", "Fixed Disp NB"), "Fixed Disp",
                              ifelse(mm$Model%in%c("Common Var COM-Poisson", "Common Var NB"), "Common Var",
                                     ifelse(mm$Model%in%c("Fixed Var COM-Poisson", "Fixed Var NB"), "Fixed Var",
                                            ifelse(mm$Model%in%c("Rho Zero Poisson", "Rho Zero NB", "Rho Zero COM-Poisson"), "Rho 0", "Deu Ruim")))))
#####################################################
####### Beta Graphic
#####################################################
mmmedia <- mm %>% 
  filter(Tipo=="Média") %>% 
  mutate(Parametros = factor(gsub("_Y\\d+", "", Parametros), 
                             levels = paste0("b",0:5),
                             labels = c("Intercepto", "Bare.ground", "Canopy.cover", "Shrub.cover", "Feral.mammal.dung", "Volume.lying.CWD")))
inf <- 1
sup <- 8
plot.mean <- function(data, inf, sup){
  data %>% 
    filter(between(Response, inf, sup)) %>% 
    ggplot() +
    geom_pointrange(aes(x = Parametros, y = Estimativas,
                        ymin = Estimativas-1.96*ErroPadrao,
                        ymax = Estimativas+1.96*ErroPadrao,
                        linetype = Model,
                        shape = NAN,
                        color = ModelSpec),
                    position = position_dodge(width = 0.6),
                    alpha = .8) +
    scale_shape_manual(values = c(1,20),
                       labels = c("Without SE", "With SE")) +
    scale_linetype_manual(values=1:12) +
    guides(linetype = guide_legend(override.aes = list(shape = NA,
                                                       direction = 'vertical'),
                                   ncol = 3),
           shape = guide_legend(override.aes = list(linetype = "blank"),
                                ncol = 1),
           color = guide_legend(ncol = 1)) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          text = element_text(size = 17),
          legend.position = "bottom",
          legend.title = element_blank()) +
    facet_wrap(~Response2, scales = "free_x", ncol = 4) +
    # coord_cartesian(ylim = c(-5,5)) +
    coord_flip(ylim = c(-20,20)) +
    labs(y = TeX('$\\hat{\\beta}\\pm{1.96 SE}$'),
         linetype = "Model")
}
pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
plot.mean(mmmedia, 1, 8)
ggsave(paste0(pathfig, "/ant_beta_initial_model_1_8.pdf"), width = 36, height = 36, units = "cm")
plot.mean(mmmedia, 9, 16)
ggsave(paste0(pathfig, "/ant_beta_initial_model_9_16.pdf"), width = 36, height = 36, units = "cm")
plot.mean(mmmedia, 17, 24)
ggsave(paste0(pathfig, "/ant_beta_initial_model_17_24.pdf"), width = 36, height = 36, units = "cm")
plot.mean(mmmedia, 25, 32)
ggsave(paste0(pathfig, "/ant_beta_initial_model_25_32.pdf"), width = 36, height = 36, units = "cm")
plot.mean(mmmedia, 33, 41)
ggsave(paste0(pathfig, "/ant_beta_initial_model_33_41.pdf"), width = 40, height = 40, units = "cm")

df <- expand.grid(X1 = 1:10, X2 = 1:10)
df$value <- df$X1 * df$X2

#####################################################
####### Dispersion
#####################################################
myHeader1 <- c(1,4,4)
names(myHeader1) <- c(" ", "COM-Poisson ($\\\\nu$)", "NB ($\\\\phi$)")

myHeader0 <- c(1,1,2,1,1,2,1)
names(myHeader0) <- c(" ","Full", "Variance", "Rho","Full", "Variance", "Rho")

###########
# Estimates
###########
mm %>% 
  filter(Tipo%in%c("Dispersão")) %>% 
  select(-Z, -`p value`, -MainModel, -Tipo, -Parametros, -NAN, -ModelSpec, -Response, -ErroPadrao) %>% 
  pivot_wider(names_from = Model,
              values_from = c("Estimativas")) %>% 
  # select(Model, Estimativas_Nmsp, ErroPadrao_Nmsp,
         # Estimativas_Nmosp, ErroPadrao_Nmosp, everything()) %>% 
  kable(caption = toupper("Dispersion parameter estimates for each initial model and outcome of ANT data"),
        label = "antdispinit",
        booktabs = T,
        digits = 2,
        row.names = F,
        escape = F,
        format = "latex",
        col.names = c("Response", rep(c("", "Common", "Fixed", "0"),2)),
        align = "c") %>% 
  add_header_above(header = myHeader0,
                   escape = F) %>% 
  add_header_above(header = myHeader1,
                   escape = F)

###########
# SE
###########
mm %>% 
  filter(Tipo%in%c("Dispersão")) %>% 
  select(-Z, -`p value`, -MainModel, -Tipo, -Parametros, -NAN, -ModelSpec, -Response, -Estimativas) %>% 
  pivot_wider(names_from = Model,
              values_from = c("ErroPadrao")) %>% 
  kable(caption = toupper("Dispersion parameter standard error (SE) estimates for each initial model and outcome of ANT data"),
        label = "antdispinit_se",
        booktabs = T,
        digits = 2,
        row.names = F,
        escape = F,
        format = "latex",
        col.names = c("Response", rep(c("", "Common", "Fixed", "0"),2)),
        align = "c") %>% 
  add_header_above(header = myHeader0,
                   escape = F) %>% 
  add_header_above(header = myHeader1,
                   escape = F) %>% 
  add_footnote(label = 'Empty cells comprises of non returned SE')


#####################################################
####### Variância
#####################################################
plot.var <- function(data, inf, sup){
  data %>% 
    filter(Tipo=="Variância") %>% 
    filter(!ModelSpec=='Common Var') %>% 
    filter(between(Response, inf, sup)) %>% 
    ggplot(aes(x = Model, y = Estimativas)) +
    geom_pointrange(aes(ymin = Estimativas-1.96*ErroPadrao,
                        ymax = Estimativas+1.96*ErroPadrao,
                        # linetype = Model,
                        shape = NAN,
                        color = ModelSpec),
                    position = position_dodge(width = 0.3),
                    alpha = .8) +
    scale_shape_manual(values = c(1,20),
                       labels = c("Without SE", "With SE")) +
    guides(shape = guide_legend(override.aes = list(linetype = "blank"),
                                nrow = 1),
           color = guide_legend(nrow = 1)) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          text = element_text(size = 16),
          legend.position = "bottom",
          legend.title = element_blank()) +
    facet_wrap(~Response2, scales = "free_x", ncol = 5) +
    coord_flip() +
    labs(y = TeX('$\\hat{\\sigma}\\pm {1.96 SE}$'),
         linetype = "Model")
}
pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
plot.var(mm, 1, 21)
ggsave(paste0(pathfig, "/nhanes_var_initial_model_1_21.pdf"), width = 31, height = 32, units = "cm")
plot.var(mm, 22, 41)
ggsave(paste0(pathfig, "/nhanes_var_initial_model_22_41.pdf"), width = 31, height = 30, units = "cm")


mm %>% 
  filter(Tipo=="Variância") %>% 
  filter(ModelSpec=='Common Var') %>% 
  select(Model, Estimativas, ErroPadrao) %>% 
  mutate(Ic.Inf = Estimativas-1.96*ErroPadrao,
         Ic.Sup = Estimativas+1.96*ErroPadrao) %>% 
  kable(digits = 4, align = c('l', rep('c', 4)),
        format = 'latex',
        label = 'antvarcommon',
        escape = F,
        booktabs = T,
        col.names = c('Model', 'Estimates', 'SE', 'Lower CI', 'Upper CI'),
        caption = "STANDARD DEVIATION ESTIMATES OF RANDOM EFFECTS AND 95% CONFIDENCE INTERVALS FOR MODELS WITH COMMON VARIANCE FOR ANT DATA") %>% 
  add_footnote(label = 'Empty cells comprises of non returned SE')

#####################################################
####### Correlation
#####################################################

dgcorr <- mm %>% 
  filter(Tipo == "Correlação")
dgcorr$Resp <- gsub("r_", "", dgcorr$Parametros)
dgcorr$Resp <- gsub("Y", "", dgcorr$Resp)
dgcorr <- dgcorr[, -c(6,7,8)]
resps2 <- stringr::str_split_fixed(dgcorr$Resp, "_", 2)
colnames(resps2) <- c("Var1", "Var2")
dgcorr <- cbind(dgcorr, resps2)
dgcorr$Var1 <- as.numeric(dgcorr$Var1)
dgcorr$Var2 <- as.numeric(dgcorr$Var2)
dgcorr$IcInf <- dgcorr$Estimativas-1.96*dgcorr$ErroPadrao
dgcorr$IcSup <- dgcorr$Estimativas+1.96*dgcorr$ErroPadrao
dgcorr$Valorp <- pnorm(abs(dgcorr$Estimativas/dgcorr$ErroPadrao), lower.tail = F)*2


# Agora operação para colocar em formato de matrix
df_to_mat <- function(dd){
  mat_est <- matrix(numeric(), nrow = 41, ncol = 41)
  mat_se <- matrix(numeric(), nrow = 41, ncol = 41)
  mat_inf <- matrix(numeric(), nrow = 41, ncol = 41)
  mat_sup <- matrix(numeric(), nrow = 41, ncol = 41)
  mat_p <- matrix(numeric(), nrow = 41, ncol = 41)
  i <- 1
  j <- 2
  for (i in 1:41){
    for (j in i:41){
      if (i == j){
        mat_est[i,j] <- 1
        mat_se[i,j] <- 0
      } else{
        dgcorrfilt <- dd[dd$Var1 %in% i & dd$Var2 %in% j, , drop = F]
        mat_est[i, j] <- dgcorrfilt[, 1]
        mat_se[i, j] <- dgcorrfilt[, 2]
        mat_inf[i, j] <- as.numeric(ifelse(dgcorrfilt[, 13]>1, 1, ifelse(dgcorrfilt[, 13]<(-1),-1, dgcorrfilt[, 13])))
        mat_sup[i, j] <- as.numeric(ifelse(dgcorrfilt[, 14]>1, 1, ifelse(dgcorrfilt[, 14]<(-1),-1, dgcorrfilt[, 14])))
        mat_p[i, j] <- dgcorrfilt[, 4]
      }
    }
  }
  ll <- list(mat_est, mat_se, mat_inf, mat_sup, mat_p)
  names(ll) <- c("est", "se", "inf", "sup", "p")
  return(ll)
}
library(corrplot)
gera_grafico <- function(data, namefile, pathfig){
  dgcorrP <- data
  corrp <- df_to_mat(dgcorrP)
  pdf(file = paste0(pathfig, namefile), width = 6, height = 5.5)
  graph1 <- corrplot(
                     corrp$est
                     # mar = c(0,0,0,0), 
                     ,tl.cex = 0.6 
                     ,cl.cex = .7
                     # pch = "*",
                     ,type = "upper"
                     ,diag = F
                     ,outline = F
                     # ,pch.cex = 1.7       # Controls the size of the star
                     # p.mat = 1-corrp$p, #Opposite
                     # ,sig.level = .95
                     )   #Opposite
  # print(graph1)
  dev.off()
}

gera_grafico2 <- function(data, namefile, pathfig){
  dgcorrP <- data
  corrp <- df_to_mat(dgcorrP)
  pdf(file = paste0(pathfig, namefile), width = 6, height = 5.5)
  graph1 <- corrplot(
    corrp$est
    # mar = c(0,0,0,0), 
    ,tl.cex = 0.6 
    ,cl.cex = .7
    ,pch = "*"
    ,type = "upper"
    ,diag = F
    ,outline = F
    ,pch.cex = 1.7       # Controls the size of the star
    ,p.mat = 1-corrp$p #Opposite
    ,sig.level = .95
  )   #Opposite
  # print(graph1)
  dev.off()
}

# Decidi excluir as estrelas de significância, pois apenas 2 modelos tinham significância estatística.
models <- unique(dgcorr$Model)
for (i in 1:length(models)) {
  nm <- models[i]
  gera_grafico(dgcorr[dgcorr$Model==nm, ], paste0("/cor_ant_",nm,"_initial.pdf"), pathfig)
}

# Refaz só para aquelas que tem coef significativo (warning é normal)
models2 <- unique(dgcorr$Model)[2:4]
for (i in 1:length(models2)) {
  nm <- models2[i]
  gera_grafico2(dgcorr[dgcorr$Model==nm, ], paste0("/cor_ant_",nm,"_initial.pdf"), pathfig)
}

# # Teste
# mmpoisson <- dgcorr[dgcorr$Model=="Fixed Disp COM-Poisson" , ]
# 
# mmpoisson[mmpoisson$`p value`<0.05,]
# 
# corrp <- df_to_mat(mmpoisson)
# graph1 <- corrplot(
#   corrp$est
#   # mar = c(0,0,0,0), 
#   ,tl.cex = 0.6 
#   ,cl.cex = .7
#   ,pch = "*"
#   ,type = "upper"
#   ,diag = F
#   ,outline = F
#   ,pch.cex = 1.7       # Controls the size of the star
#   ,p.mat = 1-corrp$p #Opposite
#   ,sig.level = .95
# )   #Opposite
# 
# 
# 
# gera_grafico(dgcorr[dgcorr$Model=="Fixed Disp COM-Poisson", ], "/cor_ant_poisson_initial.pdf", pathfig)
# gera_grafico(dgcorr[dgcorr$Model=="NB", ], "/cor_ant_NB.pdf", pathfig)
# gera_grafico(dgcorr[dgcorr$Model=="COM-Poisson", ], "/cor_ant_cmp.pdf", pathfig)
# 
# sum(dgcorr[dgcorr$Distribuição=="COM-Poisson", ]$Valorp<.05, na.rm = T)
# sum(dgcorr[dgcorr$Distribuição=="NB", ]$Valorp<.05, na.rm = T)
# sum(dgcorr[dgcorr$Distribuição=="Poisson", ]$Valorp<.05, na.rm = T)
# 
# 






