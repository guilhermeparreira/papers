rm(list = ls())
setwd("/home/guilherme/Dropbox/Underdispersed_Count/Code/Descriptive_Analysis")
pkg <- c("dplyr", "haven", "facetscales", "ggplot2", "reshape2", "kableExtra")
sapply(pkg, require, character.only = T)
dfinal <- read.table("nhanes.txt", header = T)
responses <- c("Nmsp", "Nmosp", "Nspfy")
##################
# Stats to placed in the text
##################
nrow(dfinal)     #against n    = 1280
range(dfinal$Age)
mean(dfinal$Age) #mean = 38,31
prop.table(table(dfinal$Marital))
prop.table(table(dfinal$Race))
prop.table(table(dfinal$Education))

prop.table(table(dfinal$Marital, dfinal$Nmsp), 1)

##################
# Multivariate Dispersion Index
##################
source("PaperCode.R")
gdi <- GDI(dfinal[, responses]) # = 1.092186 Almost equi dispersed
# Standard Error = .2201552
gdi_se <- sqrt(sigma_GDI(dfinal[, responses])/nrow(dfinal))

##################
# Mean/Var and correlation table  
##################
corre <- cor(dfinal[, responses], method = "spearman")
corre[lower.tri(corre, diag = T)] <- NA
mat <- cbind(corre, 
             "Mean" = t(t(sapply(dfinal[, responses], mean))),
             "Var" = t(t(sapply(dfinal[, responses], var))))
mat <- cbind(mat, 
             "DI" = mat[,5]/mat[,4])
mat <- as.data.frame(mat)
mat$GDI <- c(NA, "1.092(.22)", NA)
colnames(mat)[4:7] <- c("")

options(knitr.kable.NA = '')
kable(mat, format = "latex", booktabs = T, digits = 3, escape = F,
      align = c("c"),
      label = "desc_nhanes",
      caption = "MEAN, VARIANCE, DI AND CORRELATION FOR NHANES RESPONSE VARIABLES") %>% 
  kable_styling() %>%
  add_header_above(c(" "=1,"Spearman Correlation $\\\\rho$"=3,
                     "Mean"=1, "Variance" = 1, "DI"=1, "GDI(SE)"=1),
                   escape = F)
##################
# Barplot (Frequency) 
##################

long <- melt(dfinal[, responses])
ggplot(data = long,
       aes(x = value)) +
  geom_bar(stat = "count") +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        strip.text.x = element_text(size = 6)) +
  labs(x = "Occurrence",
       y = "Frequency")


#################
# Check of model resuls
#################
library(dplyr)
dfinal %>% 
  group_by(Marital) %>% 
  summarise(across(c(Nmsp, Nmosp, Nspfy), list(mean)))
dfinal %>% 
  group_by(Race) %>% 
  summarise(across(c(Nmsp, Nmosp, Nspfy), list(mean)))
dfinal %>% 
  group_by(Marital) %>% 
  summarise(across(c(Nspfy), list(mean)))

coef(glm(Nmsp ~ Race + Marital, data = dfinal, family = poisson(link = "log")))
coef(glm(Nmosp ~ Race + Marital, data = dfinal, family = poisson(link = "log")))
coef(glm(Nspfy ~ Race + Marital + Education, data = dfinal, family = poisson(link = "log")))


coef(glm(Nmsp ~ Race, data = dfinal, family = poisson(link = "log")))
coef(glm(Nmosp ~ Race, data = dfinal, family = poisson(link = "log")))
coef(glm(Nspfy ~ Race, data = dfinal, family = poisson(link = "log")))

coef(glm(Nmsp ~ Marital, data = dfinal, family = poisson(link = "log")))
coef(glm(Nmosp ~ Marital, data = dfinal, family = poisson(link = "log")))
coef(glm(Nspfy ~ Marital, data = dfinal, family = poisson(link = "log")))

coef(glm(Nspfy ~ Education, data = dfinal, family = poisson(link = "log")))




# scales_y <- list(
#   `Ndoc` = scale_x_continuous(breaks = seq(0, 9, 2)),
#   `Nndoc` = scale_x_continuous(limits = c(0, 11), breaks = seq(0, 12, 2)),
#   `Nadm` = scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)),
# )
# facet_grid_sc(cols = vars(variable), 
#               scales = list(x = scales_y))
