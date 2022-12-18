# https://www.dell.com/community/Notebooks-Laptops/Problemas-no-Dell-G3-e-G7/td-p/6243192
rm(list=ls())
library(mcglm);library(ggplot2)
library(dplyr);library(reshape2)
library(facetscales);library(knitr)
library(kableExtra);library(GGally)
# Vars Description #############################################################
data.des <- data.frame(Name = c("sex", "age", "income", "levyplus", "freepoor", "freerepa", "illnes", "actdays", "hscore", "chcond"),
           Description = c("Factor with levels male and female",
                           "Respondent's age in years divided by 100",
                           "Respondent's annual income in Australian dollars divided by 1000",
                           "Coded factor. If respondent is covered by private health insurance fund for private patients in public hospital with doctor of choice (1) or otherwise (0)",
                           "Coded factor. If respondent is covered by government because low income, recent immigrant, unemployed (1) or otherwise (0)",
                           "Coded factor. If respondent is covered free by government because of old",
                           "Number of illnesses in past 2 weeks, with 5 or more illnesses coded as 5",
                           "Number of days of reduced activity in the past two weeks due to illness or injury",
                           "hscore & Respondent's general health questionnaire score using Goldberg's method. High score indicates poor health",
                           "Factor with three levels. If respondent has chronic condition(s) and is limited in activity (limited), or if the respondent has chronic condition(s) but is not limited in activity (nonlimited) or otherwise (otherwise, reference level)"
                           ))
kable(data.des, label = "varsAHS" ,format = "latex", booktabs = T, digits = 3, escape = T,
      row.names = F,
      align = c("l"),
      caption = "COVARIATES COLLECTED IN THE AUSTRALIA HEALTH SURVEY (AHS)") %>% 
  kable_styling() %>% 
  column_spec(2, width = "11.5cm")


# Code #########################################################################
responses <- names(ahs)[grep("N", names(ahs))]
covariates <- names(ahs)[1:10]
long <- melt(ahs[, responses])
ahs$levyplus <- factor(ahs$levyplus, labels = c("No", "Doctor of Choice"))
ahs$freepoor <- factor(ahs$freepoor, labels = c("Not covered", "Covered gov"))
ahs$freerepa <- factor(ahs$freerepa, labels = c("Not", "Covered old"))
# Barplot ###########################################################################
scales_y <- list(
  `Ndoc` = scale_x_continuous(breaks = seq(0, 9, 2)),
  `Nndoc` = scale_x_continuous(limits = c(0, 11), breaks = seq(0, 12, 2)),
  `Nadm` = scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)),
  `Nhosp` = scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)),
  `Nmed` = scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, 2))
)
# pdf(file = "desc_ahs.pdf", height = 2)
ggplot(data = long,
       aes(x = value)) +
  geom_bar(stat = "count") +
  theme_bw() +
  labs(x = "Occurrence",
       y = "Frequency") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 9)) +
  facet_grid_sc(cols = vars(variable), 
                scales = list(x = scales_y))

# facet_wrap(~variable, scales = "free", nrow = 1) +
# scale_x_continuous(breaks = 0:10,
# labels = 0:10)
dev.off()
ggsave("desc_ahs.pdf", width = 15, height = 5, units = "cm")

# Descriptive statistics ###########################################################################
# GDI
source("~/Dropbox/SuperDisperso_AJS/Code/PaperCode.R")
# GDI = 17.94421 
gdi <- as.numeric(GDI(ahs[, responses]))
# Standard Error = 1.986486
gdi_se <- as.numeric(sqrt(sigma_GDI(ahs[, responses])/nrow(ahs)))

corre <- cor(ahs[, responses], method = "spearman")
corre[lower.tri(corre, diag = T)] <- NA
mat <- cbind(corre, 
             "Mean" = t(t(sapply(ahs[, responses], mean))),
             "Var" = t(t(sapply(ahs[, responses], var))))
mat <- cbind(mat, 
             "DI" = mat[,7]/mat[,6])
mat2 <- cbind(corre, 
              "Mean" = t(t(sapply(ahs[, responses], mean))))
mat <- round(mat, 3)
mat <- as.data.frame(mat)
mat$GDI <- c(NA,NA,"17.944 (1.99)",NA,NA)
mat <- mat[,-1]
colnames(mat)[6:8] <- ""
options(knitr.kable.NA = '')
kable(mat, format = "latex", booktabs = T, digits = 3, escape = F,
      align = c("c"),
      label = "descAHS",
      caption = "DESCRIPTIVE MEASUREMENTS FOR AHS RESPONSE VARIABLES") %>% 
  kable_styling() %>%
  add_header_above(c(" "=1,"Spearman Correlation $\\\\rho$"=4,"Mean"=1,"Variance"=1, "DI" = 1, "GDI(SE)"),
                   escape = F)

kable(mat2, format = "latex", booktabs = T, digits = 3, escape = F,
      align = c("c")) %>% 
  kable_styling(latex_options = "scale_down", position = "center", font_size = 6, full_width = F) %>%
  add_header_above(c(" "=1,"Correlação de Spearman $\\\\rho$"=5,"Média"=1),
                   escape = F)

#######################
# Attempts/Bad Charts
#######################

# Scaterplot Matrix ############################################################
ahs[, responses] <- lapply(ahs[, responses], as.numeric)
corre <- cor(ahs[, responses], method = "spearman")
library(corrplot)
whiteblack <- c("white", "black")
png(file = "cor_ahs.png", units = "cm", width = 20, height = 20, res = 150)
corrplot::corrplot(corre, mar = c(0,0,1,0), 
                   method = c("circle"),
                   addCoef.col = "black",
                   type = "upper",
                   diag = F,
                   outline = T,
                   tl.cex = 1.4, cl.cex = 1)
# corrplot(corre, col = whiteblack, bg = "gold2")
dev.off()

# ggpairs(ahs[, responses],
#         diag  = list(continuous = "barDiag"),
#         upper = list(continuous = "cor"))
# 
# PerformanceAnalytics::chart.Correlation(ahs[, responses])
# Response vs Covariates #######################################################
covariates_fac <- covariates[sapply(ahs[, covariates], is.factor)]
covariates_num <- covariates[sapply(ahs[, covariates], is.numeric)]
# Numeric
long_num <- melt(ahs[, c(covariates_num, "Ndoc")], id.vars = "Ndoc")
# long_num <- melt(ahs[, c("actdays", "Ndoc")], id.vars = "Ndoc")
g1 <- ggplot(data = long_num, aes(x = log(value+.01), y = Ndoc)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE, n = 20, span = .4) +
  facet_wrap(~variable, scales = "free_x", nrow = 1) +
  scale_y_continuous(labels = 0:9,
                     breaks = 0:9);g1

# Factor
long_fac <- melt(ahs[, c(covariates_fac, "Ndoc")], id.vars = "Ndoc")
g2 <- ggplot(data = long_fac, aes(x = value, y = Ndoc)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_x", nrow = 1) +
  scale_y_continuous(labels = 0:9,
                     breaks = 0:9)
gridExtra::grid.arrange(g1,g2)
sapply(ahs[, responses], table)
GGally::ggbivariate(ahs[, -c(12:15)],
                    outcome = "Ndoc",
                    explanatory = NULL)

suppressWarnings(p_(pm))
# Under, Equi and Over dispersion study ############################################################
data_compare <- data.frame(eq = rpois(1000, 8),
                           sup = rnbinom(1000, mu = 8, size = 4),
                           sub = rcomp(1000, mu = 8, nu = 8)) # nu = 1 is poisson, nu > 1 is subdispersion, and nu < 1 is overdispersion )
data_compare <- data.frame(eq = rpois(1000, 1.5),
                           sup = rnbinom(1000, mu = 1.33, size = .05),
                           sup2 = rnbinom(1000, mu = 1.2, size = 1.2)
                           # sub = rcomp(1000, mu = 8, nu = 8)
                           ) # nu = 1 is poisson, nu > 1 is subdispersion, and nu < 1 is overdispersion )
long <- melt(data_compare)
ggplot(data = long,
       aes(x = value)) +
  geom_bar(stat = "count") +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  theme_bw() +
  labs(x = "Occurrence",
       y = "Frequency")
# Conclusion
# For a overdispersion data, the barplot looks like a Gamma data
# For equidispersed data, it looks like normal
# For under dispersion data, it looks like normal, but only a small range of values are possible (if compared to equidispersed data)
