# Website
# https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Questionnaire&CycleBeginYear=2007
############ SEQN is the ID
############ Sexual Behavior table (Response)
# https://wwwn.cdc.gov/Nchs/Nhanes/2007-2008/SXQ_E.htm#SXQ470
setwd("/home/guilherme/Dropbox/Underdispersed_Count/Code/Descriptive_Analysis")
pkg <- c("dplyr", "haven", "ggplot2", "reshape2")
sapply(pkg, require, character.only = T)
dsexual <- read_xpt(file = "SXQ_E.XPT")
vars <- c("SEQN","SXQ450", "SXQ470", "SXQ590")
# SXQ450 - # of male sex partners/year            (700 strange)                                                  # In the past 12 months, with how many males have you had vaginal, anal, or oral sex?
# SXQ470 - # of male oral sex partners/year       (700 strange 99999 excluded)                                   # With how many of these males have you had only oral sex?
# SXQ590 - # of sex partners five years older/year in the past year (77777 & 99999 excluded) (It seems correct)  # Of the persons you had sex with in the past 12 months, how many were five or more years older than you?
dsexual <- dsexual[, vars]
dsexual <- na.omit(dsexual) #1293 (1280)
dsexual <- dsexual[dsexual$SXQ450 < 700 & !is.na(dsexual$SXQ450), ] #1292
dsexual <- dsexual[dsexual$SXQ470 < 700 & !is.na(dsexual$SXQ470), ] #1290
dsexual <- dsexual[dsexual$SXQ590 < 7 & !is.na(dsexual$SXQ590), ]   #1281
responses <- c("Nmsp", "Nmosp", "Nspfy")
names(dsexual) <- c("id", responses)
############ Demographic Variables & Sample Weights
# https://wwwn.cdc.gov/Nchs/Nhanes/Search/DataPage.aspx?Component=Demographics&CycleBeginYear=2007
ddemo <- read_xpt(file = "DEMO_E.XPT")
vars <- c("SEQN","RIAGENDR", "RIDRETH1", "DMDEDUC2", "DMDMARTL", "DMDHRAGE")
ddemo <- ddemo[, vars]
# ddemo <- na.omit(ddemo)
ddemo <- na.omit(ddemo)
ddemo <- ddemo[ddemo$RIAGENDR==2, ]
ddemo <- ddemo[ddemo$DMDEDUC2<=5, ]
ddemo <- ddemo[, -2]
names(ddemo) <- c("id", "Race", "Education", "Marital", "Age")
ddemo$Race <- factor(ifelse(ddemo$Race==3, "White", "Others"))
ddemo$Marital <- factor(ifelse(ddemo$Marital==1, "Married", "Others"),
                        levels = c("Others", "Married"))
dfinal <- merge(dsexual, ddemo, by = "id", all.x = T, sort = F)
write.table(dfinal, file = "nhanes.txt")
# Deleting cases where at least one of the variables is missing, we have a total of
# 1280 female respondents whose age range from 20 to 59 years. The average age of these
# respondents is 38.31 years. The data consist of 43% white and 57.5% married women.
# https://www.tandfonline.com/author/Famoye%2C+Felix
nrow(dfinal)
nrow(ddemo)
nrow(dsexual)
sapply(dsexual, table)
sapply(ddemo, table)
# The predictor variables considered for this illustration are
# race (1 = white, 0 = others), 
# educational level (range from 1 = below 9th grade to 5 = college graduate), (numeric)
# marital status (1 = married, 0 = others), 
# and age.
