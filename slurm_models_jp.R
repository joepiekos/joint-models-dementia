library(tidyverse)
library(survival)
library(lubridate)
library(GLMMadaptive)
library(JMbayes2)
library(mice)
library(nlme)
library(splitTools)

set.seed(123)

#this is the file run via SLURM with the stats department cluster

ORCHARD_2admits_train <- read_csv("ORCHARD_2admits_train.csv")
ORCHARD_2admits_test <- read_csv("ORCHARD_2admits_test.csv")
ORCHARD_2admits_train$Age.x <- scale(ORCHARD_2admits_train$Age.x)
ORCHARD_2admits_test$Age.x <- scale(ORCHARD_2admits_test$Age.x)

#these models removed admissions within 3 months of the first admission, not used in project
# ORCHARD_2admits_3months_train <- read_csv("ORCHARD_2admits_3months_train.csv")
# ORCHARD_2admits_3months_test <- read_csv("ORCHARD_2admits_3months_test.csv")
# ORCHARD_2admits_3months_train$Age.x <- scale(ORCHARD_2admits_3months_train$Age.x)
# ORCHARD_2admits_3months_test$Age.x <- scale(ORCHARD_2admits_3months_test$Age.x)


print("loading longitudinal models")

#if longitudinal models already fit then uncomment these

# HFRS.mixmod9 <- readRDS("HFRS.rds")
# CCI.mixmod1 <- readRDS("CCI.rds")
# infec.mixmod1 <- readRDS("infec.rds")
# delirium.mixmod1 <- readRDS("delirium.rds")
# SEND.mixmod1 <- readRDS("SEND.rds")
# SIRS.mixmod1 <- readRDS("SIRS.rds")
# Falls.mixmod1 <- readRDS("Falls.rds")


ORCHARD_2admits_train$HFRS.x <- round(ORCHARD_2admits_train$HFRS.x)
ORCHARD_2admits_test$HFRS.x <- round(ORCHARD_2admits_test$HFRS.x)
ORCHARD_2admits_train$infection_d1.x <- as.factor(ORCHARD_2admits_train$infection_d1.x)
ORCHARD_2admits_test$infection_d1.x <- as.factor(ORCHARD_2admits_test$infection_d1.x)
ORCHARD_2admits_train$has_delirium.x <- as.factor(ORCHARD_2admits_train$has_delirium.x)
ORCHARD_2admits_test$has_delirium.x <- as.factor(ORCHARD_2admits_test$has_delirium.x)
ORCHARD_2admits_train$Fall_history.x <- as.factor(ORCHARD_2admits_train$Fall_history.x)
ORCHARD_2admits_test$Fall_history.x <- as.factor(ORCHARD_2admits_test$Fall_history.x)

HFRS.mixmod9 <- mixed_model(HFRS.x ~ lower, data = ORCHARD_2admits_train, random = ~ 1 + lower|id, family = poisson(), verbose = TRUE, control = list(iter_EM=0))
# HFRS.mixmod9_3months <- mixed_model(HFRS.x ~ lower, data = ORCHARD_2admits_3months_train, random = ~ 1 + lower|id, family = poisson(), verbose = TRUE, control = list(iter_EM=0))
saveRDS(HFRS.mixmod9,file = "HFRS.rds")
# saveRDS(HFRS.mixmod9_3months,file = "HFRS_3months.rds")
CCI.mixmod1 <- mixed_model(CCI.x ~ lower, data = ORCHARD_2admits_train, random = ~ 1 + lower|id, family = GLMMadaptive::negative.binomial(), verbose = TRUE, control = list(iter_EM=0))
saveRDS(CCI.mixmod1, file = "CCI.rds")
infec.mixmod1 <- mixed_model(infection_d1.x ~ lower, data = ORCHARD_2admits_train, random = ~ 1 + lower|id, family = binomial(), verbose = TRUE, control = list(iter_EM=0))
# infec.mixmod1_3months <- mixed_model(infection_d1.x ~ lower, data = ORCHARD_2admits_3months_train, random = ~ 1 + lower|id, family = binomial(), verbose = TRUE, control = list(iter_EM=0))
saveRDS(infec.mixmod1,file = "infec.rds")
# saveRDS(infec.mixmod1_3months,file = "infec_3months.rds")
delirium.mixmod1 <- mixed_model(has_delirium.x ~ lower, data = ORCHARD_2admits_train, random = ~ 1 + lower|id, family = binomial(), verbose = TRUE, control = list(iter_EM=0))
# delirium.mixmod1_3months <- mixed_model(has_delirium.x ~ lower, data = ORCHARD_2admits_3months_train, random = ~ 1 + lower|id, family = binomial(), verbose = TRUE, control = list(iter_EM=0))
saveRDS(delirium.mixmod1, file = "delirium.rds")
# saveRDS(delirium.mixmod1_3months, file = "delirium_3months.rds")
SEND.mixmod1 <- mixed_model(SEND_TT_score.x ~ lower, data = ORCHARD_2admits_train, random = ~ 1 + lower|id, family = poisson(), verbose = TRUE, control = list(iter_EM=0))
saveRDS(SEND.mixmod1, file = "SEND.rds")
SIRS.mixmod1 <- mixed_model(SIRS.x ~ lower, data = ORCHARD_2admits_train, random = ~ 1 + lower|id, family = poisson(), verbose = TRUE, control = list(iter_EM=0))
saveRDS(SIRS.mixmod1, file = "SIRS.rds")
Falls.mixmod1 <- mixed_model(Fall_history.x ~ lower, data = ORCHARD_2admits_train, random = ~ 1 + lower|id, family = binomial(), verbose = TRUE, control = list(iter_EM=0))
saveRDS(Falls.mixmod1, file = "Falls.rds")

print("loaded longitudinal models successfully")

print("fitting cox model")

surv_ORCHARD_2admits_train <- ORCHARD_2admits_train[!duplicated(ORCHARD_2admits_train$id),]
# surv_ORCHARD_2admits_3months_train <- ORCHARD_2admits_3months_train[!duplicated(ORCHARD_2admits_3months_train$id),]

surv_ORCHARD_2admits_test <- ORCHARD_2admits_test[!duplicated(ORCHARD_2admits_test$id),]
# surv_ORCHARD_2admits_3months_test <- ORCHARD_2admits_3months_test[!duplicated(ORCHARD_2admits_3months_test$id),]

#interval censored parametric survival models can't handle values of zero, package issue. so change values of zero to very small values.
#surv_ORCHARD_2admits_train$lower[surv_ORCHARD_2admits_train$lower == 0] <- 0.00001
#surv_ORCHARD_2admits_test$lower[surv_ORCHARD_2admits_test$lower == 0] <- 0.00001
#surv_ORCHARD_2admits_train$upper[surv_ORCHARD_2admits_train$upper == 0] <- 0.00001
#surv_ORCHARD_2admits_test$upper[surv_ORCHARD_2admits_test$upper == 0] <- 0.00001
#cox_model_interval <- survreg(Surv(time = Lower ,time2 = Upper, event = status_interval, type = "interval") ~ Age.x + Sex.x + dem_IMD.x, data = surv_ORCHARD_2admits_train, x = TRUE,
 #                             model = TRUE, dist = "weibull",control = list(maxiter=500))

cox_model_right <- coxph(Surv(time = survtime, event = status_rightcens, type = "right") ~ Age.x + Sex.x + dem_IMD.x, data = surv_ORCHARD_2admits_train, x = TRUE, model = TRUE)
# cox_model_right_3months <- coxph(Surv(time = survtime, event = status_rightcens, type = "right") ~ Age.x + Sex.x + dem_IMD.x, data = surv_ORCHARD_2admits_3months_train, x = TRUE, model = TRUE)

print("cox done")

print("fitting joint model")
# jointFit_rcens <- jm(cox_model_right, list(HFRS.mixmod9,infec.mixmod1,delirium.mixmod1), time_var = "lower", data_Surv = surv_ORCHARD_2admits_train, id_var = "id", 
#                control = list(n_chains = 1, cores = 1, n_iter = 10000, n_burnin=500), priors = list(penalty_alphas = "horseshoe",penalty_gammas = "horeshoe"))
jointFit_rcens_3months <- jm(cox_model_right_3months, list(HFRS.mixmod9_3months,infec.mixmod1_3months,delirium.mixmod1_3months), time_var = "lower", data_Surv = surv_ORCHARD_2admits_3months_train, id_var = "id", 
                control = list(n_chains = 1, cores = 1, n_iter = 3500, n_burnin=500), priors = list(penalty_alphas = "horseshoe",penalty_gammas = "horeshoe"))

saveRDS(jointFit_rcens, file = "jointFit_rcens_3months.rds")

# jointFit_int <- jm(cox_model_interval, list(HFRS.mixmod9,CCI.mixmod1,SEND.mixmod1,SIRS.mixmod1,infec.mixmod1,delirium.mixmod1,Falls.mixmod1), time_var = "lower", data_Surv = surv_ORCHARD_2admits_train, id_var = "id", 
#                control = list(n_chains = 1, cores = 1, n_iter = 20000, n_burnin=3500), priors = list(penalty_alphas = "ridge"))
# saveRDS(jointFit_int, file = "jointFit_correctData_intcens.rds")

summary(jointFit)

print("done")
