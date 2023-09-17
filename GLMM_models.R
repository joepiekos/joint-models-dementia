library(tidyverse)
library(ggsurvfit)
library(survival)
library(gte)
library(icenReg)
library(GLMMadaptive)
library(lattice)
set.seed(123)

ORCHARD <- read_csv("ORCHARD_analysis_copy.csv")
ORCHARD <- dplyr::rename(ORCHARD, c("id"="pn_simple", "time"="start.x_yr", "time_end"="end_match.x_yr"))

admissions_per_patient <- ORCHARD %>%
  group_by(id) %>%
  summarise(num_admissions = n())
ORCHARD <- admissions_per_patient %>% left_join(ORCHARD, by = "id")
ORCHARD_2admits <- ORCHARD %>% filter(num_admissions > 1)
ORCHARD_3admits <- ORCHARD %>% filter(num_admissions > 2)

set.seed(123)
ids_2admits <- splitTools::partition(
  ORCHARD_2admits$id,
  p = c(train = 0.8, test = 0.2),
  type = "grouped"
)
ORCHARD_2admits_train <- ORCHARD_2admits[ids_2admits$train,]
ORCHARD_2admits_test <- ORCHARD_2admits[ids_2admits$test,]

ids_3admits <- splitTools::partition(
  ORCHARD_3admits$id,
  p = c(train = 0.8, test = 0.2),
  type = "grouped"
)
ORCHARD_3admits_train <- ORCHARD_3admits[ids_3admits$train,]
ORCHARD_3admits_test <- ORCHARD_3admits[ids_3admits$test,]
#### HFRS ####
ORCHARD_2admits_train$HFRS.x <- round(ORCHARD_2admits_train$HFRS.x)
ORCHARD_2admits_test$HFRS.x <- round(ORCHARD_2admits_test$HFRS.x)
ORCHARD_3admits_train$HFRS.x <- round(ORCHARD_3admits_train$HFRS.x)
ORCHARD_3admits_test$HFRS.x <- round(ORCHARD_3admits_test$HFRS.x) #rounding HFRS because I use Poisson model to model it.

print("HFRS models:")
#HFRS.mixmod5 <- mixed_model(HFRS.x ~ time, data = ORCHARD_3admits_train,random = ~ 1 + time| id, zi_fixed = ~1, zi_random = ~1|id,family = hurdle.negative.binomial(), verbose = TRUE)
#saveRDS(HFRS.mixmod5,file = "HFRS_hurdle_negbin_timerand_1zifix_1zirand.rds")

#HFRS.mixmod6 <- mixed_model(HFRS.x ~ time, data = ORCHARD_3admits_train,random = ~ 1 + time| id, zi_fixed = ~1, zi_random = ~1|id,family = zi.negative.binomial(), verbose = TRUE)
#saveRDS(HFRS.mixmod6,file = "HFRS_zi_negbin_timerand_1zifix_1zirand.rds")

#HFRS.mixmod7 <- mixed_model(HFRS.x ~ time, data = ORCHARD_3admits_train,random = ~ 1 + time| id, zi_fixed = ~time, zi_random = ~time|id,family = zi.poisson(), verbose = TRUE)
#saveRDS(HFRS.mixmod7,file = "HFRS_zi_pois_timerand_timezifix_timezirand.rds")

# HFRS.mixmod8 <- mixed_model(HFRS.x ~ time, data = ORCHARD_3admits_train,random = ~ 1 + time| id, zi_fixed = ~time, zi_random = ~1|id,family = hurdle.negative.binomial(), verbose = TRUE)
# saveRDS(HFRS.mixmod8,file = "HFRS_hurdle_negbin_timerand_timezifix_1zirand.rds")

#HFRS.mixmod9 <- mixed_model(HFRS.x ~ time, data = ORCHARD_2admits_train,random = ~ 1 + time| id, zi_fixed = ~time, zi_random = ~1|id,family = hurdle.negative.binomial(), verbose = TRUE)
#saveRDS(HFRS.mixmod9,file = "HFRS_mixmod9")
# 
# HFRS.mixmod8.pen <- mixed_model(HFRS.x ~ time, data = ORCHARD_3admits_train,random = ~ 1 + time| id, zi_fixed = ~time, zi_random = ~1|id,family = hurdle.negative.binomial(), 
#                                 penalized = TRUE,verbose = TRUE)
# saveRDS(HFRS.mixmod8.pen,file = "HFRS_hurdle_negbin_pen_timerand_timezifix_1zirand.rds")

# HFRS.mixmod5.2admits <- mixed_model(HFRS.x ~ time, data = ORCHARD_2admits_train,random = ~ 1 + time| id, zi_fixed = ~1, zi_random = ~1|id,family = hurdle.negative.binomial(), verbose = TRUE)
# saveRDS(HFRS.mixmod5.2admits,file = "HFRS_hurdle_negbin_2admits_timerand_1zifix_1zirand.rds")
print("done")

#### CCI ####
print("CCI models:")
CCI.mixmod1 <- mixed_model(CCI.x ~ time, data = ORCHARD_2admits_train,random = ~ 1 + time|id, zi_fixed = ~time, zi_random = ~1|id,family = hurdle.negative.binomial(), verbose = TRUE)
saveRDS(CCI.mixmod1,file = "CCI_mixmod1.rds")
print("done")

#### infection_d1 ####
print("infection models:")
infec.mixmod1 <- mixed_model(infection_d1.x ~ time, data = ORCHARD_2admits_train,random = ~ 1|id,family = binomial(), verbose = TRUE, control = list(iter_EM=0))
saveRDS(infec.mixmod1,file = "infec_mixmod1.rds")

# infec.mixmod3 <- mixed_model(infection_d1.x ~ time, data = ORCHARD_2admits_train,random = ~ 1 + time|id, zi_fixed = ~time,family = zi.binomial(), verbose = TRUE, control = list(iter_EM=0))
# saveRDS(infec.mixmod3,file = "infec_mixmod3.rds")
print("done")

#### delirium ####
print("delirium models:")
delirium.mixmod1 <- mixed_model(has_delirium.x ~ time, data = ORCHARD_2admits_train,random = ~ 1|id,family = binomial(), verbose = TRUE)
saveRDS(delirium.mixmod1,file = "delirium_mixmod1.rds")
# 
# delirium.mixmod2 <- mixed_model(has_delirium.x ~ time, data = ORCHARD_3admits_train,random = ~ 1 + time|id,family = binomial(), verbose = TRUE)
# saveRDS(delirium.mixmod2,file = "delirium_mixmod2.rds")
# 
# delirium.mixmod3 <- mixed_model(has_delirium.x ~ time, data = ORCHARD_2admits_train,random = ~ 1 + time|id, zi_fixed = ~time,family = zi.binomial(), verbose = TRUE)
# saveRDS(delirium.mixmod3,file = "delirium_mixmod3.rds")
print("done")

#### SEND SCORE ####
print("SEND models:")
SEND.mixmod1 <- mixed_model(SEND_TT_score.x ~ time, data = ORCHARD_2admits_train,random = ~ 1 + time|id, zi_fixed = ~time, zi_random = ~1|id,family = hurdle.negative.binomial(), verbose = TRUE)
saveRDS(SEND.mixmod1,file = "SEND_mixmod1.rds")
print("done")

#### falls history ####
print("fall history models:")
Falls.mixmod1 <- mixed_model(Fall_history.x, ~ time, data = ORCHARD_2admits_train,random = ~ 1|id,family = binomial(), verbose = TRUE, control = list(iter_EM=0))
saveRDS(Falls.mixmod1, file = "Falls_mixmod1.rds")
print("done")

### SIRS ####
# print("SIRS models:")
# # SIRS.mixmod1 <- mixed_model(SIRS.x, ~ time, data = ORCHARD_2admits_train, random = ~ 1 + time|id, family = GLMMadaptive::negative.binomial(), verbose = TRUE)
# # saveRDS(SIRS.mixmod1,file = "SIRS_mixmod1.rds")
# 
# SIRS.mixmod2 <- mixed_model(SIRS.x, ~ time, data = ORCHARD_2admits_train, random = ~ 1 + time|id, zi_fixed = ~time, zi_random = ~1|id,family = hurdle.negative.binomial(), verbose = TRUE)
# saveRDS(SIRS.mixmod2,file = "SIRS_mixmod2.rds")
print("done")
