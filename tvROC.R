#tvROC accuracy measure

library(JMbayes2)
library(tidyverse)

print("reading ORCHARD")

ORCHARD_2admits_test <- read_csv("ORCHARD_2admits_test.csv")

print("read ORCHARD, reading model")

jmFit <- readRDS("jointFit_3chains_rcens_penalty.rds")

print("model loaded")
print("doing tvROC()")

#roc <- tvROC(jmFit, newdata = ORCHARD_2admits_test, Tstart = 1, Dt = 1.95)

#saveRDS(roc, file = "roc_rcens_penalty.rds")

print("skipping first tvROC")
print("next tvROC")

roc2 <- tvROC(jmfit, newdata = ORCHARD_2admits_test, Tstart = 1, Dt = 1.95, eventData_fun = function(x) x %>% group_by(id) %>% filter(if (any(status_interval==2)) status_interval == 2 else row_number() == n()))

saveRDS(roc2, file = "roc_eventdataFun_rcens_penalty2.rds2")

