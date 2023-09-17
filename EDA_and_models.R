library(tidyverse)
library(ggsurvfit)
library(survival)
library(gte)
library(icenReg)
library(GLMMadaptive)
library(lattice)
library(latex2exp)
library(Unicode)
library(gtsummary)
library(confintr)
library(gridExtra)
library(grid)

set.seed(123) #to ensure consistent validation and training set throughout development

ORCHARD <- read_csv("ORCHARD_analysis_copy.csv") #read in processed data set
ORCHARD <- dplyr::rename(ORCHARD, c("id"="pn_simple")) #rename id for interpretability

admissions_per_patient <- ORCHARD %>%
  group_by(id) %>%
  summarise(num_admissions = n())
ORCHARD <- admissions_per_patient %>% left_join(ORCHARD, by = "id")
ORCHARD_2admits <- ORCHARD %>% filter(num_admissions > 1)
#ORCHARD_3admits <- ORCHARD %>% filter(num_admissions > 2)


#### Kaplan Meier survival curves ####
ORCHARD_KM <- ORCHARD[!duplicated(ORCHARD$id),]
jpeg("kaplan-meier-curve.jpeg", units = "in", width = 5.5, height = 5, res = 800)
pdf("kaplan-meier-curve.pdf")
#par(cex.axis=1.8,cex.lab=1.5,mfrow=c(1,1),mar = c(4.75,4.5,2,0.7))
survfit2(Surv(survtime, status) ~ 1, data = ORCHARD_KM) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Probability of dementia"
  ) + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22))   
dev.off()

#### Turnbull estimate (NPMLE) for interval censored data ####
ORCHARD <- dplyr::rename(ORCHARD, c("id"="pn_simple", "time"="start.x_yr", "time_end"="end_match.x_yr"))
ORCHARD_has_dementia <- ORCHARD_2admits %>% filter(has_dementia.x == 1) %>% select(start_match.x,start.x) %>% rename("lower"="start_match.x","upper"="start.x")
ORCHARD_censored <- ORCHARD_2admits %>% group_by(id) %>% filter(!any(has_dementia.x == 1)) %>% select(survtime) %>% mutate(upper = Inf) %>%rename("lower"="survtime")
ORCHARD_censored <- ORCHARD_censored[!duplicated(ORCHARD_censored$id),] %>% ungroup() %>% select(lower,upper)
ORCHARD_turnbull <- rbind(ORCHARD_has_dementia,ORCHARD_censored)
ORCHARD_turnbull <- na.omit(ORCHARD_turnbull)
survival_prob <- seq(from=1,to=0, by = -1/100)  # Choose appropriate time points
groupedFit1_np <- ic_np(cbind(lower,upper) ~ 0, data = ORCHARD_turnbull)
groupedFit1_weib <- ic_par(cbind(lower,upper) ~ 0, data = ORCHARD_turnbull, model = "ph", dist = "weibull")
groupedFit1_exp <- ic_par(cbind(lower,upper) ~ 0, data = ORCHARD_turnbull, model = "ph", dist = "exponential")
groupedFit1_lnorm <- ic_par(cbind(lower,upper) ~ 0, data = ORCHARD_turnbull, model = "ph", dist = "lnorm")
groupedFit1_loglogistic <- ic_par(cbind(lower,upper) ~ 0, data = ORCHARD_turnbull, model = "ph", dist = "loglogistic")
groupedFit1_gamma <- ic_par(cbind(lower,upper) ~ 0, data = ORCHARD_turnbull, model = "ph", dist = "gamma")
groupedFit1_generalgamma <- ic_par(cbind(lower,upper) ~ 0, data = ORCHARD_turnbull, model = "ph", dist = "generalgamma")

time_points_np <- getFitEsts(groupedFit1_np, p=survival_prob)
time_points_weib <- getFitEsts(groupedFit1_weib, p=survival_prob)
time_points_exp <- getFitEsts(groupedFit1_exp, p=survival_prob)
time_points_lnorm <- getFitEsts(groupedFit1_lnorm, p=survival_prob)
time_points_loglogistic <- getFitEsts(groupedFit1_loglogistic, p=survival_prob)
time_points_gamma <- getFitEsts(groupedFit1_gamma, p=survival_prob)
time_points_generalgamma <- getFitEsts(groupedFit1_generalgamma, p=survival_prob)

#plots parametric survival curves alongside turnbull estimate.
png("2admits_fulldata_survival_curves.png", width = 7,height = 7, units = "in", res = 500)
lines(groupedFit1_np,lwd=2)
lines(rev(time_points_weib),survival_prob,col="blue",lwd=2)
# lines(rev(time_points_exp),survival_prob,col="red",lwd=2)
lines(rev(time_points_lnorm),survival_prob,col="green",lwd=2)
# lines(rev(time_points_loglogistic),survival_prob,col="orange",lwd=2)
lines(rev(time_points_generalgamma),survival_prob,col="red",lwd=2)
lines(rev(time_points_gamma),survival_prob,col="gold",lwd=2)
legend("topright",  # Adjust the position as needed
       legend = c("Nonparametric", "Weibull", "Log-Normal", "Generalised Gamma", "Gamma"),
       col = c("black", "blue", "green", "red", "gold"),
       lwd = 2,
       title = "Survival curve form"
)
median_survival_time <- getFitEsts(groupedFit1_np, p=0.5)
survival_prob2 <- rev(seq(from=0,to=1, by = 1/2500))
time_points_np2 <- getFitEsts(groupedFit1_np, p=rev(survival_prob2))
turnbull_fit_data <- data.frame(time = time_points_np2, survival_probability = survival_prob2)
line_segment_data <- data.frame(x1 = c(0,median_survival_time),x2 = c(median_survival_time,median_survival_time),y1 = c(0.5,0.5),
                                y2 = c(0.5,0))
def_breaks <- round(labeling::extended(min(turnbull_fit_data$time), max(turnbull_fit_data$time), m = 5))

gg <- ggplot(turnbull_fit_data, aes(x = time, y = survival_probability)) +   geom_hline(yintercept = seq(0, 1, by = 0.125), color = "gray", linetype = "dashed",alpha=0.9) +  # Add horizontal gridlines
  geom_vline(xintercept = seq(0,1000, by = 125), color = "gray", linetype = "dashed",alpha=0.9) +  # Add vertical gridlines
  geom_line(size = 1) +
  geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), 
               data = line_segment_data, linetype = "solid", size = 1, color = "red") +
  scale_x_continuous(breaks = round(c(def_breaks, median_survival_time)),
                     expand = c(0.01, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, NA)) +
  labs(x = "Time (days from first admission)",y = "S(t) = P( T > t )") +
  theme_minimal() +  # Apply a minimal theme
  theme(legend.position = "top",  # Place the legend at the top of the plot
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))  # Remove minor gridlines

png("2admits_turnbull_median.png", width = 7,height = 7, units = "in", res = 500)
print(gg)
dev.off()

#### Longitudinal profiles of data ####

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
) #nice colour palette

random_patient_ids <- ORCHARD_2admits_train %>% filter(num_admissions == 5) %>%
  distinct(id) %>%
  sample_n(15)

ORCHARD_25patients <- ORCHARD_2admits_train %>% filter(id %in% random_patient_ids$id)
ORCHARD_25patients$id <- as.factor(ORCHARD_25patients$id)

png("HFRS_longitudinal profile.png", width = 7,height = 7, units = "in", res = 500)
ggplot(ORCHARD_25patients, aes(x = lower, y = HFRS.x, color = id)) +
  geom_line(size = 1.1) + # Use geom_line to connect data points
  labs(x = "Time", y = "HFRS") +  # Label the axes
  scale_color_manual(values = c25) +  # Apply the custom color palette
  theme_minimal()  +
  guides(color = FALSE) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))
dev.off()


### Data splitting ####
set.seed(123)
ids_2admits <- splitTools::partition(
  ORCHARD_2admits$id,
  p = c(train = 0.8, test = 0.2),
  type = "grouped"
)
ORCHARD_2admits_train <- ORCHARD_2admits[ids_2admits$train,]
ORCHARD_2admits_test <- ORCHARD_2admits[ids_2admits$test,]
write_csv(ORCHARD_2admits_train, file = "ORCHARD_2admits_train.csv")
write_csv(ORCHARD_2admits_test, file = "ORCHARD_2admits_test.csv")

# ids_2admits_3months <- splitTools::partition(
#   ORCHARD_2admits_3months$id,
#   p = c(train = 0.7, test = 0.3),
#   type = "grouped"
# )
# ORCHARD_2admits_3months_train <- ORCHARD_2admits_3months[ids_2admits_3months$train,]
# ORCHARD_2admits_3months_test <- ORCHARD_2admits_3months[ids_2admits_3months$test,]
# write_csv(ORCHARD_2admits_3months_train, file = "ORCHARD_2admits_3months_train.csv")
# write_csv(ORCHARD_2admits_3months_test, file = "ORCHARD_2admits_3months_test.csv")

#### summary table ####

ORCHARD_2admits_train <- read_csv("ORCHARD_2admits_train.csv")
ORCHARD_2admits_test <- read_csv("ORCHARD_2admits_test.csv")

surv_ORCHARD_2admits_train <- ORCHARD_2admits_train[!duplicated(ORCHARD_2admits_train$id),] #reduces dataset down to one row per patient, so is in the form for a survival submodel.
surv_ORCHARD_2admits_test <- ORCHARD_2admits_test[!duplicated(ORCHARD_2admits_test$id),]

surv_ORCHARD_2admits <- ORCHARD_2admits[!duplicated(ORCHARD_2admits$id),]
surv_ORCHARD_2admits_3months <- ORCHARD_2admits_3months_train[!duplicated(ORCHARD_2admits_3months_train$id),]

summarise_ORCHARD <- function(df_to_summarise){
  num_patients <- df_to_summarise %>% nrow(.)
  num_dementia_patients <- df_to_summarise %>% filter(status_interval==3) %>% nrow(.)
  percentage_dementia_patients <- num_dementia_patients/num_patients
  num_women <- table(df_to_summarise$Sex.x) %>% as.numeric(.) %>% .[1]
  percentage_women <- prop.table(table(df_to_summarise$Sex.x)) %>% as.numeric(.) %>% .[1]
  median_demIMD <- median(df_to_summarise$dem_IMD.x)
  IQR_demIMD <- c(quantile(df_to_summarise$dem_IMD.x,1/4),quantile(df_to_summarise$dem_IMD.x,3/4))
  num_demIMD_1to3 <- df_to_summarise %>% filter(between(dem_IMD.x,1,3)) %>% nrow(.)
  percentage_demIMD_1to3 <- num_demIMD_1to3/num_patients
  num_demIMD_8to10 <- df_to_summarise %>% filter(between(dem_IMD.x,8,10)) %>% nrow(.)
  percentage_demIMD_8to10 <- num_demIMD_8to10/num_patients
  median_baseline_age <- median(df_to_summarise$Age.x)
  IQR_Age <- c(quantile(df_to_summarise$Age.x,1/4),quantile(df_to_summarise$Age.x,3/4))
  median_follow_up_time <- df_to_summarise %>% pull(survtime) %>% median(.)
  IQR_follow_up <- c(quantile(df_to_summarise$survtime,1/4),quantile(df_to_summarise$survtime,3/4))
  median_admissions <- median(df_to_summarise$num_admissions)
  IQR_admissions<- c(quantile(df_to_summarise$num_admissions,1/4),quantile(df_to_summarise$num_admissions,3/4))
  num_intCens <- df_to_summarise %>% filter(status_interval == 3) %>% nrow(.)
  num_rCens <- df_to_summarise %>% filter(status_interval == 0) %>% nrow(.)
  summary_df <- data.frame("Number of patients"=num_patients,"Number of dementia cases"=num_dementia_patients,
                           "Percentage with dementia"=percentage_dementia_patients,"Number of women"=num_women,
                           "Percentage women"=percentage_women,"Median deprivation index"=median_demIMD,
                           "IQR deprivation"=paste(IQR_demIMD,collapse = ", "),"High deprivation (1 to 3)"=num_demIMD_1to3,"Percentage high deprivation (1 to 3)"=percentage_demIMD_1to3,"Low deprivation (8 to 10)"=num_demIMD_8to10,"Median survival time of dementia patients"=median_dementia_survival_time,
                           "IQR follow up"=paste(IQR_follow_up, collapse = ", "),"Percentage high deprivation (8 to 10)"=percentage_demIMD_8to10,"Median age"=median_baseline_age,
                           "IQR age"=paste(IQR_Age, collapse = ", "),"Median admissions"=median_admissions,
                           "IQR admissions"=paste(IQR_admissions,collapse= ", "),"Number interval censored"=num_rCens,"Number right censored"=num_intCens)
  return(as.data.frame(t(summary_df)))
  }

summary_ORCHARD_2admits_train <- summarise_ORCHARD(surv_ORCHARD_2admits_train) %>% rename("ORCHARD_train"="V1")
summary_ORCHARD_2admits_test <- summarise_ORCHARD(surv_ORCHARD_2admits_test) %>% rename("ORCHARD_test"="V1")
summary_allORCHARD_2admits <- cbind(summary_ORCHARD_2admits_train,summary_ORCHARD_2admits_test)


tst <- survfit(Surv(survtime, status_rightcens) ~ 1, data = surv_ORCHARD_2admits)
to_remove <- surv_ORCHARD_2admits_3months_train %>% filter(Upper == Inf) %>% slice_sample(prop=0.6) %>% pull(id)
surv_ORCHARD_2admits_turnbull <- surv_ORCHARD_2admits_3months_train %>% filter(!id %in% to_remove) %>% select(Lower,Upper)
surv_ORCHARD_2admits_turnbull$Lower <- ifelse(surv_ORCHARD_2admits_turnbull$Lower == 0,0.0001,surv_ORCHARD_2admits_turnbull$Lower)
tst <- ic_np(cbind(Lower,Upper) ~ 0, data = surv_ORCHARD_2admits_turnbull)
plot(tst)

median_survival_time <- getFitEsts(tst, p=0.5)
survival_prob2 <- rev(seq(from=0,to=1, by = 1/2500))
time_points_np2 <- getFitEsts(tst, p=rev(survival_prob2))
turnbull_fit_data <- data.frame(time = time_points_np2, survival_probability = survival_prob2)
line_segment_data <- data.frame(x1 = c(0,median_survival_time),x2 = c(median_survival_time,median_survival_time),y1 = c(0.5,0.5),
                                y2 = c(0.5,0))
def_breaks <- round(labeling::extended(0, 3, m = 5))

gg <- ggplot(turnbull_fit_data, aes(x = time, y = survival_probability)) +   geom_hline(yintercept = seq(0, 1, by = 0.125), color = "gray", linetype = "dashed",alpha=0.9) +  # Add horizontal gridlines
  geom_vline(xintercept = seq(0,3, by = 0.333), color = "gray", linetype = "dashed",alpha=0.9) +  # Add vertical gridlines
  geom_line(size = 1) +
  geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), 
               data = line_segment_data, linetype = "solid", size = 1, color = "red") +
  scale_x_continuous(breaks = round(c(def_breaks, median_survival_time)),
                     expand = c(0.01, 0), limits = c(0, NA)) +
    scale_y_continuous(expand = c(0.01, 0), limits = c(0, NA)) +
  labs(x = "Time (years from first admission)",y = "S(t) = P( T > t )") +
  theme_minimal() +  # Apply a minimal theme
  theme(legend.position = "top",  # Place the legend at the top of the plot
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))  # Remove minor gridlines

png("2admits_turnbull_median.png",res=400,units="in",width=6.3,height=6.3)
gg
dev.off()
survival_prob <- seq(from=1,to=0, by = -1/100)  # Choose appropriate time points
groupedFit1_weib <- ic_par(cbind(Lower,Upper) ~ 0, data = surv_ORCHARD_2admits_turnbull, model = "ph", dist = "weibull")
groupedFit1_exp <- ic_par(cbind(Lower,Upper) ~ 0, data = surv_ORCHARD_2admits_turnbull, model = "ph", dist = "exponential")
groupedFit1_lnorm <- ic_par(cbind(Lower,Upper) ~ 0, data = surv_ORCHARD_2admits_turnbull, model = "ph", dist = "lnorm")
groupedFit1_loglogistic <- ic_par(cbind(Lower,Upper) ~ 0, data = surv_ORCHARD_2admits_turnbull, model = "ph", dist = "loglogistic")
groupedFit1_gamma <- ic_par(cbind(Lower,Upper) ~ 0, data = surv_ORCHARD_2admits_turnbull, model = "ph", dist = "gamma")
groupedFit1_generalgamma <- ic_par(cbind(Lower,Upper) ~ 0, data = surv_ORCHARD_2admits_turnbull, model = "ph", dist = "generalgamma")

time_points_np <- getFitEsts(tst, p=survival_prob)
time_points_weib <- getFitEsts(groupedFit1_weib, p=survival_prob)
time_points_exp <- getFitEsts(groupedFit1_exp, p=survival_prob)
time_points_lnorm <- getFitEsts(groupedFit1_lnorm, p=survival_prob)
time_points_loglogistic <- getFitEsts(groupedFit1_loglogistic, p=survival_prob)
time_points_gamma <- getFitEsts(groupedFit1_gamma, p=survival_prob)
time_points_generalgamma <- getFitEsts(groupedFit1_generalgamma, p=survival_prob)

png("2admits_fulldata_survival_curves.png", width = 7,height = 7, units = "in", res = 500)
plot(tst,lwd=2,cex.axis = 1.5,cex.lab=1.2,xlab="Time (years)")
lines(rev(time_points_weib),survival_prob,col="blue",lwd=2)
lines(rev(time_points_exp),survival_prob,col="red",lwd=2)
#lines(rev(time_points_lnorm),survival_prob,col="green",lwd=2)
# lines(rev(time_points_loglogistic),survival_prob,col="orange",lwd=2)
lines(rev(time_points_generalgamma),survival_prob,col="green",lwd=2)
#lines(rev(time_points_gamma),survival_prob,col="gold",lwd=2)
legend("bottomleft",  # Adjust the position as needed
       legend = c("Nonparametric (Turnbull)", "Weibull", "Generalised Gamma", "Exponential"),
       col = c("black", "blue", "green", "red"),
       lwd = 2,
       title = "Survival curve form"
)
dev.off()

#### missing data summary ####

#run data_cleaning.R to generate missing percentages
missing_data_table <- missing_data_summary %>% format(scientific = F) %>% mutate(MissingPercentage = as.numeric(MissingPercentage)) %>%
                      filter(Column %in% c("pn_simple","start.x","end_match.x","has_dementia.x","has_delirium.x",
                                           "Age.x","Sex.x","dem_IMD.x","Fall_history.x", "CCI.x",
                                           "SIRS.x","HFRS.x","infection_d1.x","SEND_TT_score.x")) %>% mutate(MissingPercentage = round(MissingPercentage,4))
missing_data_table$MissingPercentage <- round(missing_data_table$MissingPercentage, 4)

missing_data_table$Column <- c("id","Admission start time","Admission end time","Dementia","Delirium",
                               "Age","Sex","Index of Multiple Deprivation (IMD)","Fall history", "Charlson Comorbidity Index (CCI)",
                               "Systemic inflammatory response syndrome (SIRS)","Hospital Frailty Risk Score (HFRS)"," Primary diagnosis is infection","SEND score (vital signs)")
rownames(missing_data_table) <- missing_data_table$Column
missing_data_table$Column <- NULL
colnames(missing_data_table) <- c("% of data missing")
xtable(missing_data_table)

#### plot showing the process of deriving time variables ####

dateExample_ORCHARD <- ORCHARD_2admits_train %>% filter(id %in% c(3236))
date_table1 <- dateExample_ORCHARD %>% select(id,has_dementia.x,dt_admission.x,dt_discharge.x) %>% mutate_if(is.numeric, round, digits = 4) %>% rename("ID"="id","Dementia?"="has_dementia.x",
                                                                                                                                                       "Admission date"="dt_admission.x","Discharge date"="dt_discharge.x")
date_table2 <- dateExample_ORCHARD %>% select(id,has_dementia.x,start.x,end.x) %>% mutate_if(is.numeric, round, digits = 4) %>% rename("ID"="id","Dementia?"="has_dementia.x","Admission date (days)"="start.x","Discharge date (days)"="end.x")
date_table3 <- dateExample_ORCHARD %>% select(id,has_dementia.x,start.x,end_match.x) %>% mutate_if(is.numeric, round, digits = 4) %>% rename("ID"="id","Dementia?"="has_dementia.x","Admission date (days)"="start.x","Matched discharge date (days)"="end_match.x")
date_table4 <- dateExample_ORCHARD %>% select(id,has_dementia.x,start.x_yr,end_match.x_yr) %>% mutate_if(is.numeric, round, digits = 4) %>% rename("ID"="id","Dementia?"="has_dementia.x","Admission date (years)"="start.x_yr"," Matched discharge date (years)"="end_match.x_yr")
date_table5 <- dateExample_ORCHARD %>% select(id,has_dementia.x,start.x_yr,end_match.x_yr,survtime,Lower,Upper,lower) %>% mutate_if(is.numeric, round, digits = 4) %>% rename("ID"="id","Dementia?"="has_dementia.x","Admission date \n(years)"="start.x_yr"," Matched discharge \ndate (years)"="end_match.x_yr",
                                                                                                                                                                              "Event/censoring time \n(years)*"="survtime","Lower interval \n censored time \n(years)"="Lower","Upper interval \n censored time \n(years)"="Upper","Observation time \n(years)"="lower")
date_tablegrob1 <- tableGrob(date_table1)
ggsave('date_tablegrob1.png', date_tablegrob1)
date_tablegrob2 <- tableGrob(date_table2)
ggsave('date_tablegrob2.png', date_tablegrob2, units = c("in"), width = 10,height = 4)
date_tablegrob3 <- tableGrob(date_table3)
ggsave('date_tablegrob3.png', date_tablegrob3)
date_tablegrob4 <- tableGrob(date_table4)
ggsave('date_tablegrob4.png', date_tablegrob4)
date_tablegrob5 <- tableGrob(date_table5)
ggsave('date_tablegrob5.png', date_tablegrob5, units = c("in"), width = 15,height = 4)
grid.arrange(date_tablegrob1,date_tablegrob2,date_tablegrob3,date_tablegrob4,date_tablegrob5, ncol = 1)

#### survtime threshold counts ####

# Create a vector of thresholds within the range [0, 3]
thresholds <- seq(0, 3, by = 0.1)  # Adjust the step size as needed

# Initialize an empty vector to store counts
count_above_threshold_all <- numeric(length(thresholds))
count_above_threshold_train <- numeric(length(thresholds))
count_above_threshold_test <- numeric(length(thresholds))

surv_ORCHARD_2admits <- ORCHARD_2admits %>% #change depending on if using right or interval censoring
  group_by(id) %>% filter(if (any(status_interval == 2)) status_interval == 2 else row_number() == n())

# Loop through each threshold and count how many patients have survtime above it
for (i in 1:length(thresholds)) {
  count_above_threshold_all[i] <- sum(surv_ORCHARD_2admits$survtime > thresholds[i],na.rm = TRUE)
  count_above_threshold_train[i] <- sum(ORCHARD_2admits_train$survtime > thresholds[i],na.rm = TRUE)
  count_above_threshold_test[i] <- sum(ORCHARD_2admits_test$survtime > thresholds[i], na.rm = TRUE)
}

survtime_thresholds <- data.frame(thresholds, count_above_threshold_all)
survtime_thresholds_tbl <- survtime_thresholds %>% filter(thresholds %in% c(0,0.5,1,1,5,2,2.5,3))
survtime_thresholds_tbl[1,2] <- 9155

mytheme <- gridExtra::ttheme_default(core = list(padding = unit(c(9, 3), "mm"), fg_params=list(cex = 0.95)),colhead = list(fg_params=list(cex = 0.9)))

png("survtime_thresholds_alldata.png",units="in",width=8,height=5.5,res=400)
tbl <- tableGrob(survtime_thresholds_tbl, cols = c("Threshold","Number of patients"),theme = mytheme, rows = NULL)
# tbl <- gtable::gtable_add_grob(tbl,
#                              grobs = rectGrob(gp=gpar(fill=NA,
#                                                       lwd=2)), t = 1,b = nrow(tbl), l = 1, r = ncol(tbl))
# tbl <- gtable::gtable_add_grob(tbl,
#                                grobs = rectGrob(gp=gpar(fill=NA,
#                                                         lwd=2)), t = 1,b = nrow(tbl), l = 2, r = ncol(tbl))
# tbl <- gtable::gtable_add_grob(tbl,
#                                grobs = rectGrob(gp=gpar(fill=NA,
#                                                         lwd=2)), t = 2,b = nrow(tbl), l = 1, r = ncol(tbl))
ggplot(survtime_thresholds, aes(x = thresholds, y = count_above_threshold_all)) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "red", size = 2.2, shape = 19) +
  labs(
    x = "Threshold (years)",
    y = "Number of Patients",
    title = "Total Patients with Event/Censoring Time Above Thresholds"
  ) +
  scale_y_continuous(breaks = seq(0, max(survtime_thresholds$count_above_threshold_all), by = 2000)) +
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 15),plot.title = element_text(size = 13.5)) +
  annotation_custom(tbl,xmin = 2,xmax = 3,ymin=7000,ymax=8000)
dev.off()

### joint model summaries ####

jointFit_correctData_intcens <- readRDS("jointFit_correctData_intcens.rds")
jointFit_correctData_rcens <- readRDS("jointFit_correctData_rcens.rds") #for model estimates

summary(jointFit_correctData_intcens)
summary(jointFit_correctData_rcens)

jointFit_rcens_3months_infecDeliriumCCI <- readRDS("jointFit_rcens_3months_infecDeliriumCCI.rds")
jointFit_rcens_3months_infecDelirium <- readRDS("jointFit_rcens_3months_infecDelirium.rds")
jointFit_rcens_3months_10kIter <- readRDS("jointFit_rcens_3months_10kIter.rds") # good traceplots

png("deliriumHFRS_MCMCtraceplot.png", units = "in", width = 10, height = 5, res = 400)
ggtraceplot.jm1(jointFit_rcens_3months_10kIter, "alphas", grid = TRUE, gridcols = 2,gridrows = 2,size = 0.6)
dev.off()

png("other_MCMCtraceplot.png", units = "in", width = 10, height = 9.8, res = 400)
ggtraceplot.jm1(jointFit_correctData_intcens, "alphas", grid = TRUE, gridcols = 2,gridrows = 4,size = 0.6)
dev.off()

png("survivalVars_MCMCtraceplot.png", units = "in", width = 10, height = 5, res = 400)
ggtraceplot.jm1(jointFit_correctData_intcens, "gammas", grid = TRUE, gridcols = 2,gridrows = 2,size = 0.6)
dev.off()
################################ GLMMs for each longitudinal variable ###############################################################
#### HFRS ####
ORCHARD_2admits_train$HFRS.x <- round(ORCHARD_2admits_train$HFRS.x)
ORCHARD_2admits_test$HFRS.x <- round(ORCHARD_2admits_test$HFRS.x)
ORCHARD_3admits_train$HFRS.x <- round(ORCHARD_3admits_train$HFRS.x)
ORCHARD_3admits_test$HFRS.x <- round(ORCHARD_3admits_test$HFRS.x)

# see bottom of page for previous HFRS codes

HFRS.mixmod9 <- readRDS("HFRS_hurdle_negbin_2admits_timerand_timezifix_1zirand.rds")
HFRS.mixmod9.DHARMa <- resids_plot(HFRS.mixmod9, ORCHARD_2admits_train$HFRS.x)

# HFRS.mixmod8.pen <- readRDS("HFRS_hurdle_negbin_pen_timerand_timezifix_1zirand.rds")
# HFRS.mixmod8.pen.DHARMa <- resids_plot(HFRS.mixmod8.pen, ORCHARD_3admits_train$HFRS.x)

#### CCI ####
summary(ORCHARD_2admits_train$CCI.x)

CCI.mixmod1 <- readRDS("CCI_mixmod1.rds")
CCI.mixmod1.DHARMa <- resids_plot(HFRS.mixmod8.pen, ORCHARD_3admits_train$HFRS.x)

#### infec ####

infec.mixmod1 <- readRDS("infec_mixmod1.rds")
preds_infec1 <- predict(infec.mixmod1, newdata = ORCHARD_2admits_test[ORCHARD_2admits_test$id %in% c(341,7107), ], 
                       type = "subject_specific",se.fit = TRUE, return_newdata = TRUE)

#### delirium ####

delirium.mixmod1 <- readRDS("delirium_mixmod1.rds")
delirium.mixmod1.DHARMa <- resids_plot(delirium.mixmod1, ORCHARD_2admits_train$has_delirium.x)

#### SEND ####

SEND.mixmod1 <- readRDS("SEND_mixmod1.rds")
SEND.mixmod1.DHARMa <- resids_plot(SEND.mixmod1, ORCHARD_2admits_train$SEND_TT_score.x)

#### SIRS ####
SIRS.mixmod1 <- readRDS("SIRS_mixmod1.rds")
SIRS.mixmod1.DHARMa <- resids_plot(SIRS.mixmod1, ORCHARD_2admits_train$SIRS.x)

#### Falls model ####
Falls.mixmod1 <- readRDS("Falls_mixmod1.rds")
Falls.mixmod1.DHARMa <- resids_plot(Falls.mixmod1, ORCHARD_2admits_train$Fall_history.x)

### DHARMa residuals ###

ORCHARD_2admits_3months_train <- read_csv("ORCHARD_2admits_3months_train.csv")
ORCHARD_2admits_train <- read_csv("ORCHARD_2admits_train.csv")


GLMM_delirium_3months <- readRDS("delirium_3months.rds")
GLMM_HFRS_3months <- readRDS("HFRS_3months.rds")
GLMM_infec_3months <- readRDS("infec_3months.rds")
GLMM_CCI <- readRDS("CCI.rds")
GLMM_Falls <- readRDS("Falls.rds")
GLMM_SIRS <- readRDS("SIRS.rds")

png("DHARMa_residual_plots_HFRS.png", units = "in", width = 8,height = 4,res = 400)
resids_plot(GLMM_HFRS_3months, ORCHARD_2admits_3months_train$HFRS.x)
dev.off()
png("DHARMa_residual_plots_delirium.png", units = "in", width = 8,height = 4,res = 400)
resids_plot(GLMM_delirium_3months, ORCHARD_2admits_3months_train$has_delirium.x)
dev.off()
png("DHARMa_residual_plots_infection.png", units = "in", width = 8,height = 4,res = 400)
resids_plot(GLMM_infec_3months, ORCHARD_2admits_3months_train$infection_d1.x)  
dev.off()
png("DHARMa_residual_plots_CCI.png", units = "in", width = 8,height = 4,res = 400)
resids_plot(GLMM_CCI, ORCHARD_2admits_train$CCI.x)
dev.off()
png("DHARMa_residual_plots_Falls.png", units = "in", width = 8,height = 4,res = 400)
resids_plot(GLMM_Falls, ORCHARD_2admits_train$Fall_history.x)
dev.off()
png("DHARMa_residual_plots_SIRS.png", units = "in", width = 8,height = 4,res = 400)
resids_plot(GLMM_SIRS, ORCHARD_2admits_train$SIRS.x)
dev.off()
