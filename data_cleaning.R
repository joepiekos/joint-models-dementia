library(tidyverse)
library(survival)
library(lubridate)
library(GLMMadaptive)
library(JMbayes2)
library(mice)
library(lme4)

#run whole document to create relevant survival analysis variables and do preprocessing

ORCHARD_2017_2019_original_aug <- read_csv("ORCHARD_2017-2019_aug2023.csv")
ORCHARD_2017_2019 <- ORCHARD_2017_2019_original_aug #working copy of dataframe

#remove patients with dementia at baseline
ORCHARD_baseline_dementia_removed2 <- ORCHARD_2017_2019 %>%
    group_by(pn_simple) %>%
    mutate(has_diagnosis = ifelse(row_number() == 1 & (cog_dementiahx == "Yes" | coded_dementia == "Yes"), TRUE, FALSE )) %>% #grepl(codes_to_check_regex, codes)), TRUE, FALSE)
    ungroup() %>%
    filter(!(pn_simple %in% pn_simple[has_diagnosis])) %>%
    dplyr::select(-has_diagnosis) #573 extra people removed, people who had cognitive screen as dementia-yes but didn't have ICD-dementia.

#parse dates  
ORCHARD_baseline_dementia_removed2$dt_admission <- ymd_hms(ORCHARD_baseline_dementia_removed2$dt_admission)
ORCHARD_baseline_dementia_removed2$dt_discharge <- ymd_hms(ORCHARD_baseline_dementia_removed2$dt_discharge)
  
#convert dates into days since first admission
ORCHARD_baseline_dementia_removed2 <- ORCHARD_baseline_dementia_removed2 %>% filter(...1 != 4062 & ...1 != 17155) %>% #removing specific rows that are erroneous, go back to original data and check if you want.
                                      group_by(pn_simple) %>% arrange(dt_admission) %>%
                                      mutate(start = as.numeric(difftime(dt_admission,min(dt_admission),units="days")),end = as.numeric(difftime(dt_discharge,min(dt_admission),units="days"))) 
  
ORCHARD_baseline_dementia_removed2 <- ORCHARD_baseline_dementia_removed2 %>% relocate(start,end,.after = pn_simple)
ORCHARD_baseline_dementia_removed2 <- ORCHARD_baseline_dementia_removed2 %>% group_by(pn_simple) %>%
                                      mutate(start_match = lag(end, default = 0)) %>% mutate(end_match = ifelse(is.na(lead(start)), end, lead(start))) #join disjoint admission intervals
ORCHARD_baseline_dementia_removed2 <- ORCHARD_baseline_dementia_removed2 %>% relocate(start_match,.before = end) %>% relocate(end_match,.after=start)

ORCHARD_baseline_dementia_removed2$pn_simple <- as.factor(ORCHARD_baseline_dementia_removed2$pn_simple) #convert id variable into factor
  
#create single dementia outcomes by combining cognitive screen and ICD-10 code for dementia.
ORCHARD_baseline_dementia_removed2 <- ORCHARD_baseline_dementia_removed2 %>% mutate(has_dementia = ifelse(((cog_dementiahx == "Yes") | (coded_dementia == "Yes") | (cog_dementiahx == "Uncertain" & coded_dementia == "Yes")),1,0))
#if no cognitive screen or ICD available then assume no dementia present:
ORCHARD_baseline_dementia_removed2["has_dementia"][is.na(ORCHARD_baseline_dementia_removed2["has_dementia"])] <- 0
ORCHARD_baseline_dementia_removed2 <- ORCHARD_baseline_dementia_removed2 %>% relocate(has_dementia, .after = end)

#same thing for delirium, combining ICD and cognitive screen  
ORCHARD_baseline_dementia_removed2 <- ORCHARD_baseline_dementia_removed2 %>% mutate(has_delirium = ifelse(((cog_hasdelirium == "Yes") | (coded_delirium == "Yes") | (cog_hasdelirium == "Uncertain" & coded_delirium == "Yes")),1,0))
ORCHARD_baseline_dementia_removed2["has_delirium"][is.na(ORCHARD_baseline_dementia_removed2["has_delirium"])] <- 0
ORCHARD_baseline_dementia_removed2 <- ORCHARD_baseline_dementia_removed2 %>% relocate(has_delirium, .after = cog_hasdelirium) %>% relocate (coded_delirium, .before = has_delirium)

#next block code is to remove admissions that happen after first diagnosis of dementia  
grouped_df <- ORCHARD_baseline_dementia_removed2 %>% arrange(pn_simple, start) %>% group_by(pn_simple) 
lose_diagnosis <- grouped_df %>% filter(any(has_dementia == 1) && any(has_dementia == 0)) %>% filter(min(which(has_dementia == 1)) < max(which(has_dementia == 0))) %>% relocate(coded_dementia,.before=start)
first_occurence <- grouped_df %>% filter(coded_dementia == "Yes") %>% slice(1)
filtered_df <- ORCHARD_baseline_dementia_removed2 %>% left_join(first_occurence, by = "pn_simple") %>% filter(is.na(start.y) | start.x <= start.y)
ORCHARD_analysis_copy <- filtered_df %>% dplyr::select(matches(".x$") | matches("pn_simple")) %>% relocate(coded_dementia.x, .before = has_dementia.x) %>% relocate(pn_simple, .after = pn_long.x) %>% relocate(cog_dementiahx.x, .after = coded_dementia.x)
  

ORCHARD_analysis_copy$Sex.x <- as.factor(ORCHARD_analysis_copy$Sex.x)

#create status variable which is 1 if patient has any diagnosis of dementia, used for implementation of survival model.
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>% group_by(pn_simple) %>% mutate(status = ifelse(any(has_dementia.x == 1),1,0)) %>% relocate(status, .after = has_dementia.x)
#survtime is the censoring/event time if we are doing a non-interval censored model
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>% group_by(pn_simple) %>% mutate(survtime = ifelse(any(has_dementia.x == 1),max(start.x),max(end_match.x))) %>% relocate(survtime, .after = end_match.x)

#second removal of admissions that happen after first diagnosis, to ensure everything is gone
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>%
                          group_by(pn_simple) %>%
                          arrange(pn_simple, start.x) %>%  # Replace "AdmissionDate" with the actual date column
                          mutate(HasDementiaOccurred = cumsum(has_dementia.x == 1),
                          FirstDementiaOccurrence = ifelse(HasDementiaOccurred == 1, 1, 0)) %>%
                          filter(HasDementiaOccurred == 0 | FirstDementiaOccurrence == 1) %>%
                          ungroup() %>%
                          dplyr::select(-HasDementiaOccurred, -FirstDementiaOccurrence)
  
  
# Calculate the percentage of missing data in each column
missing_percentages <- sapply(ORCHARD_analysis_copy, function(x) {
                                  missing_percentage <- sum(is.na(x)) / length(x) * 100
                                  return(missing_percentage)
                             })
  
# Create a data frame to display the missing data results
missing_data_summary <- data.frame(Column = names(ORCHARD_analysis_copy), MissingPercentage = missing_percentages)

#convert falls history and infection into numerical indicators from text.
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>% mutate(Fall_history.x = ifelse(Fall_history.x=="Yes",1,0))
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>% mutate(infection_d1.x = ifelse(infection_d1.x == "No",0,1))
  
#separate dataframe containing the variables to be included in the multiple imputation model
ORCHARD_modelling_copy <- ORCHARD_analysis_copy %>% ungroup() %>% dplyr::select(start.x,end_match.x,survtime,Age.x,Fall_history.x,NEWS1_score.x,SEND_TT_score.x,
                              has_dementia.x,must_height.x,CCI.x,SIRS.x,dem_IMD.x,HFRS.x,infection_d1.x,Sex.x)
  
vars <- c("start.x","end_match.x","survtime","Age.x","Fall_history.x","NEWS1_score.x","SEND_TT_score.x",
            "has_dementia.x","must_height.x","CCI.x","SIRS.x","dem_IMD.x","HFRS.x","infection_d1.x","Sex.x")

#function to impute variables and add the imputed columns back into the main data:
mult_imp_ORCHARD <- function(dataset,vars) {
    imp <- mice(ORCHARD_modelling_copy,m=1,maxit=2)
    c3 <- complete(imp, 1)
    for (variable in vars) {
      print(variable)
      print(sum(is.na(ORCHARD_analysis_copy[[variable]])))
      ORCHARD_analysis_copy[[variable]] <- c3[[variable]]
      print(sum(is.na(ORCHARD_analysis_copy[[variable]])))
    }
    imp_output <- list("imp"=imp, "ORCHARD_analysis_copy"=ORCHARD_analysis_copy)
    return(imp_output)
  }
  
imp_output <- mult_imp_ORCHARD(ORCHARD_modelling_copy,vars)
imp <- imp_output$imp
ORCHARD_analysis_copy <- imp_output$ORCHARD_analysis_copy
imp$loggedEvents #check for imputation errors

#check variables have been imputed, should show zero missing data now for imputed variables
missing_percentages_postimp <- sapply(ORCHARD_analysis_copy, function(x) {
    missing_percentage <- sum(is.na(x)) / length(x) * 100
    return(missing_percentage)
  })
  
# Create a dataframe to display the results
missing_data_summary_postimp <- data.frame(Column = names(ORCHARD_analysis_copy), MissingPercentage = missing_percentages_postimp)

#scale time variables from days to years for model convergence
ORCHARD_analysis_copy$start.x_yr <- ORCHARD_analysis_copy$start.x/365
ORCHARD_analysis_copy$end_match.x_yr <- ORCHARD_analysis_copy$end_match.x/365
ORCHARD_analysis_copy$end.x_yr <- ORCHARD_analysis_copy$end.x/365
ORCHARD_analysis_copy$survtime <- ORCHARD_analysis_copy$survtime/365
  
#create further status variables, for interval censored and right censored models.
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>% mutate(pn_simple = as.numeric(pn_simple)) %>% arrange(pn_simple)
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>% group_by(pn_simple) %>% mutate(status_interval = ifelse(any(has_dementia.x == 1), 3, 0)) %>% relocate(status_interval,.after=status)
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>% group_by(pn_simple) %>% mutate(status_rightcens = ifelse(any(has_dementia.x == 1), 1, 0)) %>% relocate(status_rightcens,.after=status_interval)

#lower and upper are the matched admission times from earlier, but renamed to make more relevant to interval censoring. lower will be used as observation time in GLMMs
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>% mutate(lower = start.x_yr,upper = ifelse(has_dementia.x == 1,end.x_yr,end_match.x_yr)) %>% relocate(lower,.after=end.x) %>% relocate(upper,.after=lower)

#Lower and Upper (capitalised) are the equivalent of survtime for interval censored models i.e. the times between which we know the patient developed dementia.
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>% group_by(pn_simple) %>% mutate(Lower = ifelse(any(status_interval==0),survtime,max(lower[lower != max(lower)])), Upper =  ifelse(any(status_interval==0),Inf,max(lower)))
ORCHARD_analysis_copy <- ORCHARD_analysis_copy %>% relocate(c(Lower,Upper), .after = survtime) %>% relocate(c(start_match.x,end.x,coded_dementia.x,cog_dementiahx.x,status), .after = status_rightcens) %>%
                                                     select(!c(...1.x,...2.x,...3.x,pn_long.x))
  
#write processed data to CSV
write_csv(ORCHARD_analysis_copy, file = "ORCHARD_analysis_copy.csv")

