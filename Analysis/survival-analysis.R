#### loading required packages ####

library(tidyverse)
# loads tidyverse for a bunch of useful data manipulation
# and plotting packages
library(survival)
# loads the survival package for survival analysis
library(ggsurvfit)
# loads the ggsurvfit package for plotting K-M survival curves
library(gtsummary)
# loads the gtsummary package for creating nice tables for
# regression model output

#### data ####

cancer_data <- read_csv("Data/bc_data.csv")
# loads data related to survival of patients with
# breat cancer

names(cancer_data)
# lists the column names of the cancer_data data
# frame. There are 5 columns, dep is a patient's
# deprivation score, region is the area they live,
# agediag is a patient's age at diagnosis, dead is
# whether the patient died during the study, survtime
# is the length of time a patient was known to survive
# after diagnosis

str(cancer_data)
# shows the data types of each column. dep, region 
# and dead are character variables, it makes more
# sense to recode these variables as factors

cancer_data$dep <- factor(cancer_data$dep, 
                          levels = c("leastdep", 2, 3, 4, "mostdep"))
# converts dep variable to a factor with suitable levels
cancer_data$region <- factor(cancer_data$region)
# converts the region variable to a factor
cancer_data$dead <- factor(cancer_data$dead, 
                           levels = c("alive/ce", "dead"))
# converts dep to a factor with suitable levels

str(cancer_data)
# shows the data types of each column to check the
# variables have been converted to factors correctly

summary(cancer_data)
# provides summary counts and statistics for the whole
# dataset

#### visualising the data ####

ggplot(cancer_data, aes(x = dep)) +
  geom_bar(fill = "turquoise") +
  labs(x = "Deprivation Score", y = "Number of Patients",
       title = "Distribution of deprivation score for breast cancer patients in the data")
# creates a bar plot visualising the number of patients
# of each deprivation score in the data

ggplot(cancer_data, aes(x = dep, fill = region)) +
  geom_bar() +
  labs(x = "Deprivation Score", y = "Number of Patients",
       title = "Distribution of deprivation score for breast cancer patients in each region in the data") +
  facet_wrap(region ~.)
# creates a bar plot visualising the number of patients
# of each deprivation score in each region in the data

ggplot(cancer_data, aes(y = agediag)) +
  geom_boxplot(fill = "indianred2") +
  labs(y = "Age at Diagnosis", 
       title = "Distribution of age at diagnosis for breast cancer patients in the data")
# creates a boxplot to visualise the distribution of 
# agediag for patients in the data

ggplot(cancer_data, aes(x = dep, y = agediag, fill = dep)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, col = "blue") +
  labs(x = "Deprivation index", y = "Age at diagnosis", 
       title = "Distribution of age at diagnosis by deprivation index for breast cancer patients")
# creates boxplots of agediag by dep to visualise distribution
# of agediag in different dep groups for breast cancer patients
# in the data. The mean is also marked as a blue diamond on
# each boxplot

ggplot(cancer_data, aes(x = region, y = agediag, fill = region)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, col = "blue") +
  labs(x = "Region", y = "Age at diagnosis", 
       title = "Distribution of age at diagnosis of breast cancer patients by region")
# plots boxplots of age at diagnosis by region to visualise 
# the distribution of agediag in each region for breast cancer
# patients. Mean is marked as blue diamond on each boxplot

ggplot(cancer_data, aes(x = survtime, fill = dead)) +
  geom_histogram(bins = 10, position = "dodge") +
  labs(x = "Survival time", 
       title = "Histogram of survival time in study for breast cancer patients")
# plots a histogram showing the distribution of survival time
# for breast cancer patients. The blue histogram shows the
# distribution of survival time for patients who's deaths were
# observed during the study and the red histogram is for
# patients who survived or for whose death was censored

ggplot(cancer_data, aes(x = dep, y = survtime, fill = dep)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, col = "blue") +
  labs(x = "Deprivation index", y = "Survival time", 
       title = "Distribution of survival for patients with breast cancer by deprivation index")
# creates boxplots showing how survival time is distributed for 
# patients with each deprivation score. means are marked as blue
# diamonds

ggplot(cancer_data, aes(x = dep, fill = dead)) +
  geom_bar() +
  labs(x = "deprivation index", 
       title = "Number of patients in each deprivation group by area, coloured by survival status") +
  facet_wrap(~region)
# creates a bar plot showing distribution of deprivation by region.
# survival status of patients is also coloured on each bar, blue
# is the proportion of patients whose deaths were observed, red
# is the patients who survived or deaths were censored

#### Kaplan-Meier method ####

cancer_data$dead1 <- as.numeric(cancer_data$dead) - 1
# creates a new variables in the dataset, taking a value of
# 0 for a patient who survived/death is uncensored and 1
# for a patient whose death was observed

Surv(cancer_data$survtime, cancer_data$dead1)[1:10]
# shows the survival times of the first 10 patients in the
# data. "+" indicates a censored observation, ":dead"
# indicates a patient whose death was observed

survfit_cd = survfit(Surv(survtime, dead1) ~ 1, data = cancer_data)
# creates a survival object needed to creates survival
# curves using the Kaplan-Meier method for the breast
# cancer survival times

survfit2(Surv(survtime, dead1) ~ 1, data = cancer_data) %>% 
  ggsurvfit() +
  labs(x = "Years", y = "Survival Probability",
       title = "Kaplan-Meier survival curve for patients with breast cancer data") +
  add_confidence_interval()
# creates a Kaplan-Meier curve to show the estimated survival
# probability of a patient with breast cancer at each point
# in time. Plot also has a 95% confidence interval, to show
# the 95% confidence interval for patient survival probability

#### Cox proportional hazards model ####

patient_id = 1:69199
# creates a list of numbers from 1 to 69199 to use as
# a patient id number

cancer_data$id <- patient_id
# adds the patient id number to cancer_data to distinguish
# each patient

cox_bc <- coxph(Surv(survtime, dead) ~ agediag + dep + region, 
                data = cancer_data, id = id)
# fits a coxph survival model to the data, using agediag, dep and
# region to predict hazard function for patients with breast cancer

tbl_regression(cox_bc)
# creates a nice table to show the statistics for the cox
# proportionate hazards model variables

tbl_regression(cox_bc, exp = TRUE)
# creates the same table as above except the output
# has been exponentiated to give statistics for the
# hazards ratio rather than log hazard ratio
