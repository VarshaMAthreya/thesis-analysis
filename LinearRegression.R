library(car)
library(dplyr)
library(ggplot2)
library("ggpubr")
library(ggpmisc)
options(scipen=999) ##To ensure all values are in decimals (not exponentials)

# Load your CSV file
# mtb <- read.csv('MTB_Analysis/Thesis_Analysis_Latest_11.7.csv')

mtb <- read.csv('FinalPaper/ThesisAnalysis_All_4.3.25.csv')

mtb_age <- mtb %>% mutate(age_group = case_when(Age <= 40 ~ "YNH (<=40 y)",
                                                Age > 40 ~ "ONH (>40 y)",
                                                TRUE ~ NA_character_))

# Predicting Age from AC, WC Features ----

Age_AC <- lm(Age ~ ACFeature_Final, mtb, na.action = na.exclude)
print(summary(Age_AC))

Age_WC <- lm(Age ~ WCFeature_Final, mtb, na.action = na.exclude)
print(summary(Age_WC))

AC_WC_Age <- lm(Age ~ ACFeature_Final + WCFeature_Final, mtb, na.action = na.exclude)
print(summary(AC_WC_Age))

AC_PredictedAge <- predict(Age_AC)
WC_PredictedAge <- predict(Age_WC)

# Predicting SPIN scores from AC, WC Features ----

## MST 
AC_MST <- lm(MSTThreshold ~ ACFeature_Final + LapseRate, mtb, na.action = na.exclude)
print(summary(AC_MST))

WC_MST <- lm(MSTThreshold ~ WCFeature_Final + LapseRate, mtb, na.action = na.exclude)
print(summary(WC_MST))

AC_WC_MST <- lm(MSTThreshold ~ LapseRate + ACFeature_Final + WCFeature_Final, mtb, na.action = na.exclude)
print(summary(AC_WC_MST))

Anova(AC_WC_MST)

AC_PredictedMST <- predict(AC_MST)
WC_PredictedMST <- predict(WC_MST)
AC_WC_PredictedMST <- predict(AC_WC_MST)

## MRT
AC_MRT <- lm(MRTThreshold ~ ACFeature_Final +  LapseRate, mtb, na.action = na.exclude)
print(summary(AC_MRT))

WC_MRT <- lm(MRTThreshold ~ WCFeature_Final + LapseRate, mtb, na.action = na.exclude)
print(summary(WC_MRT))

AC_WC_MRT <- lm(MRTThreshold ~ LapseRate + ACFeature_Final + WCFeature_Final, mtb, na.action = na.exclude)
print(summary(AC_WC_MRT))

Anova(AC_WC_MRT)

AC_PredictedMRT <- predict(AC_MRT)
WC_PredictedMRT <- predict(WC_MRT)
AC_WC_MRT <- predict(AC_WC_MRT)

# Predicting SPIN Scores from Age  ----
MST_age <- lm(MSTThreshold ~ Age + LapseRate, mtb, na.action=na.exclude)
print(summary(MST_age))

MRT_age <- lm(MRTThreshold ~ Age + LapseRate, mtb, na.action=na.exclude)
print(summary(MRT_age))

AgePredictedMST <- predict(MST_age)
AgePredictedMRT <- predict(MRT_age)

# Effect of peripheral hearing on AC, WC (Hearing thresholds and MEMR) ----

aud_AC <- lm(ACFeature_Final ~ RPTA + LPTA, mtb, na.action = na.exclude)
print(summary(aud_AC))

aud_WC <- lm(WCFeature_Final ~ RPTA + LPTA, mtb, na.action = na.exclude)
print(summary(aud_WC))

memr_AC <- lm(ACFeature_Final ~ HPNMEMR + WBMEMR, mtb, na.action = na.exclude)
print(summary(memr_AC))

memr_WC <- lm(WCFeature_Final ~ HPNMEMR + WBMEMR, mtb, na.action = na.exclude)
print(summary(memr_WC))

# Effect of peripheral hearing on age (Hearing thresholds and MEMR) ----

aud_age <- lm(Age ~ RPTA + LPTA, mtb, na.action = na.exclude)
print(summary(aud_age))

memr_age <- lm(Age ~ HPNMEMR + WBMEMR, mtb, na.action = na.exclude)
print(summary(memr_age))

# Plotting (Predicted vs Actual) ----
dat = data.frame(Predicted = AgePredictedMST,
                 Observed = mtb$MSTThreshold, 
                 Age = mtb$Age)

dat_age <- dat %>% mutate(age_group = case_when(Age <= 40 ~ "YNH (<=40 y)",
                                                      Age > 40 ~ "ONH (>40 y)",
                                                      TRUE ~ NA_character_))
# Order of the age groups
Ages_order <- c("YNH (<=40 y)", "ONH (>40 y)")
dat_age$age_group <- factor(dat_age$age_group, levels = Ages_order)

# Plot predicted vs actual SPIN  
p <- ggscatter(dat_age, x="Predicted", y="Observed", 
          add="reg.line",
          add.params = list(color ="black"),
          color="age_group",palette = "Dark2",
          shape="age_group", use="complete.obs",
          show.legend.text = FALSE) + 
  stat_cor(method="pearson",label.sep = "\n",na.rm=TRUE)+
  # scale_x_continuous(limits = c(-5,0), breaks = seq(-5, 0, by = 1)) +
  # scale_y_continuous(limits = c(-5,0), breaks = seq(-5,0, by = 1)) +
  # scale_x_continuous(limits = c(-50,-20), breaks = seq(-50, -20, by = 5)) +
  # scale_y_continuous(limits = c(-50,-20), breaks = seq(-50, -20, by = 5)) +
  labs(x ="Predicted MST Threshold", y="Actual MST Threshold",
       title = "MST ~ Age + Lapse Rate") +
  theme_bw()
p <- ggpar(p, legend = "top",
           legend.title = "Age group")
print(p)

ggsave("./FinalPaper/Figures/LinearRegression/MSTPredictedbyAge.png", plot=p, dpi = 400, width = 5, height = 4, bg = "transparent")