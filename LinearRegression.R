library(car)
library(dplyr)
library(ggplot2)
library(ggpubr)
# options(scipen=999) ##To ensure all values are in decimals (not exponentials)

# Load file
mtb <- read.csv('./ThesisAnalysis_DataPresent_WithoutOutliers_FINAL.csv')

# Age and GDT, CMR ----
GDT_Age <- lm(Age ~ GDTBeh, mtb, na.action = na.exclude)
print(summary(GDT_Age))

CMR_Age <- lm(Age ~ CMRCoh, mtb, na.action = na.exclude)
print(summary(CMR_Age))

GDT_CMR_Age <- lm(Age ~ GDTBeh + CMRCoh, mtb, na.action = na.exclude)
print(summary(GDT_CMR_Age))

Anova(GDT_CMR_Age, type = 2)

AgePredictedAC_WC <- predict(AC_WC_Age)
AgePredictedAC <- predict(AC_Age)
AgePredictedWC <- predict(WC_Age)

# Age and AC, WC Feature ----
AC_Age <- lm(Age ~ ACFeature_Final, mtb, na.action = na.exclude)
print(summary(AC_Age))

WC_Age <- lm(Age ~ WCFeature_Final, mtb, na.action = na.exclude)
print(summary(WC_Age))

AC_WC_Age <- lm(Age ~ ACFeature_Final + WCFeature_Final, mtb, na.action = na.exclude)
print(summary(AC_WC_Age))

Anova(AC_WC_Age, type = 2)

AgePredictedAC_WC <- predict(AC_WC_Age)
AgePredictedAC <- predict(AC_Age)
AgePredictedWC <- predict(WC_Age)

# MST from AC, WC Features ----
AC_MST <- lm(MSTThreshold ~ LapseRate + ACFeature_Final, mtb, na.action = na.exclude)
print(summary(AC_MST))
Anova(AC_MST)

WC_MST <- lm(MSTThreshold ~ LapseRate + WCFeature_Final, mtb, na.action = na.exclude)
print(summary(WC_MST))
Anova(WC_MST)

AC_WC_MST <- lm(MSTThreshold ~ ACFeature_Final + WCFeature_Final + LapseRate, mtb, na.action = na.exclude)
print(summary(AC_WC_MST))

Anova(AC_WC_MST, type = 2)

MSTPredictedAC <- predict(AC_MST)
MSTPredictedWC <- predict(WC_MST)
MSTPredictedAC_WC <- predict(AC_WC_MST)

# MRT from AC, WC Feature ----
AC_MRT <- lm(MRTThreshold ~ LapseRate + ACFeature_Final, mtb, na.action = na.exclude)
print(summary(AC_MRT))

WC_MRT <- lm(MRTThreshold ~ LapseRate + WCFeature_Final, mtb, na.action = na.exclude)
print(summary(WC_MRT))

AC_WC_MRT <- lm(MRTThreshold ~ LapseRate + ACFeature_Final + WCFeature_Final, mtb, na.action = na.exclude)
print(summary(AC_WC_MRT))

MRTPredictedAC <- predict(AC_MRT)
MRTPredictedWC <- predict(WC_MRT)
MRTPredictedAC_WC <- predict(AC_WC_MRT)

# MST from GDT, CMR ----

CMR_MST <- lm(MSTThreshold ~ LapseRate + CMRCoh, mtb, na.action = na.exclude)
print(summary(CMR_MST))
Anova(CMR_MST)

GDT_MST <- lm(MSTThreshold ~ LapseRate + GDTBeh, mtb, na.action = na.exclude)
print(summary(GDT_MST))
Anova(GDT_MST)

CMR_GDT_MST <- lm(MSTThreshold ~ LapseRate + CMRCoh + GDTBeh, mtb, na.action = na.exclude)
print(summary(CMR_GDT_MST))
Anova(CMR_GDT_MST, type = 2)

MSTPredictedCMR <- predict(CMR_MST)
MSTPredictedGDT <- predict(GDT_MST)
MSTPredictedCMR_GDT <- predict(CMR_GDT_MST)

# MRT from GDT, CMR ----

CMR_MRT <- lm(MRTThreshold ~ LapseRate + CMRCoh, mtb, na.action = na.exclude)
print(summary(CMR_MRT))

GDT_MRT <- lm(MRTThreshold ~ LapseRate + GDTBeh, mtb, na.action = na.exclude)
print(summary(GDT_MRT))

CMR_GDT_MRT <- lm(MRTThreshold ~ LapseRate + CMRCoh + GDTBeh, mtb, na.action = na.exclude)
print(summary(CMR_GDT_MRT))

MRTPredictedCMR <- predict(CMR_MRT)
MRTPredictedGDT <- predict(GDT_MRT)
MRTPredictedCMR_GDT <- predict(CMR_GDT_MRT)

# SPIN Scores from Age ----
MST_age <- lm(Age ~ LapseRate + MSTThreshold, mtb, na.action=na.exclude)
summary(MST_age)
Anova(MST_age)

MRT_age <- lm(Age ~ LapseRate + MRTThreshold, mtb, na.action=na.exclude)
summary(MRT_age)

MST_MRT_age <- lm(Age ~ LapseRate + MSTThreshold + MRTThreshold, mtb, na.action=na.exclude)
summary(MST_MRT_age)
Anova(MST_MRT_age)

AgePredictedMST <- predict(MST_age)
AgePredictedMRT <- predict(MRT_age)

# Effect of peripheral hearing on GDT, CMR ----

aud_CMR <- lm(CMRCoh ~ RPTA + LPTA, mtb, na.action = na.exclude)
print(summary(aud_CMR))
Anova(aud_CMR)

aud_GDT <- lm(GDTBeh ~ RPTA + LPTA, mtb, na.action = na.exclude)
print(summary(aud_GDT))

memr_CMR <- lm(CMRCoh ~ HPNMEMR + WBMEMR, mtb, na.action = na.exclude)
print(summary(memr_CMR))

memr_GDT <- lm(GDTBeh ~ HPNMEMR + WBMEMR, mtb, na.action = na.exclude)
print(summary(memr_GDT))

# Effect of peripheral hearing on AC, WC ----

aud_AC <- lm(ACFeature_Final ~ RPTA + LPTA, mtb, na.action = na.exclude)
print(summary(aud_AC))

aud_WC <- lm(WCFeature_Final ~ RPTA + LPTA, mtb, na.action = na.exclude)
print(summary(aud_WC))

memr_AC <- lm(ACFeature_Final ~ HPNMEMR + WBMEMR, mtb, na.action = na.exclude)
print(summary(memr_AC))

memr_WC <- lm(WCFeature_Final ~ HPNMEMR + WBMEMR, mtb, na.action = na.exclude)
print(summary(memr_WC))

# Effect of peripheral hearing on age ----

aud_age <- lm(Age ~ RPTA + LPTA, mtb, na.action = na.exclude)
print(summary(aud_age))

memr_age <- lm(Age ~ HPNMEMR + WBMEMR, mtb, na.action = na.exclude)
print(summary(memr_age))

# PLOTTING ----

dat = data.frame(Predicted = AgePredictedAC_WC,
                 Observed = mtb$Age, 
                 Age = mtb$Age)

dat_age <- dat %>% mutate(age_group = case_when(Age <= 40 ~ "YNH (<=40 y)",
                                                      Age > 40 ~ "ONH (>40 y)",
                                                      TRUE ~ NA_character_))
# Order of the age groups
Ages_order <- c("YNH (<=40 y)", "ONH (>40 y)")
dat_age$age_group <- factor(dat_age$age_group, levels = Ages_order)

#plot predicted vs actual SPIN  
p <- ggscatter(dat_age, x="Predicted", y="Observed", 
          add="reg.line",
          add.params = list(color ="black"),
          color="age_group",palette = "Dark2",
          shape="age_group", use="complete.obs",
          show.legend.text = FALSE) + 
  stat_cor(method="pearson",label.sep = "\n",na.rm=TRUE)+
  scale_x_continuous(limits = c(10, 80), breaks = seq(20, 80, by = 20)) +
  scale_y_continuous(limits = c(10, 80), breaks = seq(20, 80, by = 20)) +
  labs(x ="Predicted Age (years)", y="Actual Age (years)",
       title = "Age ~ AC Feature + WC Feature") +
  theme_bw() + 
  theme(
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

p <- ggpar(p, legend = "bottom",
           legend.title = "Age group")
print(p)

ggsave("./Figures/LinearRegression/AgePredictedbyAC_WC.png", plot=p, dpi = 400, width = 4, height = 4)

# Getting combined plots for paper ----

library(magick)
library(grid)
library(gridExtra)

img1 <- image_read("./Figures/LinearRegression/AgePredictedbyACandWCFeature_6.10.png")
img2 <- image_read("./Figures/LinearRegression/AgePredictedbyWCFeature_6.10.png")
img3 <- image_read("./Figures/LinearRegression/AgePredictedbyACFeature_6.10.png")
# img4 <- image_read("./Figures/LinearRegression/MSTPredictedbyCMRCoherent_6.10.png")

# Convert to raster grobs
g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
g2 <- rasterGrob(as.raster(img2), interpolate = TRUE)
g3 <- rasterGrob(as.raster(img3), interpolate = TRUE)
# g4 <- rasterGrob(as.raster(img4), interpolate = TRUE)

# Define layout: g1 spans the full top row; g2 and g3 go on the second row
layout <- rbind(c(1, 1),
                c(2, 3))

combined <- arrangeGrob(g1, g2, g3, layout_matrix = layout)

combined <- arrangeGrob(g1, g2, g3, g4, ncol=1)

# Save to file
png("./Figures/LinearRegression/MRTByBehavioral_Combined.png", width = 3, height = 6, units = "in", res = 500)
grid.draw(combined)
dev.off()

## Plots with spacers ----

library(grid)
library(gridExtra)
library(magick)

# Load images and convert to raster grobs
img1 <- image_read("./Figures/LinearRegression/AgePredictedbyACandWCFeature_6.10.png")
img2 <- image_read("./Figures/LinearRegression/AgePredictedbyWCFeature_6.10.png")
img3 <- image_read("./Figures/LinearRegression/AgePredictedbyACFeature_6.10.png")
img4 <- image_read("./Figures/LinearRegression/MSTPredictedbyACFeature_6.10.png")

g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
g2 <- rasterGrob(as.raster(img2), interpolate = TRUE)
g3 <- rasterGrob(as.raster(img3), interpolate = TRUE)
g4 <- rasterGrob(as.raster(img4), interpolate = TRUE)

# Create visible separator (e.g., light grey bar or black line)
separator <- rectGrob(gp = gpar(fill = "grey80", col = NA), height = unit(2, "mm"))

# Create two rows of plots
top_row <- arrangeGrob(g1,ncol = 1)
bottom_row <- arrangeGrob(g2, g3, ncol = 2)

# Combine with separator only between rows
combined <- arrangeGrob(
  top_row,
  separator,
  bottom_row,
  ncol = 1,
  heights = unit.c(unit(1, "null"), unit(2, "mm"), unit(1, "null"))
)

# Save to file
png("./Figures/LinearRegression/AgeByEEG_Combined.png", width = 4, height = 4, units = "in", res = 500)
grid.draw(combined)
dev.off()
