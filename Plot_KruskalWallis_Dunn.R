library(dplyr)
library(FSA)         # For Dunn's Test
library(rstatix)     # For formatting output (optional)
library(tidyr)
library(ggpubr)      # Optional, for plots
options(scipen = 999) 

# Load your CSV
mtb <- read.csv('./ThesisAnalysis_All_WithoutOutliers_5.30.csv')

# Create age groups
mtb_age <- mtb %>%
  mutate(age_group = case_when(
    Age <= 40 ~ "YNH (<=40 y)",
    Age > 40 ~ "ONH (>40 y)",
    TRUE ~ NA_character_
  )) %>%
  mutate(age_group = factor(age_group, levels = c("YNH (<=40 y)", "ONH (>40 y)")))

# Define dependent variables

dv_vars = c("RPTA", "LPTA")

##Summary statistics

pta_summary <- mtb_age %>%
  group_by(age_group) %>%
  summarise(
    RPTA_mean = mean(RPTA, na.rm = TRUE),
    RPTA_sd = sd(RPTA, na.rm = TRUE),
    RPTA_median = median(RPTA, na.rm = TRUE),
    RPTA_min = min(RPTA, na.rm = TRUE),
    RPTA_max = max(RPTA, na.rm = TRUE),
    
    LPTA_mean = mean(LPTA, na.rm = TRUE),
    LPTA_sd = sd(LPTA, na.rm = TRUE),
    LPTA_median = median(LPTA, na.rm = TRUE),
    LPTA_min = min(LPTA, na.rm = TRUE),
    LPTA_max = max(LPTA, na.rm = TRUE)
  )

# dv_vars <-  c("Incoh12P1Lat4chan", "Incoh12P1Amp4chan",
#               "Incoh12P2Lat4chan", "Incoh12P2Amp4chan", "Incoh12SS_4chan",
#               "Incoh20P1Lat4chan", "Incoh20P1Amp4chan",
#               "Incoh20P2Lat4chan", "Incoh20P2Amp4chan", "Incoh20SS_4chan",
#               "Coh12Lat4chan", "Coh12Amp4chan", "Coh12SS_4chan",
#               "Coh20Lat4chan", "Coh20Amp4chan", "Coh20SS_4chan")
# # "SSDiff12", "SSDiff20")

# dv_vars <-  c("GapP14chanAmp16", "GapP14chanLat16",
#               "GapP24chanAmp16", "GapP24chanLat16",
#               "GapP14chanAmp32", "GapP14chanLat32",
#               "GapP24chanAmp32", "GapP24chanLat32",
#               "GapP14chanAmp64", "GapP14chanLat64",
#               "GapP24chanAmp64", "GapP24chanLat64")

# Filter to complete cases
dat_mtb <- mtb_age %>%
  select(all_of(c(dv_vars, "age_group"))) %>%
  filter(complete.cases(.))

run_kruskal_dunn <- function(data, dv_vars, group_var = "age_group", p_adj_method = "bonferroni") {
  results <- data.frame(
    Variable = character(),
    KW_ChiSq = numeric(),
    KW_p = numeric(),
    KW_signif = character(),
    Dunn_comparison = character(),
    Dunn_Z = numeric(),
    Dunn_p = numeric(),
    Dunn_p_adj = numeric(),
    Dunn_signif = character(),
    stringsAsFactors = FALSE
  )
  
  for (dv in dv_vars) {
    formula <- as.formula(paste(dv, "~", group_var))
    kw_test <- kruskal.test(formula, data = data)
    
    kw_sig <- cut(kw_test$p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                  labels = c("***", "**", "*", "ns"))
    
    if (kw_test$p.value < 0.05) {
      dunn <- dunnTest(formula, data = data, method = p_adj_method)$res
      
      dunn_result <- dunn %>%
        mutate(
          Variable = dv,
          KW_ChiSq = as.numeric(kw_test$statistic),
          KW_p = kw_test$p.value,
          KW_signif = kw_sig,
          Dunn_signif = cut(P.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                            labels = c("***", "**", "*", "ns"))
        ) %>%
        rename(
          Dunn_comparison = Comparison,
          Dunn_Z = Z,
          Dunn_p = P.unadj,
          Dunn_p_adj = P.adj
        ) %>%
        select(Variable, KW_ChiSq, KW_p, KW_signif,
               Dunn_comparison, Dunn_Z, Dunn_p, Dunn_p_adj, Dunn_signif)
      
      results <- bind_rows(results, dunn_result)
      
    } else {
      results <- bind_rows(results, data.frame(
        Variable = dv,
        KW_ChiSq = as.numeric(kw_test$statistic),
        KW_p = kw_test$p.value,
        KW_signif = kw_sig,
        Dunn_comparison = paste(levels(data[[group_var]])[1], "vs", levels(data[[group_var]])[2]),
        Dunn_Z = NA,
        Dunn_p = NA,
        Dunn_p_adj = NA,
        Dunn_signif = "ns"
      ))
    }
  }
  
  return(results)
}

# Run test
kruskal_dunn_results <- run_kruskal_dunn(dat_mtb, dv_vars)

# PLOTTING ----
create_kruskal_plot <- function(data, title, facet_titles, stat.test,
                                y_label, y_limits = NULL, y_breaks = NULL) {
  p <- ggplot(data, aes(x = age_group, y = value, fill = Variable)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, 
                 outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.3, size = 1.2, aes(shape = age_group)) +
    stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0,
                       hide.ns = FALSE, size = 6) +
    scale_fill_brewer(palette = "Dark2") +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    labs(x = NULL, y = y_label, title = title) +
    guides(fill = "none", alpha = "none") +
    facet_wrap(~ Variable, scales = "fixed", labeller = as_labeller(facet_titles)) +
    theme_bw() +
    theme(
      strip.text.x = element_text(size = 15, face = "bold"),
      plot.title = element_text(size = 20, face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      strip.placement = "outside",
      panel.spacing.x = unit(1.5, "lines"),
      panel.border = element_rect(color = "grey60")
    )
  return(p)
}

# Your selected variables
vars_to_plot <- c("Incoh12SS_4chan", "Coh12SS_4chan")

# Reshape to long format
long_df_subset <- dat_mtb %>%
  select(age_group, all_of(vars_to_plot)) %>%
  pivot_longer(cols = -age_group, names_to = "Variable", values_to = "value") %>%
  filter(!is.na(value))

long_df_subset$age_group <- factor(long_df_subset$age_group, levels = c("YNH (<=40 y)", "ONH (>40 y)"))

# Filter Dunn test results for these variables
stat.test_subset <- kruskal_dunn_results %>%
  filter(Variable %in% vars_to_plot) %>%
  mutate(
    group1 = "YNH (<=40 y)",
    group2 = "ONH (>40 y)",
    y.position = 3.5,
    p.signif = Dunn_signif
  )

facet_titles_subset <- c(
  "Incoh12SS_4chan" = "Incoherent",
  "Coh12SS_4chan" = "Coherent"
)

subset_plot <- create_kruskal_plot(
  data = long_df_subset,
  title = "A. Steady-State | 12-Tone",
  facet_titles = facet_titles_subset,
  stat.test = stat.test_subset,
  # y_label = "Latency (ms)",
  y_label = "Amplitude (µV)",
  y_limits = c(-3, 3.5),
  y_breaks = seq(-3, 3.5, by = 1)
)

print(subset_plot)

ggsave('./Figures/Kruskal-Wallis/SS12_KW.png', subset_plot, dpi = 500, width = 8, height = 5)


# Getting combined plots for paper ----

library(magick)
library(grid)
library(gridExtra)

img1 <- image_read("./Figures/Kruskal-Wallis/SS12_KW.png")
img2 <- image_read("./Figures/Kruskal-Wallis/SS20_KW.png")
# img3 <- image_read("./Figures/CMRBehInCoh_MW.png")

# Convert to raster grobs
g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
g2 <- rasterGrob(as.raster(img2), interpolate = TRUE)
# g3 <- rasterGrob(as.raster(img3), interpolate = TRUE)

# Define layout: g1 spans the full top row; g2 and g3 go on the second row
# layout <- rbind(c(1, 1),
#                 c(2, 3))
# 
# combined <- arrangeGrob(g1, g2, g3, layout_matrix = layout)

combined <- arrangeGrob(g1, g2, ncol=2)

# Save to file
png("./Figures/Kruskal-Wallis/SS_Combined_KW.png", width = 5, height = 1.5, units = "in", res = 500)
grid.draw(combined)
dev.off()

## Plots with spacers ----

library(grid)
library(gridExtra)
library(magick)

# Load images and convert to raster grobs
img1 <- image_read("./Figures/Mann-Whitney/CMRBeh_MW.png")
img2 <- image_read("./Figures/Mann-Whitney/CMRBehCoh_MW.png")
img3 <- image_read("./Figures/Mann-Whitney/CMRBehInCoh_MW.png")
# img4 <- image_read("./Figures/Kruskal-Wallis/GapP2AmplitudeCluster_KW.png")

g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
g2 <- rasterGrob(as.raster(img2), interpolate = TRUE)
g3 <- rasterGrob(as.raster(img3), interpolate = TRUE)
# g4 <- rasterGrob(as.raster(img4), interpolate = TRUE)

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
png("./Figures/Mann-Whitney/CMR_Combined.png", width = 3, height = 3, units = "in", res = 500)
grid.draw(combined)
dev.off()

