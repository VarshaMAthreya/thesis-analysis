library(dplyr)
library(tidyr)
library(readr)

# Load your CSV file
mtb <- read.csv('./ThesisAnalysis_All_4.3.25.csv')

# Divide by age group
mtb_age <- mtb %>%
  mutate(age_group = case_when(
    Age <= 40 ~ "YNH (\u2264 40 y)",
    Age >= 41 ~ "ONH (\u2265 41 y)",
    TRUE ~ NA_character_
  ))

# Order of the age groups
Ages_order <- c("YNH (\u2264 40 y)", "ONH (\u2265 41 y)")
mtb_age$age_group <- factor(mtb_age$age_group, levels = Ages_order)

### -------------------------------------------
### SHAPIRO-WILK: Exclude Subject, Age, Gender
### -------------------------------------------

exclude_vars <- c("Subject", "Age", "Gender")
test_vars <- mtb_age %>%
  select(-all_of(exclude_vars)) %>%
  select(where(is.numeric))

shapiro_results <- sapply(test_vars, function(x) {
  if (length(na.omit(x)) >= 3 && length(na.omit(x)) <= 5000) {
    shapiro.test(na.omit(x))$p.value
  } else {
    NA
  }
})

shapiro_df <- data.frame(
  Variable = names(shapiro_results),
  Shapiro_p = round(shapiro_results, 4),
  Normal = ifelse(shapiro_results > 0.05, "Yes", "No")
) %>% arrange(Shapiro_p)

print(shapiro_df)

### -------------------------------------------
### OUTLIER DETECTION: Exclude Subject, Age, Gender
### -------------------------------------------

# Keep only relevant columns (exclude Age & Gender)
numeric_for_outlier <- mtb_age %>%
  select(Subject, age_group, where(is.numeric)) %>%
  select(-Age)

# Pivot to long format
long_data <- numeric_for_outlier %>%
  pivot_longer(cols = -c(Subject, age_group), names_to = "variable", values_to = "value")

# Flag outliers using z-score
long_flagged <- long_data %>%
  group_by(age_group, variable) %>%
  mutate(z = scale(value), is_outlier = abs(z) >= 2.5) %>%
  ungroup()

# Optional: view outliers
outliers_by_subject <- long_flagged %>% filter(is_outlier == TRUE)
print(outliers_by_subject)

### -------------------------------------------
###Q-Q Plots

library(ggpubr)

# Select numeric variables excluding Subject, Age, Gender
qq_vars <- mtb_age %>%
  select(where(is.numeric)) %>%
  select(-Age)

# Loop through each variable and plot Q-Q plots
for (varname in names(qq_vars)) {
  p <- ggqqplot(qq_vars, x = varname,
                title = paste("Q-Q Plot:", varname),
                color = "black") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  print(p)}

### -------------------------------------------
### CLEAN AND SAVE

# Replace outliers with NA
long_cleaned <- long_flagged %>%
  mutate(value_cleaned = as.numeric(ifelse(is_outlier, NA, value)))  # force numeric

# Pivot back to wide format
cleaned_wide <- long_cleaned %>%
  select(Subject, variable, value_cleaned) %>%
  pivot_wider(names_from = variable, values_from = value_cleaned)

# Rejoin with Age, Gender, and age_group from original
final_cleaned <- mtb_age %>%
  select(Subject, Age, Gender, age_group) %>%
  left_join(cleaned_wide, by = "Subject")

# Save cleaned data
write_csv(final_cleaned, "ThesisAnalysis_All_WithoutOutliers_5.30.csv")
