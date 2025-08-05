# Master Thesis 
# Loading Short LEA Peptides on Different Nanoparticle Systems: Characterization and Impact on Plant Stress Tolerance
# Biological Functions Engineering, Kyushu Institute of Technology, Japan
# Author: Ahmed Mustafa Kamel Saber
# Supervisor: Dr. Shinya Ikeno
# The data used in this script is provided upon formal requests to 
# [saber.ahmed-mustafa818@mail.kyutech.jp] or [ahmedmkamel.research@gmail.com] and
# [ikeno@life.kyutech.jp]
# Required Packages ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(emmeans)
library(multcompView)
library(multcomp)

# Bar plots ----
# Chlorophyll-CS-LEA-K-100----
df_100 <- read.csv("barplot-data-100.csv", header = TRUE)
# Rename first column to "Sample"
colnames(df_100)[1] <- "Sample"

# Pivot to long format
df_long_100 <- pivot_longer(df_100, cols = -Sample, names_to = "Condition", values_to = "Value")
trial_100 <- read.csv("df_long_100.csv")
# Basic grouped bar plot
ggplot(df_long_100, aes(x = Condition, y = Value, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Control (Arabidopsis)" = "gray40",
                               "Recombinant LEA-K (Arabidopsis)" = "black")) +
  theme_minimal() +
  labs(title = "Chlorophyll Levels (Chitosan-LEA-K-100 (µg/mL))",
       x = "Condition", y = "Chlorophyll (mg/gFW)")
# To calculate the mean and SD values

colnames(trial_100)[3] <- "Mean"
ggplot(trial_100, aes(x = Condition, y = Mean, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("Wild-type (Arabidopsis)" = "gray40",
                               "Recombinant LEA-K (Arabidopsis)" = "black")) +
  labs(y = "Chlorophyll (mg/gFW)", title = "Chlorophyll Levels (Chitosan-LEA-K-100 (µg/mL))")+
  theme_minimal()

# Chlorophyll-CS-LEA-K-250----
df_250 <- read.csv("barplot-data-250.csv", header = TRUE)
# Rename first column to "Sample"
colnames(df_250)[1] <- "Sample"

# Pivot to long format
df_long_250 <- pivot_longer(df_250, cols = -Sample, names_to = "Condition", values_to = "Value")
write.csv(df_long_250, "df_long_250.csv")
trial_250 <- read.csv("df_long_250.csv", header = TRUE)
# Basic grouped bar plot
ggplot(df_long_250, aes(x = Condition, y = Value, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Control (Arabidopsis)" = "gray40",
                               "Recombinant LEA-K (Arabidopsis)" = "black")) +
  theme_minimal() +
  labs(title = "Chlorophyll Levels (Chitosan-LEA-K-250 (µg/mL))",
       x = "Condition", y = "Chlorophyll (mg/gFW)")
# To calculate the mean and SD values

#colnames(trial_250)[3] <- "Mean"
ggplot(trial_250, aes(x = Condition, y = Mean, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("Wild-type (Arabidopsis)" = "gray40",
                               "Recombinant LEA-K (Arabidopsis)" = "black")) +
  labs(y = "Chlorophyll (mg/gFW)", title = "Chlorophyll Levels (Chitosan-LEA-K-250 (µg/mL))")+
  theme_minimal()

# Chlorophyll-LDH-LEA-Nanobiocomposites----
df_LDH <- read.csv("barplot-data-LDH.csv", header = TRUE)

# Pivot to long format
df_long_LDH <- pivot_longer(df_LDH, cols = -Sample, names_to = "Condition", values_to = "Value")
write.csv(df_long_LDH,"df_long_LDH.csv")
trial_LDH <- read.csv("df_long_LDH.csv")
# Basic grouped bar plot
ggplot(df_long_LDH, aes(x = Condition, y = Value, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), show.legend = TRUE) +
  scale_fill_manual(
    values = c("Wild-type (Arabidopsis)" = "black")
  ) +
  theme_minimal() +
  labs(
    title = "Chlorophyll Levels (LDH-LEA Nanobiocomposites)",
    x = "Condition", y = "Chlorophyll (mg/gFW)",
    fill = "Wild-type (Arabidopsis)"  # Legend title
  )

# To calculate the mean and SD values
colnames(trial_LDH)[3] <- "Mean"
ggplot(trial_LDH, aes(x = Condition, y = Mean, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7, show.legend = TRUE) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("Wild-type (Arabidopsis)" = "gray40")) +
  labs(y = "Chlorophyll (mg/gFW)", 
       title = "Chlorophyll Levels (LDH-LEA Nanobiocomposites) on Wild-type Arabidopsis")+
  theme_minimal()

# Statistical tests ----
# Run Two-way ANOVA
chlorophyll_data_100 <- read.csv("chlorophyll-data-100.csv", header = TRUE)

table(chlorophyll_data_100$Condition, chlorophyll_data_100$Treatment)
# Ensure categorical variables are factors
str(chlorophyll_data_100)
chlorophyll_data_100$Treatment <- as.factor(chlorophyll_data_100$Treatment)
chlorophyll_data_100$Condition <- as.factor(chlorophyll_data_100$Condition)
chlorophyll_data_100$Chlorophyll <- as.numeric(as.character(chlorophyll_data_100$Chlorophyll))
# Two-way ANOVA (replicates are automatically accounted for)

anova_simple_100 <- aov(Chlorophyll ~ Treatment + Condition, data = chlorophyll_data_100)
summary(anova_simple_100)
anova_data <- anova(anova_simple_100)
View(anova_data)
emmeans(anova_simple_100, pairwise ~ Treatment | Condition)

# Tukey post-hoc test
library(emmeans)
emmeans(anova_simple_100, pairwise ~ Treatment)
emmeans(anova_simple_100, pairwise ~ Condition)
emmeans(anova_simple_100, pairwise ~ Treatment | Condition)
emmeans(anova_simple_100, pairwise ~ Condition | Treatment)
emmeans(anova_simple_100, pairwise ~ Treatment * Condition)
#emmeans(anova_result, pairwise ~ Chlorophyll | Treatment)

results <- emmeans(anova_simple_100, pairwise ~ Treatment | Condition)
summary(results$contrasts)  # Pairwise comparisons
summary(results$emmeans)    # Group means with SE and CI
write.csv(as.data.frame(summary(results$contrasts)), "pairwise_comparisons.csv")
write.csv(as.data.frame(summary(results$emmeans)), "group_means.csv")

# Group Letters (Compact Letter Display)
library(multcompView)
library(multcomp)

# For each condition
cld_letters <- cld(results$emmeans, by = "Condition", Letters = letters)
cld_letters <- cld(results$emmeans, Letters = letters)
cld_letters

library(ggplot2)

# Convert emmeans result to data frame
df_plot <- as.data.frame(results$emmeans)

# Plot
ggplot(df_plot, aes(x = Treatment, y = emmean, fill = Treatment)) +
  geom_col(position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  facet_wrap(~Condition) +
  labs(title = "Estimated Marginal Means of Chlorophyll",
       y = "Estimated Mean Chlorophyll",
       x = "Treatment") +
  theme_minimal() +
  theme(legend.position = "none")

# Average across all conditions
combined_emm <- emmeans(anova_simple_100, ~ Treatment)
combined_cld <- cld(combined_emm, Letters = letters)

library(ggplot2)

ggplot(combined_cld, aes(x = Treatment, y = emmean, fill = Treatment)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  geom_text(aes(label = .group, y = emmean + SE + 0.05), vjust = 0) +
  labs(title = "Chlorophyll Levels Across All Treatments",
       y = "Estimated Mean Chlorophyll",
       x = "Treatment") +
  scale_fill_manual(values = c("Control (Arabidopsis)" = "gray40",
                               "Recombinant LEA-K (Arabidopsis)" = "black")) +
  theme_minimal() +
  theme(legend.position = "none")

# Grouped Bar chart
ggplot(cld_letters, aes(x = Treatment, y = emmean, fill = Condition)) +
  geom_col(position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.2, position = position_dodge(0.7)) +
  geom_text(aes(label = .group, y = emmean + SE + 0.05),
            position = position_dodge(0.7), vjust = 0) +
  labs(title = "Chlorophyll Levels by Treatment and Condition",
       y = "Estimated Mean Chlorophyll",
       x = "Treatment") +
  scale_fill_manual(values = c("Control (Arabidopsis)" = "gray40",
                               "Recombinant LEA-K (Arabidopsis)" = "black")) +
  theme_minimal()

# Define p values and add it to the bar plot----

# Read your CSV
df_pval <- read.csv("pairwise_comparisons.csv")

# Define significance stars
df_pval$p.signif <- cut(df_pval$p.value,
                        breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, 1),
                        labels = c("***", "**", "*", ".", "ns"))

# Split the contrast column to get group1 and group2
contrast_split <- strsplit(as.character(df_pval$contrast), " - ")
df_pval$group1 <- sapply(contrast_split, `[`, 1)
df_pval$group2 <- sapply(contrast_split, `[`, 2)

# Optional: define a y-axis position for the asterisks
df_pval$y.position <- seq(1.2, 3.0, length.out = nrow(df_pval))  # customize for your plot

# Final data frame for stat_pvalue_manual
stat.test <- df_pval[, c("group1", "group2", "p.value", "p.signif", "y.position")]
names(stat.test)[3] <- "p.adj"  # ggpubr expects p.adj

library(ggpubr)

ggplot(trial_100, aes(x = Condition, y = Mean, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("Wild-type (Arabidopsis)" = "gray40",
                               "Recombinant LEA-K (Arabidopsis)" = "black")) +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01) +
  labs(y = "Chlorophyll (mg/gFW)", title = "Chlorophyll Levels (Chitosan-LEA-K-100 (µg/mL))")+
  theme_minimal()

ggplot(df_summary, aes(x = Sample, y = Mean, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                width = 0.2, position = position_dodge()) +
  scale_fill_manual(values = c("Wild-type" = "gray40", "LEA-K" = "black")) +

  theme_minimal()
