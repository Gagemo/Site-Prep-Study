################################################################################
################################################################################
#########################  FWF - Site Prep Study  ##############################
#########################           Cover         ##############################
#########################  University of Florida  ##############################
#########################      Gage LaPierre      ##############################
#########################      2024 - 2025        ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014")

#########################      Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", "plotrix",
                      "ggsignif", "multcompView", "ggpubr", "rstatix", "labdsv",
                      "tables", "ggtext")
new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################      Loads Packages      ##############################
library(tidyverse)
library(vegan)
library(labdsv)
library(agricolae)
library(extrafont)
library(ggsignif)
library(multcompView)
library(ggpubr)
library(plotrix)
library(rstatix)
library(tables)
library(ggtext)

########################## Read in Data ########################################
Data = read.csv("Data/Post-Treatment Data.csv")

# Transform data from wide to long format
df_long <- Data %>%
  pivot_longer(
    cols = c('X0m', 'X3m', 'X6m', 'X9m', 'X12m', 'X15m', 'X18m',
             'X21m', 'X24m', 'X27m'),
    names_to = 'Quadrat_ID',
    values_to = 'Coverage'
  ) %>%
  # Remove 'X' and 'm' from Quadrat_ID and create Quadrat_ID_Numeric
  mutate(Quadrat_ID_Numeric = str_replace(Quadrat_ID, 'X', '') %>%
           str_replace('m', '')) %>%
  # Create the unique ID, now including 'Month'
  mutate(Unique_ID = paste(Block, Plot, Quadrat_ID_Numeric, sep = '-'))

# Filter out *only* blank entries (empty strings) in Coverage
# We will handle "NA" and other NA values by converting them to 0
df_cleaned <- df_long %>%
  filter(Coverage != "")

# Convert Coverage to numeric. This will turn any non-numeric strings (including "NA" if still present) into actual R `NA` values.
df_cleaned$Coverage <- as.numeric(df_cleaned$Coverage)

# Replace all NA values in Coverage with 0
df_cleaned <- df_cleaned %>%
  mutate(Coverage = replace_na(Coverage, 0))

# Reclassify coverage data (CV) from 1-10 scale to percent scale
# Explicitly ensure the output of case_when is numeric
df_cleaned <- df_cleaned %>%
  mutate(Coverage = as.numeric(case_when(
    Coverage == 10 ~ 97.5,
    Coverage == 1 ~ 0.1,
    Coverage == 2 ~ 0.5,
    Coverage == 3 ~ 1.5,
    Coverage == 4 ~ 3.5,
    Coverage == 5 ~ 7.5,
    Coverage == 6 ~ 17.5,
    Coverage == 7 ~ 37.5,
    Coverage == 8 ~ 62.5,
    Coverage == 9 ~ 85,
    Coverage == 0 ~ 0,
    TRUE ~ Coverage)))

# Convert Treatment to a factor with specified levels
df_cleaned$Treatment = factor(df_cleaned$Treatment,
                              levels=c('Control', 'Fill Sand',
                                       'Preemergent','Soil Inversion'))

#################### Species abundances ########################################
# Creates and joins  data month to make long data format #
November <- filter(df_cleaned, Month == "November")
Jan <- filter(df_cleaned, Month == "Jan")
April <- filter(df_cleaned, Month == "April")

##################################  November ###################################
BG = filter(November, Species == "Bare Ground")

# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = BG)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_BG = BG %>% kruskal_test(Coverage ~ Treatment) %>%
  add_significance()
anova_BG

lm(formula = Coverage ~ Treatment, BG)
tukey_BG <- BG %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey_BG

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=BG)
tmp

# Corrected plot code
p_values <- tukey_BG$p.adj
names(p_values) <- paste(tukey_BG$group1, tukey_BG$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- BG %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

Nov_BG_Box =
  ggplot(BG, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, label = letters), size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova_BG, detailed = TRUE),
       caption = get_pwc_label(tukey_BG)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        strip.text.x =
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% coverage", title = "Bare Ground")
Nov_BG_Box

##################################  Bahia ##############################
PN = filter(November, Species == "Paspalum notatum")

# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = PN)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = PN %>% kruskal_test(Coverage ~ Treatment) %>%
  add_significance()
anova

lm(formula = Coverage ~ Treatment, PN)
tukey <- PN %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=PN)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- PN %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

Nov_PN_Box =
  ggplot(PN, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, label = letters), 
            size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova, detailed = TRUE),
       caption = get_pwc_label(tukey)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        strip.text.x =
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% coverage", title = "Paspalum notatum")
Nov_PN_Box


################################################################################
##################################  Jan ########################################
################################################################################

BG = filter(Jan, Species == "Bare Ground")

# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = BG)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_BG = BG %>% kruskal_test(Coverage ~ Treatment) %>%
  add_significance()
anova_BG

lm(formula = Coverage ~ Treatment, BG)
tukey_BG <- BG %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey_BG

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=BG)
tmp

# Corrected plot code
p_values <- tukey_BG$p.adj
names(p_values) <- paste(tukey_BG$group1, tukey_BG$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- BG %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

Jan_BG_Box =
  ggplot(BG, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, label = letters), size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova_BG, detailed = TRUE),
       caption = get_pwc_label(tukey_BG)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        strip.text.x =
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% coverage", title = "Bare Ground")
Jan_BG_Box

##################################  Bahia ##############################
PN = filter(Jan, Species == "Paspalum notatum")

# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = PN)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = PN %>% kruskal_test(Coverage ~ Treatment) %>%
  add_significance()
anova

lm(formula = Coverage ~ Treatment, PN)
tukey <- PN %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=PN)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- PN %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

Jan_PN_Box =
  ggplot(PN, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, label = letters), size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova, detailed = TRUE),
       caption = get_pwc_label(tukey)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, colour = "black"),
        text=element_text(size=16),
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.ticks.x=element_blank(), 
        legend.position = "none") +
  labs(x = "", y = "% coverage", title = "Paspalum notatum")
Jan_PN_Box

#################### Oenothera laciniata #######################################

OL = filter(Jan, Species == "Oenothera laciniata")
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = OL)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = OL %>% kruskal_test(Coverage ~ Treatment) %>%
  add_significance()
anova

lm(formula = Coverage ~ Treatment, OL)
tukey <- OL %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=OL)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- OL %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

Jan_OL_Box =
  ggplot(OL, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, label = letters), size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova, detailed = TRUE)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = ggtext::element_markdown(hjust = 0.5, color = "black"),
        text=element_text(size=16),
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% coverage", title = "**3-month**<br>Oenothera laciniata") +
  ylim(0,70)
Jan_OL_Box
Jan_OL_Box

#################### Richardia spp. #######################################

RS = filter(Jan, Species == "Richardia spp.")
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = RS)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = RS %>% kruskal_test(Coverage ~ Treatment) %>%
  add_significance()
anova

lm(formula = Coverage ~ Treatment, RS)
tukey <- RS %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=RS)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- RS %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

Jan_RS_Box =
  ggplot(RS, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, 
                                   label = letters), size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova, detailed = TRUE)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, colour = "black"),
        text=element_text(size=16),
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% coverage", title = "Richardia spp.")
Jan_RS_Box

################################################################################
##################################  April ######################################
################################################################################

BG = filter(April, Species == "Bare Ground")

# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = BG)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance and create letters
tukey_BG <- BG %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()

# Get the p-values and convert to compact letters
p_values <- tukey_BG$p.adj

# Create the names using the group1 and group2 columns
names(p_values) <- paste(tukey_BG$group1, tukey_BG$group2, sep = "-")

# Correct way to get the letters and create the data frame
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)

# Create a new summary table for plotting the letters
label_data <- BG %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

# The plot code will now use the new 'label_data' and 'geom_text'
April_BG_Box = ggplot(BG, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill = Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  # Use geom_text to add the letters from the new data frame
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, label = letters), size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova_BG, detailed = TRUE),
       caption = get_pwc_label(tukey_BG)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text = element_text(size = 16),
        axis.title.x = element_text(size = 15, face="bold", colour = "black"),
        axis.title.y = element_text(size = 15, face="bold", colour = "black"),
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        strip.text.x = element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% coverage", title = "Bare Ground")
April_BG_Box

##################################  Bahia ##############################
PN = filter(April, Species == "Paspalum notatum")

# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = PN)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = PN %>% kruskal_test(Coverage ~ Treatment) %>%
  add_significance()
anova

lm(formula = Coverage ~ Treatment, PN)
tukey <- PN %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=PN)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- PN %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

April_PN_Box =
  ggplot(PN, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, label = letters), size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova, detailed = TRUE),
       caption = get_pwc_label(tukey)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, colour = "black"),
        text=element_text(size=16),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% coverage", title = "Paspalum notatum ")
April_PN_Box

#################### Oenothera laciniata #######################################

OL = filter(April, Species == "Oenothera laciniata")
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = OL)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = OL %>% kruskal_test(Coverage ~ Treatment) %>%
  add_significance()
anova

lm(formula = Coverage ~ Treatment, OL)
tukey <- OL %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=OL)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- OL %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

April_OL_Box =
  ggplot(OL, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, label = letters), size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova, detailed = TRUE)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = ggtext::element_markdown(hjust = 0.5, color = "black"),
        text=element_text(size=16),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% coverage", title = "**6-month**<br>Oenothera laciniata") +
  ylim(0,70)
April_OL_Box

#################### Richardia spp. #######################################

RS = filter(April, Species == "Richardia spp.")
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = RS)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = RS %>% kruskal_test(Coverage ~ Treatment) %>%
  add_significance()
anova

lm(formula = Coverage ~ Treatment, RS)
tukey <- RS %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=RS)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- RS %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

April_RS_Box =
  ggplot(RS, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, 
                                   label = letters), size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova, detailed = TRUE)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, colour = "black"),
        text=element_text(size=16),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% coverage", title = "Richardia spp. ") +
  ylim(0, 45)
April_RS_Box

#################### Indigofera hirsuta #######################################

IH = filter(April, Species == "Indigofera hirsuta")
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = IH)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = IH %>% kruskal_test(Coverage ~ Treatment) %>%
  add_significance()
anova

lm(formula = Coverage ~ Treatment, IH)
tukey <- IH %>%
  dunn_test(Coverage ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=IH)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- IH %>%
  group_by(Treatment) %>%
  summarise(Coverage = max(Coverage) + 5) %>%
  left_join(letters_df, by = "Treatment")

April_IH_Box =
  ggplot(IH, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Coverage, 
                                   label = letters), size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova, detailed = TRUE),
       caption = get_pwc_label(tukey)) +
  scale_color_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels=c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values=c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        strip.text.x =
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% coverage", title = "Indigofera hirsuta")
April_IH_Box

################################################################################

Coverage = ggarrange(Jan_OL_Box, April_OL_Box,
                     Jan_RS_Box, April_RS_Box, 
                     Jan_PN_Box, April_PN_Box,
                     nrow = 3, ncol = 2)
Coverage
ggsave("Figures/Coverage.png", width = 14, height = 12, dpi=600)
