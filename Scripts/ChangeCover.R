################################################################################
################################################################################
#########################  FWF - Site Prep Study  ##############################
#########################    Change in Cover      ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2024 - 2025        ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", "plotrix", 
                      "ggsignif", "multcompView", "ggpubr", "rstatix", "labdsv",
                      "tables")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
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

# Orders years and treatments so that they display in same sequence in graphs #
df_cleaned$Month = factor(df_cleaned$Month, levels=c('November','April'))

# Convert Treatment to a factor with specified levels
df_cleaned$Treatment = factor(df_cleaned$Treatment,
                              levels=c('Control', 'Fill Sand', 
                                       'Preemergent','Soil Inversion'))

#Renames values in Treatment treatments for heat map later #
#Data$Treatment <- recode(Data$Treatment, 
#                    C ="Control", W = "Wiregrass", BH = "Broomsedge High", 
#                    BM = "Broomsedge Medium", LH = "Lovegrass High", 
#                    LM = "Lovegrass Medium")

#################### Species abundances ########################################
# Creates and joins  data month to make long data format #
Two_Abundance <- df_cleaned[which(df_cleaned$Month == "November"),]
Three_Abundance <- df_cleaned[which(df_cleaned$Month == "April"),]

Abundance_w <- full_join(Two_Abundance, Three_Abundance, 
                         by = c('Unique_ID', "Treatment", 'Species'))
Abundance_w = arrange(Abundance_w, Treatment)

# Turns NA values into zeros #
Abundance_w$Coverage.x <- ifelse(is.na(Abundance_w$Coverage.x), 0, 
                                 Abundance_w$Coverage.x)
Abundance_w$Coverage.y <- ifelse(is.na(Abundance_w$Coverage.y), 0, 
                                 Abundance_w$Coverage.y)

# Change abundance to reflect percentage change from (Year 1) to (Year 2)  #
Change_Abundance <- Abundance_w %>% 
  dplyr::select(Unique_ID, Treatment, Species, 
                Coverage.x, Coverage.y) %>%
  group_by(Unique_ID, Treatment, Species) %>% 
  mutate(Change_abundance = Coverage.y - Coverage.x)

##################################  COVER CAHNGES ##############################
BG = 
  Change_Abundance[which(Change_Abundance$Species == "Bare Ground"),]
BG<-as.data.frame(BG)
BG$Treatment<-factor(BG$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = BG)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_BG = BG %>% kruskal_test(Change_abundance ~ Treatment) %>%
  add_significance()
anova_BG

# Test for Significance and create letters
tukey_BG <- BG %>%
  dunn_test(Change_abundance ~ Treatment) %>%
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
  summarise(Change_abundance = max(Change_abundance) + 5) %>%
  left_join(letters_df, by = "Treatment")

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=BG)
tmp

# The plot code will now use the new 'label_data' and 'geom_text'
April_BG_Box <- ggplot(BG, aes(x = Treatment, y = Change_abundance)) +
  geom_boxplot(aes(fill = Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data,
            aes(x = Treatment, y = Change_abundance, label = letters),
            size = 8, fontface = "bold") +
  labs(subtitle = get_test_label(anova_BG, detailed = TRUE),
       caption = get_pwc_label(tukey_BG)) +
  scale_color_manual(labels = c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                     values = c("yellow", "orange", "#66CC00", "#CC66CC")) +
  scale_fill_manual(labels = c('Control', 'Fill', 'Preemergent','Soil Inversion'),
                    values = c("yellow", "orange", "#66CC00", "#CC66CC")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text = element_text(size = 16),
        axis.title.x = element_text(size = 15, face="bold", colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        strip.text.x = element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "", title = "1- to 6-month")
April_BG_Box

##################################  COVER CAHNGES ##############################
PN = 
  Change_Abundance[which(Change_Abundance$Species == "Paspalum notatum"),]
PN<-as.data.frame(PN)
PN$Treatment<-factor(PN$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = PN)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = PN %>% kruskal_test(Change_abundance ~ Treatment) %>%
  add_significance()
anova

lm(formula = Change_abundance ~ Treatment, PN)
tukey <- PN %>%
  dunn_test(Change_abundance ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=PN)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- PN %>%
  group_by(Treatment) %>%
  summarise(Change_abundance = max(Change_abundance) + 5) %>%
  left_join(letters_df, by = "Treatment")

April_PN_Change_Box =
  ggplot(PN, aes(x = Treatment, y = Change_abundance)) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Change_abundance, label = letters), size = 8, fontface = "bold") +
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
  labs(x = "", y = "Change in % coverage", title = "Change in % coverage of Paspalum notatum 1- to 6-month")
April_PN_Change_Box

##################################  COVER CAHNGES ##############################
OL = 
  Change_Abundance[which(Change_Abundance$Species == "Oenothera laciniata"),]
OL<-as.data.frame(OL)
OL$Treatment<-factor(OL$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = OL)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = OL %>% kruskal_test(Change_abundance ~ Treatment) %>%
  add_significance()
anova

lm(formula = Change_abundance ~ Treatment, OL)
tukey <- OL %>%
  dunn_test(Change_abundance ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=OL)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- OL %>%
  group_by(Treatment) %>%
  summarise(Change_abundance = max(Change_abundance) + 5) %>%
  left_join(letters_df, by = "Treatment")

April_OL_Change_Box =
  ggplot(OL, aes(x = Treatment, y = Change_abundance)) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Change_abundance, label = letters), size = 8, fontface = "bold") +
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
  labs(x = "", y = "Change in % Coverage", title = "Oenothera laciniata")
April_OL_Change_Box

################################################################################
################################################################################
#################################  Jan  ########################################
################################################################################
################################################################################

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

# Orders years and treatments so that they display in same sequence in graphs #
df_cleaned$Month = factor(df_cleaned$Month, levels=c('November','Jan'))

# Convert Treatment to a factor with specified levels
df_cleaned$Treatment = factor(df_cleaned$Treatment,
                              levels=c('Control', 'Fill Sand', 
                                       'Preemergent','Soil Inversion'))

#Renames values in Treatment treatments for heat map later #
#Data$Treatment <- recode(Data$Treatment, 
#                    C ="Control", W = "Wiregrass", BH = "Broomsedge High", 
#                    BM = "Broomsedge Medium", LH = "Lovegrass High", 
#                    LM = "Lovegrass Medium")

#################### Species abundances ########################################
# Creates and joins  data month to make long data format #
Two_Abundance <- df_cleaned[which(df_cleaned$Month == "November"),]
Three_Abundance <- df_cleaned[which(df_cleaned$Month == "Jan"),]

Abundance_w <- full_join(Two_Abundance, Three_Abundance, 
                         by = c('Unique_ID', "Treatment", 'Species'))
Abundance_w = arrange(Abundance_w, Treatment)

# Turns NA values into zeros #
Abundance_w$Coverage.x <- ifelse(is.na(Abundance_w$Coverage.x), 0, 
                                 Abundance_w$Coverage.x)
Abundance_w$Coverage.y <- ifelse(is.na(Abundance_w$Coverage.y), 0, 
                                 Abundance_w$Coverage.y)

# Change abundance to reflect percentage change from (Year 1) to (Year 2)  #
Change_Abundance <- Abundance_w %>% 
  dplyr::select(Unique_ID, Treatment, Species, 
                Coverage.x, Coverage.y) %>%
  group_by(Unique_ID, Treatment, Species) %>% 
  mutate(Change_abundance = Coverage.y - Coverage.x)

##################################  COVER CAHNGES ##############################
BG = 
  Change_Abundance[which(Change_Abundance$Species == "Bare Ground"),]
BG<-as.data.frame(BG)
BG$Treatment<-factor(BG$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = BG)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_BG = BG %>% kruskal_test(Change_abundance ~ Treatment) %>%
  add_significance()
anova_BG

# Test for Significance and create letters
tukey_BG <- BG %>%
  dunn_test(Change_abundance ~ Treatment) %>%
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
  summarise(Change_abundance = max(Change_abundance) + 5) %>%
  left_join(letters_df, by = "Treatment")

# The plot code will now use the new 'label_data' and 'geom_text'
Jan_BG_Box = ggplot(BG, aes(x = Treatment, y = Change_abundance)) +
  geom_boxplot(aes(fill = Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  # Use geom_text to add the letters from the new data frame
  geom_text(data = label_data, aes(x = Treatment, y = Change_abundance, 
                                   label = letters), size = 8, 
            fontface = "bold") +
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
  labs(x = "", y = "Bare ground % coverage change", title = "1- to 3-month")
Jan_BG_Box

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=BG)
tmp

##################################  COVER CAHNGES ##############################
PN = 
  Change_Abundance[which(Change_Abundance$Species == "Paspalum notatum"),]
PN<-as.data.frame(PN)
PN$Treatment<-factor(PN$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = PN)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = PN %>% kruskal_test(Change_abundance ~ Treatment) %>%
  add_significance()
anova

lm(formula = Change_abundance ~ Treatment, PN)
tukey <- PN %>%
  dunn_test(Change_abundance ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=PN)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- PN %>%
  group_by(Treatment) %>%
  summarise(Change_abundance = max(Change_abundance) + 5) %>%
  left_join(letters_df, by = "Treatment")

Jan_PN_Change_Box =
  ggplot(PN, aes(x = Treatment, y = Change_abundance)) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Change_abundance, label = letters), size = 8, fontface = "bold") +
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
  labs(x = "", y = "Change in % coverage", title = "Change in % coverage of Paspalum notatum 1- to 3-month")
Jan_PN_Change_Box

##################################  COVER CAHNGES ##############################
OL = 
  Change_Abundance[which(Change_Abundance$Species == "Oenothera laciniata"),]
OL<-as.data.frame(OL)
OL$Treatment<-factor(OL$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = OL)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova = OL %>% kruskal_test(Change_abundance ~ Treatment) %>%
  add_significance()
anova

lm(formula = Change_abundance ~ Treatment, OL)
tukey <- OL %>%
  dunn_test(Change_abundance ~ Treatment) %>%
  add_significance() %>%
  add_xy_position()
tukey

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=OL)
tmp

# Corrected plot code
p_values <- tukey$p.adj
names(p_values) <- paste(tukey$group1, tukey$group2, sep = "-")
letters <- multcompLetters(p_values)$Letters
letters_df <- data.frame(Treatment = names(letters), letters = letters)
label_data <- OL %>%
  group_by(Treatment) %>%
  summarise(Change_abundance = max(Change_abundance) + 5) %>%
  left_join(letters_df, by = "Treatment")

Jan_OL_Change_Box =
  ggplot(OL, aes(x = Treatment, y = Change_abundance)) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3,
             position = position_jitterdodge(), alpha = 0.7) +
  geom_text(data = label_data, aes(x = Treatment, y = Change_abundance, label = letters), size = 8, fontface = "bold") +
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
  labs(x = "", y = "Change in % coverage", title = "Oenothera laciniata")
Jan_OL_Change_Box

################################################################################
################################################################################

ChangeCoverage_all = ggarrange(Jan_BG_Box, April_BG_Box, nrow = 1, ncol = 2)
ChangeCoverage_all

ggsave("Figures/ChangeCoverage_All.tiff", dpi = 100,
       width = 14, height = 8)
