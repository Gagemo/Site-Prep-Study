################################################################################
################################################################################
#########################    Site Prep Study      ##############################
#########################  Seedbank Germination   ##############################
######################### University of Florida   ##############################
#########################    Gage LaPierre        ##############################
#########################     2024 - 2025         ##############################
################################################################################
################################################################################

######################### Clears Environment & History #########################
rm(list = ls(all = TRUE))
cat("\014")

######################### Installs Packages #####################################
list.of.packages <- c(
  "tidyverse", "car", "ggpubr", "rstatix",
  "emmeans", "multcompView", "multcomp"
)

new.packages <- list.of.packages[
  !(list.of.packages %in% installed.packages()[, "Package"])
]

if (length(new.packages)) {
  install.packages(new.packages)
}

######################### Loads Packages ########################################
library(tidyverse)
library(car)
library(ggpubr)
library(rstatix)
library(emmeans)
library(multcompView)
library(multcomp)

######################### File and Output Settings ##############################
file_path <- "Data/Seedbank.csv"

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Tables", showWarnings = FALSE, recursive = TRUE)

######################### Read in Data ##########################################
Data <- read.csv(
  file_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = TRUE
)

######################### Prepare Data ##########################################
required_columns <- c(
  "Month", "Block", "Plot", "Sample", "Treatment"
)

missing_columns <- setdiff(required_columns, names(Data))

if (length(missing_columns) > 0) {
  stop(
    "The dataset is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Species columns are every column after the metadata columns.
species_columns <- setdiff(
  names(Data),
  c("Month", "Block", "Plot", "Sample", "Treatment")
)

if (length(species_columns) == 0) {
  stop("No species columns were detected.")
}

# Convert blank, NA, N/A, and "." entries to zero and make species columns numeric.
Data <- Data %>%
  mutate(
    across(
      all_of(species_columns),
      ~ {
        x <- trimws(as.character(.x))
        x[
          is.na(x) |
            x == "" |
            toupper(x) %in% c("NA", "N/A", ".")
        ] <- "0"
        suppressWarnings(as.numeric(x))
      }
    )
  )

if (any(is.na(Data[, species_columns]))) {
  stop(
    "One or more species values could not be converted to numeric."
  )
}

# Ensure treatment and survey factors are in the desired order.
Data$Treatment <- factor(
  Data$Treatment,
  levels = c(
    "Control",
    "Fill Sand",
    "Preemergent",
    "Soil Inversion"
  )
)

Data$Month <- factor(
  Data$Month,
  levels = c(
    "Pre-Treatment",
    "Post-Treatment"
  )
)

if (any(is.na(Data$Treatment))) {
  stop(
    "One or more treatment names do not match the expected factor levels."
  )
}

if (any(is.na(Data$Month))) {
  stop(
    "One or more Month values do not match the expected survey names."
  )
}

######################### Calculate Germination per Sample ######################
# Sum all emerged seedlings across species for each seedbank sample.
Data <- Data %>%
  rowwise() %>%
  mutate(
    TotalGermination = sum(
      c_across(all_of(species_columns)),
      na.rm = TRUE
    )
  ) %>%
  ungroup()

######################### Check Sample Replication ##############################
sample_check <- Data %>%
  group_by(
    Month,
    Block,
    Plot,
    Treatment
  ) %>%
  summarise(
    Number_of_Samples = n_distinct(Sample),
    Number_of_Rows = n(),
    .groups = "drop"
  )

if (any(
  sample_check$Number_of_Samples != 4 |
  sample_check$Number_of_Rows != 4
)) {
  write.csv(
    sample_check,
    "Tables/Seedbank_Germination_Sample_Check.csv",
    row.names = FALSE
  )
  
  stop(
    "One or more plot-by-survey combinations do not contain exactly four samples."
  )
}

######################### Average Samples Within Plots ##########################
# Treatment was applied to plots, so the four seedbank samples are subsamples.
# Average total germination across the four samples within each plot.

Plot_Germination <- Data %>%
  group_by(
    Month,
    Block,
    Plot,
    Treatment
  ) %>%
  summarise(
    MeanGermination = mean(
      TotalGermination,
      na.rm = TRUE
    ),
    Number_of_Samples = n_distinct(Sample),
    .groups = "drop"
  )

write.csv(
  Plot_Germination,
  "Tables/Seedbank_Germination_Plot_Level_Data.csv",
  row.names = FALSE
)

######################### Color Palette #########################################
treat_colors <- c(
  "Control" = "yellow",
  "Fill Sand" = "orange",
  "Preemergent" = "#66CC00",
  "Soil Inversion" = "#CC66CC"
)

######################### Analysis and Plot Function ############################
analyze_germination <- function(
    df,
    survey_name,
    show_y_axis = TRUE
) {
  
  survey_data <- df %>%
    filter(Month == survey_name) %>%
    droplevels()
  
  # Include Block because the experiment used a randomized complete block design.
  model <- lm(
    MeanGermination ~ Block + Treatment,
    data = survey_data
  )
  
  # Type II ANOVA.
  anova_results <- car::Anova(
    model,
    type = 2
  )
  
  cat("\n--------------------------------------------------\n")
  cat("ANOVA Table for ", survey_name, " Seed Germination\n", sep = "")
  cat("--------------------------------------------------\n")
  print(anova_results)
  
  # Assumption checks.
  shapiro_results <- shapiro.test(residuals(model))
  levene_results <- car::leveneTest(
    MeanGermination ~ Treatment,
    data = survey_data
  )
  
  cat("\nShapiro-Wilk test of model residuals:\n")
  print(shapiro_results)
  
  cat("\nLevene test for homogeneity of variance:\n")
  print(levene_results)
  
  # Estimated marginal means and Tukey-adjusted compact letter display.
  emm <- emmeans(
    model,
    ~ Treatment
  )
  
  tukey_cld <- multcomp::cld(
    emm,
    adjust = "tukey",
    Letters = letters,
    decreasing = TRUE,
    alpha = 0.05
  ) %>%
    as.data.frame() %>%
    mutate(
      label = trimws(.group)
    )
  
  # Place letters above the highest plot-level observation in each treatment.
  max_y_values <- survey_data %>%
    group_by(Treatment) %>%
    summarise(
      max_y = max(
        MeanGermination,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  tukey_cld <- left_join(
    tukey_cld,
    max_y_values,
    by = "Treatment"
  )
  
  # Extract the treatment p-value for the plot subtitle.
  anova_table <- as.data.frame(anova_results)
  
  treatment_p <- anova_table[
    "Treatment",
    "Pr(>F)"
  ]
  
  subtitle_text <- paste0(
    "F = ",
    round(anova_table["Treatment", "F value"], 2),
    ", p ",
    ifelse(
      treatment_p < 0.001,
      "< 0.001",
      paste0("= ", round(treatment_p, 3))
    )
  )
  
  p <- ggplot(
    survey_data,
    aes(
      x = Treatment,
      y = MeanGermination,
      fill = Treatment
    )
  ) +
    geom_boxplot(
      alpha = 0.5,
      outlier.shape = NA,
      show.legend = FALSE
    ) +
    geom_jitter(
      color = "black",
      size = 3,
      alpha = 0.7,
      width = 0.2,
      show.legend = FALSE
    ) +
    geom_text(
      data = tukey_cld,
      aes(
        x = Treatment,
        y = max_y + 8,
        label = label
      ),
      inherit.aes = FALSE,
      vjust = 0,
      hjust = 0.5,
      size = 8,
      fontface = "bold"
    ) +
    scale_fill_manual(values = treat_colors) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(
      x = "",
      y = "Mean germination count per sample",
      title = survey_name,
      subtitle = subtitle_text
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        colour = "black"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        colour = "black"
      ),
      text = element_text(size = 20),
      axis.title.x = element_text(
        size = 20,
        face = "bold",
        colour = "black"
      ),
      axis.title.y = element_text(
        size = 20,
        face = "bold",
        colour = "black"
      ),
      axis.text.x = element_text(
        size = 20,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(
        size = 20,
        face = "bold",
        color = "black"
      ),
      legend.text = element_text(
        size = 20,
        face = "bold",
        color = "black"
      ),
      legend.position = "bottom"
    )
  
  if (!show_y_axis) {
    p <- p +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
      )
  }
  
  return(
    list(
      Model = model,
      ANOVA = anova_results,
      Shapiro = shapiro_results,
      Levene = levene_results,
      EMM = emm,
      CLD = tukey_cld,
      Plot = p
    )
  )
}

######################### Pretreatment Analysis #################################
pre_results <- analyze_germination(
  Plot_Germination,
  survey_name = "Pre-Treatment",
  show_y_axis = TRUE
)

pre <- pre_results$Plot

######################### Post-treatment Analysis ###############################
post_results <- analyze_germination(
  Plot_Germination,
  survey_name = "Post-Treatment",
  show_y_axis = FALSE
)

post <- post_results$Plot

######################### Save Statistical Tables ###############################
write.csv(
  as.data.frame(pre_results$ANOVA),
  "Tables/Seedbank_Germination_Pretreatment_ANOVA.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(post_results$ANOVA),
  "Tables/Seedbank_Germination_Posttreatment_ANOVA.csv",
  row.names = TRUE
)

write.csv(
  pre_results$CLD,
  "Tables/Seedbank_Germination_Pretreatment_Tukey.csv",
  row.names = FALSE
)

write.csv(
  post_results$CLD,
  "Tables/Seedbank_Germination_Posttreatment_Tukey.csv",
  row.names = FALSE
)

######################### Arrange Figures #######################################
seedbank_germ_count <- ggarrange(
  pre,
  post,
  nrow = 1,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom",
  labels = c("A", "B")
)

final_figure <- annotate_figure(
  seedbank_germ_count,
  bottom = text_grob(
    "Letters indicate Tukey-adjusted pairwise comparisons.",
    size = 15
  )
)

print(final_figure)

######################### Save Figure ###########################################
ggsave(
  "Figures/seedbank_germ_count.tiff",
  plot = final_figure,
  dpi = 300,
  width = 18,
  height = 8,
  units = "in",
  compression = "lzw"
)