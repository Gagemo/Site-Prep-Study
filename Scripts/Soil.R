################################################################################
################################################################################
#########################  FWF - Site Prep Study  ##############################
#########################      Soil Chemistry     ##############################
#########################  University of Florida  ##############################
#########################      Gage LaPierre      ##############################
#########################      2024 - 2025        ##############################
################################################################################
################################################################################

######################### Clears Environment & History #########################
rm(list = ls(all = TRUE))
cat("\014")

######################### Installs Packages #####################################
list.of.packages <- c(
  "tidyverse", "car", "ggpubr", "emmeans",
  "multcomp", "multcompView", "lme4", "lmerTest"
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
library(emmeans)
library(multcomp)
library(multcompView)
library(lme4)
library(lmerTest)

######################### Read in Soil Data #####################################
Data <- read.csv(
  "Data/Site Prep Data - Soil.csv",
  stringsAsFactors = FALSE,
  check.names = TRUE
)

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Tables", showWarnings = FALSE, recursive = TRUE)

######################### Prepare Variables #####################################

# In this dataset, Sample identifies the experimental plot:
# 1-1 through 1-12 are in Block 1
# 2-1 through 2-12 are in Block 2

Data <- Data %>%
  mutate(
    Block = factor(sub("-.*", "", Sample)),
    Plot_ID = factor(Sample),
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "Fill Sand",
        "Preemergent",
        "Soil Inversion"
      )
    ),
    Time = factor(
      Time,
      levels = c("Pre", "Post")
    )
  )

# Confirm the design structure
print(table(Data$Block, Data$Treatment, Data$Time))
print(table(Data$Plot_ID, Data$Time))

if (nlevels(Data$Block) != 2) {
  stop("The dataset should contain exactly two blocks.")
}

if (nlevels(Data$Plot_ID) != 24) {
  stop("The dataset should contain exactly 24 unique plots.")
}

if (any(table(Data$Plot_ID, Data$Time) != 1)) {
  stop("Each plot should have one Pre and one Post observation.")
}

######################### Soil Properties #######################################
soil_properties <- c(
  "pH", "P", "K", "Ca", "Mg",
  "S", "Cu", "Mn", "Zn", "B"
)

missing_properties <- setdiff(
  soil_properties,
  names(Data)
)

if (length(missing_properties) > 0) {
  stop(
    "Missing soil-property columns: ",
    paste(missing_properties, collapse = ", ")
  )
}

######################### Treatment Colors ######################################
treat_colors <- c(
  "Control" = "yellow",
  "Fill Sand" = "orange",
  "Preemergent" = "#66CC00",
  "Soil Inversion" = "#CC66CC"
)

######################### Analysis and Plot Function ############################
plot_soil_prop <- function(df, prop) {
  
  cat("\n==================================================\n")
  cat("Soil property:", prop, "\n")
  cat("==================================================\n")
  
  # Block is a fixed design effect.
  # Plot_ID is a random intercept accounting for repeated Pre/Post
  # measurements from the same experimental plot.
  formula <- as.formula(
    paste0(
      prop,
      " ~ Block + Treatment * Time + (1 | Plot_ID)"
    )
  )
  
  model <- lmer(
    formula,
    data = df,
    REML = TRUE
  )
  
  # Type III tests with Satterthwaite denominator degrees of freedom.
  # Type III is appropriate here because the model contains an interaction.
  aov_res <- anova(
    model,
    type = 3,
    ddf = "Satterthwaite"
  )
  
  print(aov_res)
  
  ######################## Assumption Checks ####################################
  
  residuals_df <- data.frame(
    Fitted = fitted(model),
    Residual = residuals(model)
  )
  
  shapiro_res <- shapiro.test(
    residuals(model)
  )
  
  cat("\nShapiro-Wilk test of model residuals:\n")
  print(shapiro_res)
  
  ######################## Estimated Marginal Means #############################
  
  emm <- emmeans(
    model,
    ~ Treatment * Time
  )
  
  # Tukey-adjusted comparisons among all eight Treatment x Time combinations.
  pairwise_res <- pairs(
    emm,
    adjust = "tukey"
  )
  
  cld_res <- multcomp::cld(
    emm,
    adjust = "tukey",
    Letters = letters,
    alpha = 0.05
  ) %>%
    as.data.frame() %>%
    mutate(
      .group = trimws(.group)
    )
  
  print(emm)
  print(pairwise_res)
  
  ######################## Label Positions ######################################
  
  observed_range <- range(
    df[[prop]],
    na.rm = TRUE
  )
  
  label_offset <- diff(observed_range) * 0.08
  
  if (label_offset == 0) {
    label_offset <- max(abs(observed_range), 1) * 0.08
  }
  
  label_data <- df %>%
    group_by(
      Treatment,
      Time
    ) %>%
    summarise(
      y_pos = max(
        .data[[prop]],
        na.rm = TRUE
      ) + label_offset,
      .groups = "drop"
    ) %>%
    left_join(
      cld_res %>%
        dplyr::select(
          Treatment,
          Time,
          .group
        ),
      by = c(
        "Treatment",
        "Time"
      )
    )
  
  ######################## Plot #################################################
  
  p <- ggplot(
    df,
    aes(
      x = Time,
      y = .data[[prop]],
      fill = Treatment
    )
  ) +
    geom_boxplot(
      alpha = 0.5,
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    geom_point(
      color = "black",
      size = 2.5,
      alpha = 0.75,
      position = position_jitterdodge(
        jitter.width = 0.12,
        dodge.width = 0.8
      )
    ) +
    geom_text(
      data = label_data,
      aes(
        x = Time,
        y = y_pos,
        label = .group,
        group = Treatment
      ),
      position = position_dodge(width = 0.8),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = treat_colors
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.05, 0.18)
      )
    ) +
    labs(
      x = "",
      y = prop,
      title = prop
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        colour = "black"
      ),
      text = element_text(size = 16),
      axis.title.x = element_text(
        size = 15,
        face = "bold",
        colour = "black"
      ),
      axis.title.y = element_text(
        size = 15,
        face = "bold",
        colour = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        face = "bold",
        color = "black"
      ),
      legend.position = "bottom"
    )
  
  ######################## Return Results #######################################
  
  return(
    list(
      Property = prop,
      Model = model,
      ANOVA = aov_res,
      Shapiro = shapiro_res,
      EMM = emm,
      Pairwise = pairwise_res,
      CLD = cld_res,
      Plot = p
    )
  )
}

######################### Run All Soil Analyses #################################
soil_results <- lapply(
  soil_properties,
  function(prop) {
    plot_soil_prop(
      Data,
      prop
    )
  }
)

names(soil_results) <- soil_properties

######################### Extract Plots #########################################
soil_plots <- lapply(
  soil_results,
  function(x) {
    x$Plot
  }
)

######################### Arrange Plot Panel ####################################
soil_panel <- ggpubr::ggarrange(
  plotlist = soil_plots,
  ncol = 2,
  nrow = 5,
  common.legend = TRUE,
  legend = "bottom"
)

print(soil_panel)

ggsave(
  "Figures/Soil_Properties_Panel_with_Letters.tiff",
  plot = soil_panel,
  width = 12,
  height = 18,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

######################### Export ANOVA Results ##################################
anova_tables <- lapply(
  soil_results,
  function(x) {
    as.data.frame(x$ANOVA) %>%
      rownames_to_column("Effect") %>%
      mutate(
        Property = x$Property,
        .before = 1
      )
  }
)

anova_table_all <- bind_rows(
  anova_tables
)

write.csv(
  anova_table_all,
  "Tables/Soil_Properties_Repeated_Measures_ANOVA.csv",
  row.names = FALSE
)

######################### Export Estimated Marginal Means #######################
emm_tables <- lapply(
  soil_results,
  function(x) {
    as.data.frame(x$EMM) %>%
      mutate(
        Property = x$Property,
        .before = 1
      )
  }
)

emm_table_all <- bind_rows(
  emm_tables
)

write.csv(
  emm_table_all,
  "Tables/Soil_Properties_Estimated_Marginal_Means.csv",
  row.names = FALSE
)

######################### Export Tukey Comparisons ##############################
pairwise_tables <- lapply(
  soil_results,
  function(x) {
    as.data.frame(x$Pairwise) %>%
      mutate(
        Property = x$Property,
        .before = 1
      )
  }
)

pairwise_table_all <- bind_rows(
  pairwise_tables
)

write.csv(
  pairwise_table_all,
  "Tables/Soil_Properties_Tukey_Comparisons.csv",
  row.names = FALSE
)

######################### Export Compact Letter Groups ##########################
cld_tables <- lapply(
  soil_results,
  function(x) {
    x$CLD %>%
      mutate(
        Property = x$Property,
        .before = 1
      )
  }
)

cld_table_all <- bind_rows(
  cld_tables
)

write.csv(
  cld_table_all,
  "Tables/Soil_Properties_Compact_Letter_Groups.csv",
  row.names = FALSE
)

######################### Model Diagnostics Plots ###############################
pdf(
  "Figures/Soil_Properties_Model_Diagnostics.pdf",
  width = 8,
  height = 8
)

for (prop in soil_properties) {
  
  model <- soil_results[[prop]]$Model
  
  par(mfrow = c(1, 2))
  
  plot(
    fitted(model),
    residuals(model),
    xlab = "Fitted values",
    ylab = "Residuals",
    main = paste(prop, "- Residuals vs fitted")
  )
  
  abline(
    h = 0,
    lty = 2
  )
  
  qqnorm(
    residuals(model),
    main = paste(prop, "- Normal Q-Q")
  )
  
  qqline(
    residuals(model)
  )
}

dev.off()