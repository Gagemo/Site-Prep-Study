################################################################################
################################################################################
#########################  FWF - Site Prep Study  ##############################
#########################        Soil Chemistry   ##############################
#########################  University of Florida  ##############################
#########################      Gage LaPierre      ##############################
#########################      2024 - 2025        ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014")

#########################      Installs Packages   ##############################
list.of.packages <- c("tidyverse", "car", "ggpubr", "rstatix",
                      "emmeans", "multcomp", "multcompView")
new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################      Loads Packages      ##############################
library(tidyverse)
library(car)
library(ggpubr)
library(rstatix)
library(emmeans)
library(multcompView)
library(multcomp)


##########################   Read in Soil Data   ##############################
Data <- read.csv("Data/Site Prep Data - Soil.csv")

# Ensure factors are set correctly
Data$Treatment <- factor(Data$Treatment,
                         levels = c("Control", "Fill Sand", "Preemergent", "Soil Inversion"))
Data$Time <- factor(Data$Time, levels = c("Pre", "Post"))

##########################   Function for ANOVA + Letters   #####################
plot_soil_prop <- function(df, prop) {
  
  # Build model with Treatment, Time, and their interaction
  formula <- as.formula(paste(prop, "~ Treatment * Time"))
  model <- lm(formula, data = df)
  
  # Type II ANOVA
  aov_res <- Anova(model, type = 2)
  print(aov_res)
  
  # Estimated marginal means by Treatment x Time
  emm <- emmeans(model, ~ Treatment * Time)
  
  # Compact letter display for interaction
  cld_res <- cld(emm, adjust = "tukey") %>%
    as.data.frame()
  
  # Position letters above max values for each group
  # Use dplyr::select to ensure the correct function is called
  label_data <- df %>%
    group_by(Treatment, Time) %>%
    summarise(y_pos = max(.data[[prop]], na.rm = TRUE) * 1.05,
              .groups = "drop") %>%
    left_join(cld_res %>% dplyr::select(Treatment, Time, .group),
              by = c("Treatment","Time"))
  
  
  # Custom colors for treatments
  treat_colors <- c("Control" = "yellow",
                    "Fill Sand" = "orange",
                    "Preemergent" = "#66CC00",
                    "Soil Inversion" = "#CC66CC")
  
  # Plot: split by Time on x-axis
  p <- ggplot(df, aes(x = Time, y = .data[[prop]], fill = Treatment)) +
    geom_boxplot(aes(y = .data[[prop]]), alpha = 0.5, outlier.shape = NA,
                 position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color = Treatment, y = .data[[prop]]),
                size = 3, alpha = 0.7,
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
    geom_text(data = label_data,
              aes(x = Time, y = y_pos, label = .group, group = Treatment),
              position = position_dodge(width = 0.8),
              inherit.aes = FALSE,
              size = 5, fontface = "bold") +
    scale_fill_manual(values = treat_colors) +
    scale_color_manual(values = treat_colors) +
    labs(x = "", y = paste(prop, "value"),
         title = paste("", prop, "")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", colour = "black"),
          text = element_text(size = 16),
          axis.title.x = element_text(size = 15, face = "bold", colour = "black"),
          axis.title.y = element_text(size = 15, face = "bold", colour = "black"),
          axis.text.x = element_text(size = 13, face = "bold", color = "black"),
          axis.text.y = element_text(size = 13, face = "bold", color = "black"),
          legend.position = "bottom")
  
  return(p)
}

##########################   Generate Plots for Key Props   #####################
soil_plots <- lapply(c("pH", "P", 
                       "K", "Ca", 
                       "Mg", "S", 
                       "Cu", "Mn", 
                       "Zn", "B"),
                     function(prop) plot_soil_prop(Data, prop))

soil_panel <- ggpubr::ggarrange(plotlist = soil_plots, ncol = 2, nrow = 5,
                                common.legend = TRUE, legend = "bottom")
ggsave("Figures/Soil_Properties_Panel_with_Letters.jpg", soil_panel,
       width = 12, height = 10)
