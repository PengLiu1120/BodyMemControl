# TacMem Recall Phase Accuracy
# ..............................................................................
# written by P.Liu
# Email: peng.liu@uni-tuebingen.de
# last updated on 09th June 2026 by P.Liu
# ..............................................................................

# ..............................................................................load packages
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggdist)
library(patchwork)

# ..............................................................................read data
setwd("/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/behavioural/exp3")

Accuracy <- read.table("data2.csv", header = TRUE, sep = ',')

Accuracy <- Accuracy %>%
  mutate(condition = case_when(
    condition == "T"  ~ "Enhancement",
    condition == "B"  ~ "Baseline",
    condition == "NT" ~ "Suppression"
  )) %>%
  mutate(condition = factor(condition, levels = c("Enhancement", "Baseline", "Suppression")))

# ..............................................................................data aggregation
Accuracy_1way_agg <- Accuracy %>%
  group_by(subject, condition) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop")

# ..............................................................................one-way repeated measures ANOVA
res.anova_1way <- Accuracy_1way_agg %>%
  anova_test(dv = accuracy, wid = subject, within = condition)

get_anova_table(res.anova_1way)

# ..............................................................................multiple comparisons
pwc_1way <- Accuracy_1way_agg %>%
  pairwise_t_test(
    accuracy ~ condition,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )

# ..............................................................................define colours per condition
condition_colours <- c(
  "Enhancement" = "#2166ac",
  "Baseline"    = "#878745",
  "Suppression" = "#d6604d"
)

# ..............................................................................helper function for raincloud plot
raincloud_plot <- function(data, x_var, y_var, fill_var, colour_var,
                           colour_map, x_label = NULL, y_label = NULL,
                           facet_formula = NULL) {
  
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]],
                        fill = .data[[fill_var]], colour = .data[[colour_var]])) +
    
    # ..........................................................................half violin
    ggdist::stat_halfeye(
      adjust        = 0.5,
      width         = 0.5,
      justification = -0.2,
      .width        = 0,
      point_colour  = NA,
      alpha         = 0.55
    ) +
    
    # ..........................................................................individual points
    geom_point(
      position = position_jitter(width = 0.05, seed = 42),
      size     = 2.0,
      alpha    = 0.65
    ) +
    
    # ..........................................................................mean and 95% CI
    stat_summary(
      fun.data  = mean_cl_normal,
      geom      = "pointrange",
      colour    = "black",
      size      = 0.5,
      linewidth = 1.0,
      position  = position_nudge(x = -0.15)
    ) +
    
    scale_fill_manual(values   = colour_map) +
    scale_colour_manual(values = colour_map) +
    
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      labels = scales::label_percent()
    ) +
    
    coord_cartesian(ylim = c(0, 1)) +
    
    labs(x = x_label, y = y_label, title = NULL) +
    
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = 1.2),
      axis.text.x      = element_text(size = 16, colour = "black"),
      axis.text.y      = element_text(size = 16, colour = "black"),
      axis.title.y     = element_text(size = 16),
      legend.position  = "none"
    )
  
  if (!is.null(facet_formula)) {
    p <- p + facet_wrap(facet_formula, scales = "free_x")
  }
  
  return(p)
}

# ..............................................................................PLOT 1 — condition raincloud
# pre-compute means for mean line
means_condition <- Accuracy_1way_agg %>%
  group_by(condition) %>%
  summarise(mean_acc = mean(accuracy, na.rm = TRUE), .groups = "drop")

p1 <- ggplot(Accuracy_1way_agg, aes(x = condition, y = accuracy,
                                    fill = condition, colour = condition)) +
  ggdist::stat_halfeye(
    adjust        = 0.5,
    width         = 0.5,
    justification = -0.2,
    .width        = 0,
    point_colour  = NA,
    alpha         = 0.55
  ) +
  geom_point(
    position = position_jitter(width = 0.05, seed = 42),
    size     = 2.0,
    alpha    = 0.65
  ) +
  geom_line(
    aes(group = subject),
    position  = position_jitter(width = 0.05, seed = 42),
    colour    = "grey60",
    alpha     = 0.3,
    linewidth = 0.4
  ) +
  stat_summary(
    fun.data  = mean_cl_normal,
    geom      = "pointrange",
    colour    = "black",
    size      = 0.5,
    linewidth = 1.0,
    position  = position_nudge(x = -0.15)
  ) +
  geom_line(
    data        = means_condition,
    aes(x       = condition, y = mean_acc, group = 1),
    colour      = "black",
    linewidth   = 1.2,
    linetype    = "solid",
    position    = position_nudge(x = -0.15),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values   = condition_colours) +
  scale_colour_manual(values = condition_colours) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::label_percent()) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 1.2),
    axis.text.x      = element_text(size = 16, colour = "black"),
    axis.text.y      = element_text(size = 16, colour = "black"),
    legend.position  = "none"
  )

ggsave("TNT_recall_condition.tiff", plot = p1, units = "in", width = 8, height = 6, dpi = 600)
cat("Saved: TNT_recall_condition.tiff\n")

# ..............................................................................PLOT 2 — body part x condition raincloud
bodypart_agg <- Accuracy %>%
  group_by(subject, condition, bodypart) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop")

bodypart_colours <- c("hand" = "#91bbd3", "foot" = "#d26d6f")

ggplot(bodypart_agg, aes(x = condition, y = accuracy,
                         fill = bodypart, colour = bodypart)) +
  ggdist::stat_halfeye(
    adjust        = 0.5,
    width         = 0.35,
    justification = -0.3,
    .width        = 0,
    point_colour  = NA,
    alpha         = 0.55
  ) +
  geom_point(
    position = position_jitter(width = 0.05, seed = 42),
    size     = 2.0,
    alpha    = 0.65
  ) +
  geom_line(
    aes(group = subject),
    position  = position_jitter(width = 0.05, seed = 42),
    colour    = "grey60",
    alpha     = 0.3,
    linewidth = 0.4
  ) +
  stat_summary(
    fun.data  = mean_cl_normal,
    geom      = "pointrange",
    colour    = "black",
    size      = 0.5,
    linewidth = 1.0,
    position  = position_nudge(x = -0.15)
  ) +
  geom_line(
    data = bodypart_agg %>%
      group_by(condition, bodypart) %>%
      summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop"),
    aes(x = condition, y = accuracy, group = 1),
    colour    = "black",
    linewidth = 1.2,
    linetype  = "solid",
    position  = position_nudge(x = -0.15),
    inherit.aes = FALSE
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.2))) +
  facet_wrap(~ bodypart, ncol = 2) +
  scale_fill_manual(values   = bodypart_colours) +
  scale_colour_manual(values = bodypart_colours) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::label_percent()) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 1.2),
    axis.text.x      = element_text(size = 16, colour = "black"),
    axis.text.y      = element_text(size = 16, colour = "black"),
    legend.position  = "none"
  )

ggsave("TNT_recall_bodypart_condition.tiff", units = "in", width = 10, height = 6, dpi = 600)
cat("Saved: TNT_recall_bodypart_condition.tiff\n")

# ..............................................................................PLOT 3 — stimulus x condition raincloud
stimulus_agg <- Accuracy %>%
  group_by(subject, condition, stimulus) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
  mutate(stimulus = case_when(
    stimulus == "rough"   ~ "rough",
    stimulus == "neutral" ~ "neutral",
    stimulus == "smooth"  ~ "smooth"
  )) %>%
  mutate(stimulus = factor(stimulus, levels = c("rough", "neutral", "smooth")))

stimulus_colours <- c("rough" = "#e44a33", "neutral" = "#becac3", "smooth" = "#7dc0f5")

ggplot(stimulus_agg, aes(x = condition, y = accuracy,
                         fill = stimulus, colour = stimulus)) +
  ggdist::stat_halfeye(
    adjust        = 0.5,
    width         = 0.35,
    justification = -0.3,
    .width        = 0,
    point_colour  = NA,
    alpha         = 0.55
  ) +
  geom_point(
    position = position_jitter(width = 0.05, seed = 42),
    size     = 2.0,
    alpha    = 0.65
  ) +
  geom_line(
    aes(group = subject),
    position  = position_jitter(width = 0.05, seed = 42),
    colour    = "grey60",
    alpha     = 0.3,
    linewidth = 0.4
  ) +
  stat_summary(
    fun.data  = mean_cl_normal,
    geom      = "pointrange",
    colour    = "black",
    size      = 0.5,
    linewidth = 1.0,
    position  = position_nudge(x = -0.15)
  ) +
  geom_line(
    data = stimulus_agg %>%
      group_by(condition, stimulus) %>%
      summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop"),
    aes(x = condition, y = accuracy, group = 1),
    colour    = "black",
    linewidth = 1.2,
    linetype  = "solid",
    position  = position_nudge(x = -0.15),
    inherit.aes = FALSE
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.2))) +
  facet_wrap(~ stimulus, ncol = 3) +
  scale_fill_manual(values   = stimulus_colours) +
  scale_colour_manual(values = stimulus_colours) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::label_percent()) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 1.2),
    axis.text.x      = element_text(size = 16, colour = "black"),
    axis.text.y      = element_text(size = 16, colour = "black"),
    legend.position  = "none"
  )

ggsave("TNT_recall_stimulus_condition.tiff", units = "in", width = 18, height = 6, dpi = 600)
cat("Saved: TNT_recall_stimulus_condition.tiff\n")

# ..............................................................................PLOT 4 — side x bodypart x stimulus x condition raincloud
interaction_agg_full <- Accuracy %>%
  group_by(subject, condition, stimulus, bodypart, side) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    stimulus = factor(case_when(
      stimulus == "rough"   ~ "rough",
      stimulus == "neutral" ~ "neutral",
      stimulus == "smooth"  ~ "smooth"
    ), levels = c("rough", "neutral", "smooth")),
    bodypart = factor(case_when(
      bodypart == "hand" ~ "hand",
      bodypart == "foot" ~ "foot"
    )),
    side = factor(case_when(
      side == "left"  ~ "left",
      side == "right" ~ "right"
    ))
  )

ggplot(interaction_agg_full, aes(x = stimulus, y = accuracy,
                                 fill = stimulus, colour = stimulus)) +
  ggdist::stat_halfeye(
    adjust        = 0.5,
    width         = 0.35,
    justification = -0.3,
    .width        = 0,
    point_colour  = NA,
    alpha         = 0.55
  ) +
  geom_point(
    position = position_jitter(width = 0.05, seed = 42),
    size     = 1.5,
    alpha    = 0.65
  ) +
  stat_summary(
    fun.data  = mean_cl_normal,
    geom      = "pointrange",
    colour    = "black",
    size      = 0.4,
    linewidth = 0.8,
    position  = position_nudge(x = -0.15)
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.2))) +
  facet_grid(bodypart + side ~ condition) +
  scale_fill_manual(values   = stimulus_colours) +
  scale_colour_manual(values = stimulus_colours) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::label_percent()) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 1.2),
    axis.text.x      = element_text(size = 14, colour = "black"),
    axis.text.y      = element_text(size = 14, colour = "black"),
    legend.position  = "none"
  )

ggsave("TNT_recall_side_bodypart_stimulus_condition.tiff", units = "in", width = 16, height = 12, dpi = 600)
cat("Saved: TNT_recall_side_bodypart_stimulus_condition.tiff\n")

# ..............................................................................PLOT 5 — body part x stimulus x condition raincloud
interaction_agg_bp <- Accuracy %>%
  group_by(subject, condition, stimulus, bodypart) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    stimulus = factor(case_when(
      stimulus == "rough"   ~ "rough",
      stimulus == "neutral" ~ "neutral",
      stimulus == "smooth"  ~ "smooth"
    ), levels = c("rough", "neutral", "smooth")),
    bodypart = factor(case_when(
      bodypart == "hand" ~ "hand",
      bodypart == "foot" ~ "foot"
    ))
  )

ggplot(interaction_agg_bp, aes(x = stimulus, y = accuracy,
                               fill = stimulus, colour = stimulus)) +
  ggdist::stat_halfeye(
    adjust        = 0.5,
    width         = 0.35,
    justification = -0.3,
    .width        = 0,
    point_colour  = NA,
    alpha         = 0.55
  ) +
  geom_point(
    position = position_jitter(width = 0.05, seed = 42),
    size     = 1.5,
    alpha    = 0.65
  ) +
  stat_summary(
    fun.data  = mean_cl_normal,
    geom      = "pointrange",
    colour    = "black",
    size      = 0.4,
    linewidth = 0.8,
    position  = position_nudge(x = -0.15)
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.2))) +
  facet_grid(bodypart ~ condition) +
  scale_fill_manual(values   = stimulus_colours) +
  scale_colour_manual(values = stimulus_colours) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::label_percent()) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 1.2),
    axis.text.x      = element_text(size = 16, colour = "black"),
    axis.text.y      = element_text(size = 16, colour = "black"),
    legend.position  = "none"
  )

ggsave("TNT_recall_bodypart_stimulus_condition.tiff", units = "in", width = 14, height = 6, dpi = 600)
cat("Saved: TNT_recall_bodypart_stimulus_condition.tiff\n")

# ..............................................................................PLOT 6 — bodypart x stimulus raincloud
stimulus_bodypart_agg <- Accuracy %>%
  group_by(subject, stimulus, bodypart) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    stimulus = factor(case_when(
      stimulus == "rough"   ~ "rough",
      stimulus == "neutral" ~ "neutral",
      stimulus == "smooth"  ~ "smooth"
    ), levels = c("rough", "neutral", "smooth")),
    bodypart = factor(case_when(
      bodypart == "hand" ~ "hand",
      bodypart == "foot" ~ "foot"
    ))
  )

ggplot(stimulus_bodypart_agg, aes(x = bodypart, y = accuracy,
                                  fill = stimulus, colour = stimulus)) +
  ggdist::stat_halfeye(
    adjust        = 0.5,
    width         = 0.35,
    justification = -0.3,
    .width        = 0,
    point_colour  = NA,
    alpha         = 0.55
  ) +
  geom_point(
    position = position_jitter(width = 0.05, seed = 42),
    size     = 2.0,
    alpha    = 0.65
  ) +
  stat_summary(
    fun.data  = mean_cl_normal,
    geom      = "pointrange",
    colour    = "black",
    size      = 0.5,
    linewidth = 1.0,
    position  = position_nudge(x = -0.15)
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.2))) +
  facet_wrap(~ stimulus, ncol = 3) +
  scale_fill_manual(values   = stimulus_colours) +
  scale_colour_manual(values = stimulus_colours) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::label_percent()) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 1.2),
    axis.text.x      = element_text(size = 16, colour = "black"),
    axis.text.y      = element_text(size = 16, colour = "black"),
    legend.position  = "none"
  )

ggsave("TNT_recall_bodypart_stimulus.tiff", units = "in", width = 14, height = 6, dpi = 600)
cat("Saved: TNT_recall_bodypart_stimulus.tiff\n")

# ..............................................................................PLOT 7 — side x bodypart x condition raincloud
side_bodypart_agg <- Accuracy %>%
  group_by(subject, condition, side, bodypart) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    bodypart = factor(case_when(
      bodypart == "hand" ~ "hand",
      bodypart == "foot" ~ "foot"
    )),
    side = factor(case_when(
      side == "left"  ~ "left",
      side == "right" ~ "right"
    ))
  )

ggplot(side_bodypart_agg, aes(x = condition, y = accuracy,
                              fill = condition, colour = condition)) +
  ggdist::stat_halfeye(
    adjust        = 0.5,
    width         = 0.35,
    justification = -0.3,
    .width        = 0,
    point_colour  = NA,
    alpha         = 0.55
  ) +
  geom_point(
    position = position_jitter(width = 0.05, seed = 42),
    size     = 1.5,
    alpha    = 0.65
  ) +
  geom_line(
    aes(group = subject),
    position  = position_jitter(width = 0.05, seed = 42),
    colour    = "grey60",
    alpha     = 0.3,
    linewidth = 0.4
  ) +
  stat_summary(
    fun.data  = mean_cl_normal,
    geom      = "pointrange",
    colour    = "black",
    size      = 0.4,
    linewidth = 0.8,
    position  = position_nudge(x = -0.15)
  ) +
  geom_line(
    data = side_bodypart_agg %>%
      group_by(condition, side, bodypart) %>%
      summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop"),
    aes(x = condition, y = accuracy, group = 1),
    colour    = "black",
    linewidth = 1.0,
    linetype  = "solid",
    position  = position_nudge(x = -0.15),
    inherit.aes = FALSE
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.2))) +
  facet_grid(bodypart ~ side) +
  scale_fill_manual(values   = condition_colours) +
  scale_colour_manual(values = condition_colours) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = scales::label_percent()) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 1.2),
    axis.text.x      = element_text(size = 16, colour = "black"),
    axis.text.y      = element_text(size = 16, colour = "black"),
    legend.position  = "none"
  )

ggsave("TNT_recall_side_bodypart_condition.tiff", units = "in", width = 14, height = 6, dpi = 600)
cat("Saved: TNT_recall_side_bodypart_condition.tiff\n")

# ..............................................................................four-way Linear Mixed Model
Accuracy_4way_agg <- Accuracy %>%
  group_by(subject, side, bodypart, stimulus, condition) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop")

model_4way <- lmer(accuracy ~ side * bodypart * stimulus * condition + (1|subject),
                   data = Accuracy)
anova(model_4way)
res.anova_4way <- anova(model_4way)

posthoc_4way_condition           <- emmeans(model_4way, specs = pairwise ~ condition, adjust = "bonferroni")
posthoc_4way_bodypart_stimulus   <- emmeans(model_4way, pairwise ~ bodypart | stimulus, adjust = "bonferroni")
posthoc_4way_side_bodypart_condition <- emmeans(model_4way, pairwise ~ side | bodypart | condition, adjust = "bonferroni")

summary(posthoc_4way_condition$contrasts)
summary(posthoc_4way_bodypart_stimulus$contrasts)
summary(posthoc_4way_side_bodypart_condition$contrasts)

# ..............................................................................three-way Linear Mixed Model
Accuracy_3way_agg <- Accuracy %>%
  group_by(subject, bodypart, stimulus, condition) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop")

model_3way <- lmer(accuracy ~ bodypart * stimulus * condition + (1|subject),
                   data = Accuracy)
anova(model_3way)
res.anova_3way <- anova(model_3way)

posthoc_3way_condition          <- emmeans(model_3way, specs = pairwise ~ condition, adjust = "bonferroni")
posthoc_3way_stimulus           <- emmeans(model_3way, specs = pairwise ~ stimulus, adjust = "bonferroni")
posthoc_3way_bodypart_stimulus  <- emmeans(model_3way, pairwise ~ bodypart | stimulus, adjust = "bonferroni")
posthoc_3way_stimulus_condition <- emmeans(model_3way, pairwise ~ stimulus | condition, adjust = "bonferroni")

summary(posthoc_3way_condition$contrasts)
summary(posthoc_3way_stimulus$contrasts)
summary(posthoc_3way_bodypart_stimulus$contrasts)
summary(posthoc_3way_stimulus_condition$contrasts)