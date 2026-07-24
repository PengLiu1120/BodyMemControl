# TacMem Recall phase encoded accuracy ANOVA & figures
# ..............................................................................
# written by P. Liu
# Email: peng.liu@uni-tuebingen.de
# last updated on 23rd July 2026 by P. Liu
# ..............................................................................
# ..............................................................................load packages

library(tidyverse)
library(rstatix)
library(lme4)
library(lmerTest)
library(emmeans)
library(effectsize)
library(ggdist)

stopifnot(requireNamespace("Hmisc", quietly = TRUE))
options(contrasts = c("contr.sum", "contr.poly"))

# ..............................................................................specify directory
setwd("/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/behavioural/exp3")
data_dir <- "/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/behavioural/exp3"
fig_dir  <- file.path(data_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ..............................................................................read data
condition_labels <- c(T = "Enhancement", B = "Baseline", NT = "Suppression")
bodypart_labels  <- c(left_hand  = "left hand",  right_hand = "right hand",
                      left_foot  = "left foot",  right_foot = "right foot")
stimulus_levels  <- c("rough", "neutral", "smooth")

accuracy_raw <- read_csv(file.path(data_dir, "data_learned.csv"), show_col_types = FALSE)

required_cols <- c("subject", "condition", "bodypart", "stimulus", "accuracy")
stopifnot(all(required_cols %in% names(accuracy_raw)))

accuracy_trial <- accuracy_raw %>%
  mutate(
    subject   = factor(subject),
    condition = factor(unname(condition_labels[condition]),
                       levels = unname(condition_labels)),
    bodypart  = factor(unname(bodypart_labels[bodypart]),
                       levels = unname(bodypart_labels)),
    stimulus  = factor(stimulus, levels = stimulus_levels)
  )

if (anyNA(select(accuracy_trial, condition, bodypart, stimulus))) {
  stop("Unrecognised level in condition, bodypart or stimulus; check the label maps.")
}

# ..............................................................................data aggregation
aggregate_accuracy <- function(data, ...) {
  data %>%
    group_by(subject, ...) %>%
    summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop")
}

acc_condition <- aggregate_accuracy(accuracy_trial, condition)
acc_cell      <- aggregate_accuracy(accuracy_trial, bodypart, stimulus, condition)
acc_body_stim <- aggregate_accuracy(accuracy_trial, stimulus, bodypart)

# ..............................................................................one-way repeated measures ANOVA
res_anova_1way <- anova_test(acc_condition,
                             dv = accuracy, wid = subject, within = condition)
print(get_anova_table(res_anova_1way))

pwc_1way <- pairwise_t_test(acc_condition, accuracy ~ condition,
                            paired = TRUE, p.adjust.method = "bonferroni")
print(pwc_1way)
print(rstatix::cohens_d(acc_condition, accuracy ~ condition, paired = TRUE))

# ..............................................................................three-way Linear Mixed Model
model_3way <- lmer(accuracy ~ bodypart * stimulus * condition + (1 | subject),
                   data = acc_cell)

res_anova_3way <- anova(model_3way)          # type III via lmerTest
print(res_anova_3way)
print(effectsize::eta_squared(model_3way, partial = TRUE))

emm_condition <- emmeans(model_3way, ~ condition)
print(pairs(emm_condition, adjust = "bonferroni"))
print(eff_size(emm_condition,
               sigma = sigma(model_3way),
               edf   = df.residual(model_3way)))

emm_body_by_stim <- emmeans(model_3way, ~ bodypart | stimulus)
print(pairs(emm_body_by_stim, adjust = "bonferroni"))

# ..............................................................................raincloud plots
condition_colours <- c("Enhancement" = "#2166ac",
                       "Baseline"    = "#878745",
                       "Suppression" = "#d6604d")

stimulus_colours  <- c("rough"   = "#e44a33",
                       "neutral" = "#becac3",
                       "smooth"  = "#7dc0f5")

theme_tacmem <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid      = element_blank(),
      axis.line       = element_line(colour = "black", linewidth = 1.2),
      axis.text       = element_text(size = base_size, colour = "black"),
      axis.title      = element_text(size = base_size),
      strip.text      = element_text(size = base_size, colour = "black"),
      legend.position = "none"
    )
}

raincloud_plot <- function(data, x, fill, colours,
                           y             = "accuracy",
                           subject_lines = FALSE,
                           mean_path     = FALSE,
                           facet         = NULL,
                           violin_width  = 0.5,
                           violin_just   = -0.2,
                           x_expand      = expansion(add = c(0.5, 0.8))) {
  
  pos_jitter <- position_jitter(width = 0.05, seed = 42)
  pos_nudge  <- position_nudge(x = -0.15)
  
  p <- ggplot(data, aes(x      = .data[[x]], y      = .data[[y]],
                        fill   = .data[[fill]], colour = .data[[fill]])) +
    ggdist::stat_halfeye(adjust        = 0.5,
                         width         = violin_width,
                         justification = violin_just,
                         .width        = 0,
                         point_colour  = NA,
                         alpha         = 0.55)
  
  if (subject_lines) {
    p <- p + geom_line(aes(group = subject), position = pos_jitter,
                       colour = "grey60", alpha = 0.3, linewidth = 0.4)
  }
  
  p <- p +
    geom_point(position = pos_jitter, size = 2, alpha = 0.65) +
    stat_summary(fun.data  = mean_cl_normal, geom = "pointrange",
                 colour    = "black", size = 0.5, linewidth = 1,
                 position  = pos_nudge)
  
  if (mean_path) {
    group_means <- data %>%
      group_by(.data[[x]]) %>%
      summarise(mean_y = mean(.data[[y]], na.rm = TRUE), .groups = "drop")
    p <- p + geom_line(data = group_means,
                       aes(x = .data[[x]], y = mean_y, group = 1),
                       colour = "black", linewidth = 1.2,
                       position = pos_nudge, inherit.aes = FALSE)
  }
  
  p <- p +
    scale_fill_manual(values   = colours) +
    scale_colour_manual(values = colours) +
    scale_x_discrete(expand    = x_expand) +
    scale_y_continuous(breaks  = seq(0, 1, by = 0.2),
                       labels  = scales::label_percent()) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = NULL) +
    theme_tacmem()
  
  if (!is.null(facet)) p <- p + facet_wrap(facet, ncol = 3)
  
  p
}

save_fig <- function(plot, filename, width, height) {
  ggsave(file.path(fig_dir, filename), plot = plot, units = "in",
         width = width, height = height, dpi = 600, compression = "lzw")
  message("Saved: ", file.path(fig_dir, filename))
}

# ..............................................................................save figures
p_condition <- raincloud_plot(acc_condition,
                              x             = "condition",
                              fill          = "condition",
                              colours       = condition_colours,
                              subject_lines = TRUE,
                              mean_path     = TRUE)

save_fig(p_condition, "tacmem_control_recall_condition_learned.tiff",
         width = 8, height = 6)

p_body_stim <- raincloud_plot(acc_body_stim,
                              x            = "bodypart",
                              fill         = "stimulus",
                              colours      = stimulus_colours,
                              facet        = ~ stimulus,
                              violin_width = 0.35,
                              violin_just  = -0.3,
                              x_expand     = expansion(add = c(0.5, 1.2)))

save_fig(p_body_stim, "tacmem_control_recall_bodypart_stimulus_learned.tiff",
         width = 18, height = 6)