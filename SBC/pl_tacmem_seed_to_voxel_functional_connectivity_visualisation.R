# TacMem Seed-to-Voxel Based Connectivity value visualisation
# Bar plots per contrast + raincloud plots per seed
# ..............................................................................
# written by P.Liu
# Email: peng.liu@uni-tuebingen.de
# last updated on 23rd July 2026 by P. Liu
# ..............................................................................
# ..............................................................................load packages
library(tidyverse)
library(ggdist)
library(RColorBrewer)
library(patchwork)

# ..............................................................................load data
data_dir <- "/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/MRI/stats"
out_dir  <- file.path(data_dir, "SBC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv(file.path(data_dir, "all_seeds_connectivity_values.csv"),
                 header = TRUE)

# ..............................................................................define colour scheme
colour_map <- c(
  "Think(1).NoThink(-1)_increased"                        = "#1d4e89",
  "Think(1).NoThink(-1)_decreased"                        = "#90bae0",
  "Think(-1).NoThink(1)_increased"                        = "#a50026",
  "Think(-1).NoThink(1)_decreased"                        = "#e8a0a0",
  "Think_Recall(1).Baseline_Recall(-1)_increased"         = "#1a7a3c",
  "Think_Recall(1).Baseline_Recall(-1)_decreased"         = "#90d4a8"
)

# ..............................................................................create colour key per row
data <- data %>%
  mutate(colour_key = paste0(contrast, "_", direction))

# ..............................................................................define readable labels
seed_labels <- c(
  "TNT_ROI_l_IPL"        = "L IPL",
  "TNT_ROI_l_dlPFC"      = "L dlPFC",
  "TNT_ROI_l_precuneus"  = "L PCu",
  "TNT_ROI_r_precuneus"  = "R PCu",
  "TNT_ROI_r_insula"     = "R Insula",
  "NTT_ROI_r_vlPFC"      = "R vlPFC",
  "NTT_ROI_r_dlPFC"      = "R dlPFC",
  "TB_ROI_l_IPL"         = "L IPL",
  "TB_ROI_r_IPL"         = "R IPL"
)

target_labels <- c(
  "L_LateralOccipital_PCu"  = "L Lat Occ/PCu",
  "R_IFG_BA44"              = "R IFG",
  "R_OccipitalPole"         = "R Occ Pole",
  "R_MFG_BA8"               = "R MFG",
  "L_IFG_OFC_Insula"        = "L IFG/OFC",
  "L_IFG_FrontalPole"       = "L IFG/FP",
  "L_IFO_Precentral"        = "L IFO/PreCG",
  "R_MFG_Premotor"          = "R MFG/PM",
  "L_IFG_MFG"               = "L IFG/MFG",
  "R_SMA_SFG"               = "R SMA",
  "L_Precentral_IFG"        = "L PreCG/IFG",
  "B_Paracingulate"         = "B Paracingulate",
  "R_FrontalPole"           = "R Frontal Pole",
  "R_SPL"                   = "R SPL",
  "L_Supramarginal"         = "L SupMar",
  "R_InferiorTemporal"      = "R ITG",
  "R_LateralOccipital"      = "R Lat Occ",
  "L_LateralOccipital"      = "L Lat Occ",
  "R_Fusiform"              = "R Fusiform"
)

data <- data %>%
  mutate(
    seed_label         = recode(seed,   !!!seed_labels),
    target_label_short = recode(target, !!!target_labels),
    seed_target        = paste0(seed_label, " → ", target_label_short)
  )

# ..............................................................................one sample t-test against zero
ttest_results <- data %>%
  group_by(contrast, seed, target, seed_target, direction, colour_key) %>%
  summarise(
    mean_conn          = mean(connectivity, na.rm = TRUE),
    se_conn            = sd(connectivity, na.rm = TRUE) / sqrt(n()),
    t_stat             = t.test(connectivity, mu = 0)$statistic,
    p_value            = t.test(connectivity, mu = 0)$p.value,
    target_label_short = first(recode(target, !!!target_labels)),
    .groups = 'drop'
  )

# ..............................................................................plot 1 — bar plots per contrast
plot_bar_contrast <- function(contrast_name, title_label) {
  
  d <- ttest_results %>%
    filter(contrast == contrast_name) %>%
    arrange(seed, direction)
  
  if (nrow(d) == 0) return(NULL)
  
  ggplot(d, aes(x = seed_target, y = mean_conn, fill = colour_key)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 1.0) +
    geom_bar(stat = "identity", width = 0.8, colour = "white", linewidth = 0.5) +
    geom_errorbar(
      aes(ymin = mean_conn - se_conn, ymax = mean_conn + se_conn),
      width = 0.5, linewidth = 1.0, colour = "black"
    ) +
    scale_fill_manual(values = colour_map) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 14)) +
    labs(
      title = NULL,
      x     = NULL,
      y     = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.line          = element_line(colour = "black", linewidth = 2),
      axis.text.x        = element_text(size = 16, angle = 35, hjust = 1),
      axis.text.y        = element_text(size = 16),
      legend.position    = "none"
    )
}

p_bar_TNT    <- plot_bar_contrast("Think(1).NoThink(-1)",                "Think > NoThink")
p_bar_NTT    <- plot_bar_contrast("Think(-1).NoThink(1)",                "NoThink > Think")
p_bar_recall <- plot_bar_contrast("Think_Recall(1).Baseline_Recall(-1)", "Think > Baseline (Recall)")

# ..............................................................................plot 2 — raincloud plots per seed
plot_raincloud_seed <- function(seed_name, title_label, contrast_name) {
  
  d <- data %>%
    filter(seed == seed_name, contrast == contrast_name)
  
  if (nrow(d) == 0) return(NULL)
  
  ggplot(d, aes(x = target_label_short, y = connectivity,
                fill = colour_key, colour = colour_key)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 1.0) +
    ggdist::stat_halfeye(
      adjust        = 0.6,
      width         = 0.5,
      justification = -0.25,
      .width        = 0,
      point_colour  = NA,
      alpha         = 0.55
    ) +
    geom_point(
      position = position_jitter(width = 0.06, seed = 42),
      size     = 2.0,
      alpha    = 0.65
    ) +
    stat_summary(
      fun.data  = mean_cl_normal,
      geom      = "pointrange",
      colour    = "black",
      linewidth = 1.0,
      size      = 1.0,
      position  = position_nudge(x = -0.18)
    ) +
    scale_fill_manual(values   = colour_map) +
    scale_colour_manual(values = colour_map) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) +
    labs(
      title = NULL,
      x     = NULL,
      y     = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.line          = element_line(colour = "black", linewidth = 2.0),
      axis.text.x        = element_text(size = 16),
      axis.text.y        = element_text(size = 16),
      axis.title.y       = element_text(size = 12),
      plot.title         = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position    = "none"
    )
}

# ..............................................................................generate raincloud plots
p_rc_l_IPL_TNT    <- plot_raincloud_seed("TNT_ROI_l_IPL",       "L IPL  |  Think > NoThink",    "Think(1).NoThink(-1)")
p_rc_l_dlPFC      <- plot_raincloud_seed("TNT_ROI_l_dlPFC",     "L dlPFC  |  Think > NoThink",  "Think(1).NoThink(-1)")
p_rc_l_PCu        <- plot_raincloud_seed("TNT_ROI_l_precuneus", "L PCu  |  Think > NoThink",    "Think(1).NoThink(-1)")
p_rc_r_PCu        <- plot_raincloud_seed("TNT_ROI_r_precuneus", "R PCu  |  Think > NoThink",    "Think(1).NoThink(-1)")
p_rc_r_insula     <- plot_raincloud_seed("TNT_ROI_r_insula",    "R Insula  |  Think > NoThink", "Think(1).NoThink(-1)")
p_rc_r_vlPFC      <- plot_raincloud_seed("NTT_ROI_r_vlPFC",     "R vlPFC  |  NoThink > Think",  "Think(-1).NoThink(1)")
p_rc_r_dlPFC      <- plot_raincloud_seed("NTT_ROI_r_dlPFC",     "R dlPFC  |  NoThink > Think",  "Think(-1).NoThink(1)")
p_rc_l_IPL_recall <- plot_raincloud_seed("TB_ROI_l_IPL",        "L IPL  |  Think > Baseline",   "Think_Recall(1).Baseline_Recall(-1)")
p_rc_r_IPL_recall <- plot_raincloud_seed("TB_ROI_r_IPL",        "R IPL  |  Think > Baseline",   "Think_Recall(1).Baseline_Recall(-1)")

# ..............................................................................save bar plots
ggsave(file.path(out_dir, "tacmem_connectivity_bar_TNT_TvsNT.tiff"),
       plot = p_bar_TNT,    units = "in", width = 10, height = 5, dpi = 600)
ggsave(file.path(out_dir, "tacmem_connectivity_bar_TNT_NTvsT.tiff"),
       plot = p_bar_NTT,    units = "in", width = 8,  height = 5, dpi = 600)
ggsave(file.path(out_dir, "tacmem_connectivity_bar_recall_TvsB.tiff"),
       plot = p_bar_recall, units = "in", width = 8,  height = 5, dpi = 600)
cat("Saved bar plots\n")

# ..............................................................................save raincloud plots
rc_plots <- list(
  "tacmem_connectivity_rc_l_IPL_TNT.tiff"    = p_rc_l_IPL_TNT,
  "tacmem_connectivity_rc_l_dlPFC_TNT.tiff"  = p_rc_l_dlPFC,
  "tacmem_connectivity_rc_l_PCu_TNT.tiff"    = p_rc_l_PCu,
  "tacmem_connectivity_rc_r_PCu_TNT.tiff"    = p_rc_r_PCu,
  "tacmem_connectivity_rc_r_insula_TNT.tiff" = p_rc_r_insula,
  "tacmem_connectivity_rc_r_vlPFC_NTT.tiff"  = p_rc_r_vlPFC,
  "tacmem_connectivity_rc_r_dlPFC_NTT.tiff"  = p_rc_r_dlPFC,
  "tacmem_connectivity_rc_l_IPL_recall.tiff" = p_rc_l_IPL_recall,
  "tacmem_connectivity_rc_r_IPL_recall.tiff" = p_rc_r_IPL_recall
)

for (rc_name in names(rc_plots)) {
  if (!is.null(rc_plots[[rc_name]])) {
    ggsave(file.path(out_dir, rc_name), plot = rc_plots[[rc_name]],
           units = "in", width = 6, height = 5, dpi = 600)
    cat("Saved:", rc_name, "\n")
  } else {
    cat("Skipped (no data):", rc_name, "\n")
  }
}

cat("\nAll connectivity plots saved to:", out_dir, "\n")