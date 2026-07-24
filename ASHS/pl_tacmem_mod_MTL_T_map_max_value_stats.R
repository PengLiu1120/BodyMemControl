# TacMem Modulation phase MTL T map value stats
# ..............................................................................
# written by P.Liu
# Email: peng.liu@uni-tuebingen.de
# last updated on 23rd July 2026 by P. Liu
# ..............................................................................
# ..............................................................................load packages
# ez has plyr and MASS in Depends, so library(ez) ATTACHES both and they mask
# dplyr::summarise, dplyr::mutate, dplyr::arrange, dplyr::rename and
# dplyr::select. tidyverse is therefore attached LAST so that dplyr wins, and
# the dplyr:: prefixes below guard against the order ever changing.
library(ez)
library(emmeans)
library(effectsize)
library(lme4)
library(tidyverse)

# ..............................................................................specify parameters
stats_dir <- "/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/MRI/stats"
out_dir   <- file.path(stats_dir, "MTL")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

subfield_order <- c("BA35", "BA36", "CA1", "CA2", "CA3",
                    "CS", "DG", "ERC", "PHC", "SUB", "Tail")

wong_palette <- c(
  "#E69F00", # orange        — BA35
  "#56B4E9", # sky blue      — BA36
  "#009E73", # bluish green  — CA1
  "#F0E442", # yellow        — CA2
  "#0072B2", # blue          — CA3
  "#D55E00", # vermillion    — CS
  "#CC79A7", # reddish purple— DG
  "#000000", # black         — ERC
  "#999999", # grey          — PHC
  "#44AA99", # teal          — SUB
  "#882255"  # dark wine     — Tail
)
subfield_colours <- setNames(wong_palette, subfield_order)

# Reference line positions (Cohen's d benchmarks)
REF_SMALL  <- 0.2
REF_MEDIUM <- 0.5
REF_LARGE  <- 0.8

hemisphere_labels <- c(left = "LH", right = "RH")

# ..............................................................................one-sample t test vs 0
# one t.test per subfield rather than three identical calls per row
compute_stats <- function(df) {
  df %>%
    dplyr::group_by(Subfield_Name) %>%
    dplyr::group_modify(~ {
      x  <- .x$Max_T[!is.na(.x$Max_T)]
      tt <- t.test(x, mu = 0)
      tibble(
        mean    = mean(x),
        sd      = sd(x),
        n       = length(x),                 # non-missing, so SEM is correct
        SEM     = sd(x) / sqrt(length(x)),
        t_stat  = unname(tt$statistic),
        df      = unname(tt$parameter),
        p_value = tt$p.value
      )
    }) %>%
    dplyr::ungroup() %>%
    # NOTE: FDR is applied across the 11 subfields WITHIN one condition map and
    # hemisphere. If the family is all 44 tests, move this to the combined table
    # in the export section instead.
    dplyr::mutate(p_fdr = p.adjust(p_value, method = "fdr"))
}

# ..............................................................................repeated-measures ANOVA + sphericity corrections
compute_anova <- function(df) {
  anova_out <- ezANOVA(
    data       = df,
    dv         = Max_T,
    wid        = Subject,
    within     = Subfield_Name,
    type       = 3,
    return_aov = TRUE
  )
  list(
    table       = as.data.frame(anova_out$ANOVA),
    sphericity  = as.data.frame(anova_out$`Mauchly's Test for Sphericity`),
    corrections = as.data.frame(anova_out$`Sphericity Corrections`)
  )
}

# ..............................................................................post-hoc pairwise comparisons via mixed model
# NOTE: this uses a different error structure from ezANOVA above (a random
# intercept rather than a subject by subfield error term), so the post-hoc
# tests are not strictly the decomposition of that omnibus F
compute_posthoc <- function(df) {
  model <- lmer(Max_T ~ Subfield_Name + (1 | Subject), data = df)
  em    <- emmeans(model, pairwise ~ Subfield_Name, adjust = "fdr")
  list(
    means     = as.data.frame(em$emmeans),
    contrasts = as.data.frame(em$contrasts)
  )
}

# ..............................................................................Cohen's d vs 0 per subfield
# CI taken from effectsize (non-central t, asymmetric) rather than the normal
# approximation d +/- 1.96 * sqrt(1/n + d^2/(2n)) used previously
compute_cohens_d <- function(df) {
  df %>%
    dplyr::group_by(Subfield_Name) %>%
    dplyr::group_modify(~ {
      es <- effectsize::cohens_d(.x$Max_T, mu = 0, ci = 0.95)
      tibble(cohens_d = es$Cohens_d, ci_low = es$CI_low, ci_upp = es$CI_high)
    }) %>%
    dplyr::ungroup()
}

# ..............................................................................check the design before ezANOVA
# ezANOVA needs a complete balanced design and fails on any missing cell, so
# report the offending subjects rather than letting it error out
check_design <- function(df, label) {
  cells <- df %>%
    dplyr::count(Subject, Subfield_Name) %>%
    dplyr::count(Subject, name = "n_subfields")
  
  incomplete <- cells %>% dplyr::filter(n_subfields != length(subfield_order))
  if (nrow(incomplete) > 0) {
    message(label, ": subject(s) without all ", length(subfield_order),
            " subfields:")
    print(incomplete)
  }
  
  duplicated_cells <- df %>%
    dplyr::count(Subject, Subfield_Name) %>%
    dplyr::filter(n > 1)
  if (nrow(duplicated_cells) > 0) {
    message(label, ": duplicated subject by subfield cells:")
    print(duplicated_cells)
  }
  
  invisible(df)
}

# ..............................................................................full pipeline for one condition map × hemisphere combination
run_pipeline <- function(data, condition, hemisphere) {
  label <- paste(condition, hemisphere, sep = " / ")
  
  df <- data %>%
    dplyr::filter(
      Condition_Map == condition,
      Hemisphere    == hemisphere,
      Subfield_Name %in% subfield_order
    ) %>%
    dplyr::mutate(
      Subject       = factor(Subject),              # ez requires a factor wid
      Subfield_Name = factor(Subfield_Name, levels = subfield_order)
    ) %>%
    droplevels()
  
  if (nrow(df) == 0) {
    stop("No rows for ", label, "; check Condition_Map and Hemisphere labels.")
  }
  check_design(df, label)
  
  tag <- function(x) {
    dplyr::mutate(x,
                  Condition  = condition,
                  Hemisphere = unname(hemisphere_labels[hemisphere]))
  }
  
  list(
    condition  = condition,
    hemisphere = hemisphere,
    stats      = tag(compute_stats(df)),
    anova      = compute_anova(df),
    posthoc    = compute_posthoc(df),
    cohens_d   = tag(compute_cohens_d(df)),
    data       = df
  )
}

# ..............................................................................assemble bilateral Cohen's d plot data for one condition map
assemble_plot_data <- function(left_res, right_res) {
  cohens_d_combined <- dplyr::bind_rows(left_res$cohens_d, right_res$cohens_d) %>%
    dplyr::mutate(
      Subfield_Name = factor(Subfield_Name, levels = subfield_order),
      Hemisphere    = factor(Hemisphere,    levels = c("LH", "RH"))
    )
  
  stats_combined <- dplyr::bind_rows(left_res$stats, right_res$stats)
  
  cohens_d_combined %>%
    dplyr::left_join(
      stats_combined %>% dplyr::select(Subfield_Name, Hemisphere, n, t_stat, df),
      by = c("Subfield_Name", "Hemisphere")
    )
}

# ..............................................................................plot
PANEL_W  <- 22      # cm — plot width
PANEL_H  <- 14      # cm — plot height
PLOT_DPI <- 600     # dpi — sufficient for print; use 600 if journal requires it

BASE_SIZE <- 20     # pt — increase if text looks small after inserting into PPT

theme_cb <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size + 2, face = "bold"),
      
      axis.text        = element_text(size = base_size,     colour = "black"),
      axis.title.x     = element_text(size = base_size + 1, margin = margin(t = 6)),
      axis.title.y     = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = 1.5),
      axis.ticks       = element_line(colour = "black", linewidth = 1.5),
      
      legend.position  = "right",
      legend.title     = element_text(size = base_size + 2, face = "bold"),
      legend.text      = element_text(size = base_size),
      legend.key.size  = unit(0.5,  "cm"),
      legend.spacing.y = unit(0.15, "cm"),
      # panels
      panel.spacing    = unit(1.2, "lines"),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

make_forest_plot <- function(plot_data, title = NULL) {
  
  # axis range derived from the data. limits = c(-0.5, 6.0) on
  # scale_x_continuous DELETED any estimate or CI bound outside that window
  # instead of clipping the view
  x_lo <- min(c(plot_data$ci_low, 0), na.rm = TRUE)
  x_hi <- max(c(plot_data$ci_upp, REF_LARGE), na.rm = TRUE)
  pad  <- 0.05 * (x_hi - x_lo)
  
  ggplot(
    plot_data,
    aes(x = cohens_d, y = Subfield_Name, colour = Subfield_Name)
  ) +
    
    geom_vline(xintercept = 0,          linetype = "solid",  colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = REF_SMALL,  linetype = "dashed", colour = "grey65", linewidth = 0.3) +
    geom_vline(xintercept = REF_MEDIUM, linetype = "dashed", colour = "grey65", linewidth = 0.3) +
    geom_vline(xintercept = REF_LARGE,  linetype = "dashed", colour = "grey65", linewidth = 0.3) +
    
    # geom_linerange replaces geom_errorbarh, which is deprecated in
    # ggplot2 >= 3.5 and only had height = 0 to remove its caps anyway
    geom_linerange(
      aes(xmin = ci_low, xmax = ci_upp),
      linewidth = 1.0,
      alpha     = 0.75
    ) +
    geom_point(size = 3.5) +
    
    scale_colour_manual(values = subfield_colours, name = "Subfield") +
    scale_y_discrete(limits = rev(subfield_order)) +
    scale_x_continuous(
      name   = "Effect size (Cohen's d)",
      breaks = scales::breaks_pretty(n = 6),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    coord_cartesian(xlim = c(x_lo - pad, x_hi + pad)) +
    facet_wrap(~ Hemisphere, ncol = 2) +
    labs(title = title) +
    theme_cb() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0,
                                margin = margin(b = 4))
    )
}

# ..............................................................................load data
T_map <- readr::read_csv(file.path(stats_dir, "hippocampus_stats_TNT.csv"),
                         show_col_types = FALSE)

missing_subfields <- setdiff(subfield_order, unique(T_map$Subfield_Name))
if (length(missing_subfields) > 0) {
  warning("Subfield(s) absent from the data: ",
          paste(missing_subfields, collapse = ", "), call. = FALSE)
}

# ..............................................................................run pipeline
# spmT_0001 (Think > No-Think), spmT_0002 (No-Think > Think)
run_grid <- tidyr::expand_grid(
  condition  = c("spmT_0001", "spmT_0002"),
  hemisphere = c("left", "right")
)

all_results <- purrr::pmap(run_grid, ~ run_pipeline(T_map, ..1, ..2)) %>%
  setNames(paste(run_grid$condition, run_grid$hemisphere, sep = "_"))

plot_data_t1 <- assemble_plot_data(all_results$spmT_0001_left,
                                   all_results$spmT_0001_right)
plot_data_t2 <- assemble_plot_data(all_results$spmT_0002_left,
                                   all_results$spmT_0002_right)

# ..............................................................................save plots
plot_t1 <- make_forest_plot(plot_data_t1, title = "Memory Enhancement")
plot_t2 <- make_forest_plot(plot_data_t2, title = "Memory Reduction")

print(plot_t1)
print(plot_t2)

ggsave(file.path(out_dir, "tacmem_cohens_d_TNT_spmT_0001_think_nothink.tiff"),
       plot = plot_t1, units = "cm", width = PANEL_W, height = PANEL_H,
       dpi = PLOT_DPI, compression = "lzw")

ggsave(file.path(out_dir, "tacmem_cohens_d_TNT_spmT_0002_nothink_think.tiff"),
       plot = plot_t2, units = "cm", width = PANEL_W, height = PANEL_H,
       dpi = PLOT_DPI, compression = "lzw")

# ..............................................................................export stats
# the ANOVA tables, sphericity corrections, post-hoc contrasts and effect sizes
# were all computed and then discarded previously; each is now written out
tag_run <- function(x, res) {
  dplyr::mutate(x,
                Condition  = res$condition,
                Hemisphere = unname(hemisphere_labels[res$hemisphere]),
                .before    = 1)
}

readr::write_csv(
  purrr::map(all_results, "stats") %>% purrr::list_rbind(),
  file.path(out_dir, "tacmem_subfield_stats_summary.csv")
)

readr::write_csv(
  purrr::map(all_results, "cohens_d") %>% purrr::list_rbind(),
  file.path(out_dir, "tacmem_subfield_cohens_d.csv")
)

readr::write_csv(
  purrr::map(all_results, ~ tag_run(.x$anova$table, .x)) %>% purrr::list_rbind(),
  file.path(out_dir, "tacmem_subfield_anova_tables.csv")
)

readr::write_csv(
  purrr::map(all_results, ~ tag_run(.x$anova$corrections, .x)) %>% purrr::list_rbind(),
  file.path(out_dir, "tacmem_subfield_sphericity_corrections.csv")
)

readr::write_csv(
  purrr::map(all_results, ~ tag_run(.x$posthoc$contrasts, .x)) %>% purrr::list_rbind(),
  file.path(out_dir, "tacmem_subfield_posthoc_contrasts.csv")
)

message("Output written to: ", out_dir)