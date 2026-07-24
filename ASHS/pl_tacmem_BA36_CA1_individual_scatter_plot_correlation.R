# TacMem BA36 and CA1 Max T value correlation matrix — hemisphere side by side
# Layout: rows = contrasts, columns = left | right hemisphere
# ..............................................................................
# written by P.Liu
# Email: peng.liu@uni-tuebingen.de
# last updated on 23rd July 2026 by P. Liu
# ..............................................................................
# ..............................................................................load packages
library(tidyverse)
library(patchwork)

# ..............................................................................load data
stats_dir <- "/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/MRI/stats"
out_dir   <- file.path(stats_dir, "MTL")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

hipp_tnt    <- readr::read_csv(file.path(stats_dir, "hippocampus_stats_TNT.csv"),
                               show_col_types = FALSE)
hipp_recall <- readr::read_csv(file.path(stats_dir, "hippocampus_stats_recall.csv"),
                               show_col_types = FALSE)

# ..............................................................................combine
hipp_all <- bind_rows(hipp_tnt, hipp_recall)

required_cols <- c("Subject", "Hemisphere", "Analysis", "Condition_Map",
                   "Subfield_Name", "Max_T")
stopifnot(all(required_cols %in% names(hipp_all)))

# ..............................................................................define contrasts
#   TNT    spmT_0001 Think > No-Think      spmT_0002 No-Think > Think
#   recall spmT_0001 Think > Baseline      spmT_0002 No-Think > Baseline
contrast_map <- c(
  "TNT_spmT_0001"    = "Enhancement",
  "TNT_spmT_0002"    = "Suppression",
  "recall_spmT_0001" = "Enhanced Recall",
  "recall_spmT_0002" = "Suppressed Recall"
)

# ..............................................................................define colours
contrast_colours <- c(
  "Enhancement"       = "#1d4e89",
  "Suppression"       = "#a50026",
  "Enhanced Recall"   = "#1a7a3c",
  "Suppressed Recall" = "#b35900"
)

stopifnot(all(contrast_map %in% names(contrast_colours)))

# ..............................................................................extract BA36 and CA1 — keep hemisphere separate
data_wide <- hipp_all %>%
  dplyr::filter(Subfield_Name %in% c("BA36", "CA1")) %>%
  dplyr::mutate(contrast_key = paste0(Analysis, "_", Condition_Map)) %>%
  dplyr::filter(contrast_key %in% names(contrast_map)) %>%
  dplyr::select(Subject, Hemisphere, contrast_key, Subfield_Name, Max_T)

# pivot_wider silently produces list columns if a cell is duplicated, so check
duplicate_cells <- data_wide %>%
  dplyr::count(Subject, Hemisphere, contrast_key, Subfield_Name) %>%
  dplyr::filter(n > 1)
if (nrow(duplicate_cells) > 0) {
  message("Duplicated subject by cell entries:")
  print(duplicate_cells)
  stop("Resolve duplicates before pivoting.")
}

data_wide <- data_wide %>%
  tidyr::pivot_wider(names_from = Subfield_Name, values_from = Max_T) %>%
  dplyr::mutate(
    # explicit lookup rather than recode(), which leaves unmatched values in
    # place instead of failing
    contrast_label = unname(contrast_map[contrast_key]),
    contrast_label = factor(contrast_label, levels = unname(contrast_map))
  )

if (anyNA(data_wide$contrast_label)) {
  stop("Unmapped contrast_key: ",
       paste(unique(data_wide$contrast_key[is.na(data_wide$contrast_label)]),
             collapse = ", "))
}
if (!all(c("BA36", "CA1") %in% names(data_wide))) {
  stop("BA36 or CA1 missing after pivot; check Subfield_Name labels.")
}

# ..............................................................................global axis ranges
# shared across every panel so the four contrasts stay comparable. Both are
# derived from the data; the y range was previously fixed at c(0, 2), which
# hid any CA1 value above 2 with no indication
x_global <- range(data_wide$BA36, na.rm = TRUE)
y_global <- range(data_wide$CA1,  na.rm = TRUE)

# ..............................................................................run correlations
# one cor.test per cell rather than three identical calls per row
cor_results <- data_wide %>%
  dplyr::group_by(contrast_key, contrast_label, Hemisphere) %>%
  dplyr::group_modify(~ {
    ok <- stats::complete.cases(.x$BA36, .x$CA1)
    if (sum(ok) < 3) {
      return(tibble(r = NA_real_, p_value = NA_real_, df = NA_real_,
                    n = sum(ok)))
    }
    ct <- cor.test(.x$BA36[ok], .x$CA1[ok], method = "pearson")
    tibble(r = unname(ct$estimate), p_value = ct$p.value,
           df = unname(ct$parameter), n = sum(ok))
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(phase = sub("_.*$", "", contrast_key), .before = 1) %>%
  # FDR runs WITHIN each phase: four cells for TNT (two contrasts by two
  # hemispheres) and four for recall. TNT and recall are separate phases and so
  # are treated as separate families.
  dplyr::group_by(phase) %>%
  dplyr::mutate(
    p_fdr      = p.adjust(p_value, method = "fdr"),
    sig_label  = dplyr::case_when(
      p_fdr < .001 ~ "***",
      p_fdr < .01  ~ "**",
      p_fdr < .05  ~ "*",
      TRUE         ~ "ns"
    ),
    annotation = sprintf("r = %.3f\np(FDR) = %.3f", r, p_fdr)
  ) %>%
  dplyr::ungroup()

print(cor_results)
readr::write_csv(cor_results,
                 file.path(out_dir, "tacmem_BA36_CA1_correlation_summary_by_hemi.csv"))

# ..............................................................................function to create one scatter panel
make_panel <- function(ckey, hemi) {
  
  cor_row <- cor_results %>%
    dplyr::filter(contrast_key == ckey, Hemisphere == hemi)
  
  d <- data_wide %>%
    dplyr::filter(contrast_key == ckey, Hemisphere == hemi)
  
  # fail with a readable message rather than propagating character(0)
  if (nrow(cor_row) != 1 || nrow(d) == 0) {
    stop("No data for ", ckey, " / ", hemi,
         "; check Analysis, Condition_Map and Hemisphere labels.")
  }
  
  clabel <- as.character(cor_row$contrast_label)
  annot  <- cor_row$annotation
  col    <- unname(contrast_colours[clabel])
  
  hemi_label <- ifelse(hemi == "left", "LH", "RH")
  
  ggplot(d, aes(x = BA36, y = CA1)) +
    geom_point(size = 2.5, alpha = 0.75, colour = col) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = col, fill = col,
                alpha = 0.15, linewidth = 0.8) +
    coord_cartesian(xlim = x_global, ylim = y_global) +
    # annotation placed relative to the shared ranges, so it follows the axes
    annotate("text",
             x     = x_global[1] + diff(x_global) * 0.05,
             y     = y_global[2] - diff(y_global) * 0.08,
             label = annot,
             hjust = 0, vjust = 1,
             size  = 5, colour = "black") +
    labs(
      title = paste0(clabel, "\n", hemi_label),
      x     = "BA36 Max T",
      y     = "CA1 Max T"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(colour = "black", linewidth = 1.0),
      plot.title        = element_text(size = 24, face = "bold", hjust = 0.5),
      axis.title.x      = element_text(size = 24),
      axis.title.y      = element_text(size = 24),
      axis.text.x       = element_text(size = 20),
      axis.text.y       = element_text(size = 20)
    )
}

# ..............................................................................build one 2x2 figure (rows = contrast, cols = hemisphere)
build_figure <- function(keys) {
  stopifnot(length(keys) == 2)
  (make_panel(keys[1], "left") | make_panel(keys[1], "right")) /
    (make_panel(keys[2], "left") | make_panel(keys[2], "right"))
}

p_tnt    <- build_figure(c("TNT_spmT_0001", "TNT_spmT_0002"))
p_recall <- build_figure(c("recall_spmT_0001", "recall_spmT_0002"))

# ..............................................................................save
correlation_figures <- list(
  "tacmem_BA36_CA1_correlation_TNT.tiff"    = p_tnt,
  "tacmem_BA36_CA1_correlation_recall.tiff" = p_recall
)

for (fname in names(correlation_figures)) {
  ggsave(file.path(out_dir, fname), plot = correlation_figures[[fname]],
         units = "in", width = 10, height = 10, dpi = 600, compression = "lzw")
  message("Saved: ", fname)
}

message("All outputs saved to: ", out_dir)