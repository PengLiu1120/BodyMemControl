# TacMem mediation analyses
#   Family A: vlPFC - right insula
#   Family B: vlPFC - SMA (control region)
#   Family C: vlPFC - PCu (control region)
# ..............................................................................
# written by P. Liu
# Email: peng.liu@uni-tuebingen.de
# last updated on 23rd July 2026 by P. Liu
# ..............................................................................
# ..............................................................................load packages
# mediation depends on MASS, which exports its own select(). Attach mediation
# FIRST so that dplyr::select wins the conflict. The dplyr:: prefixes in
# prepare_roi_data() guard against the same clash if this order ever changes.
library(mediation)
library(tidyverse)

# ..............................................................................specify directory
mri_dir    <- "/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/MRI/stats"
behav_dir  <- "/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/behavioural/exp3"
out_dir    <- file.path(mri_dir, "mediation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_sims    <- 5000        # bootstrap resamples per analysis
boot_type <- "bca"
seed      <- 42

# ..............................................................................suppression index
accuracy <- read_csv(file.path(behav_dir, "data.csv"), show_col_types = FALSE) %>%
  rename(Subject = subject)

behav_summary <- accuracy %>%
  group_by(Subject) %>%
  summarise(
    Baseline_Acc      = mean(accuracy[condition == "B"],  na.rm = TRUE),
    NoThink_Acc       = mean(accuracy[condition == "NT"], na.rm = TRUE),
    Suppression_Index = Baseline_Acc - NoThink_Acc,
    .groups           = "drop"
  ) %>%
  mutate(Subject = sprintf("sub-%02d", as.integer(Subject)))

# ..............................................................................read ROIs
prepare_roi_data <- function(csv_name, roi_map, family) {
  
  roi_raw <- read_csv(file.path(mri_dir, csv_name), show_col_types = FALSE) %>%
    mutate(Subject = as.character(Subject))
  
  dat <- behav_summary %>%
    inner_join(roi_raw, by = "Subject") %>%
    dplyr::select(Subject, Baseline_Acc, NoThink_Acc, Suppression_Index,
                  ROI, Max_T) %>%
    pivot_wider(names_from = ROI, values_from = Max_T)
  
  missing_rois <- setdiff(unname(roi_map), names(dat))
  if (length(missing_rois) > 0) {
    stop(sprintf("%s: ROI label(s) not found in %s: %s",
                 family, csv_name, paste(missing_rois, collapse = ", ")))
  }
  
  dropped <- setdiff(behav_summary$Subject, dat$Subject)
  message(sprintf("%s: %d behavioural, %d ROI, %d joined%s",
                  family, nrow(behav_summary), n_distinct(roi_raw$Subject),
                  nrow(dat),
                  if (length(dropped) > 0)
                    paste0(" | dropped: ", paste(dropped, collapse = ", ")) else ""))
  if (nrow(dat) == 0) stop(family, ": join produced no rows; check Subject formats.")
  
  dat %>%
    dplyr::rename(all_of(roi_map)) %>%
    dplyr::mutate(across(all_of(names(roi_map)),
                         ~ as.numeric(scale(.x)), .names = "{.col}_z"),
                  SI_z = as.numeric(scale(Suppression_Index)))
}

run_mediation <- function(data, x_var, m_var, y_var, label, family,
                          sims = n_sims, verbose = TRUE) {
  
  # do.call with quote(data) keeps `data` as a symbol in each model call, which
  # is what mediate() needs when it refits the models on bootstrap resamples.
  # Do not simplify these to plain lm() calls with a piped data frame.
  fit <- function(rhs) {
    do.call(lm, list(formula = as.formula(paste(y_var, "~", rhs)),
                     data    = quote(data)))
  }
  
  model_total <- fit(x_var)                              # c path
  model_y     <- fit(paste(x_var, "+", m_var))           # b and c' paths
  model_m     <- do.call(lm, list(                       # a path
    formula = as.formula(paste(m_var, "~", x_var)),
    data    = quote(data)))
  
  set.seed(seed)
  med_result <- mediation::mediate(
    model.m      = model_m,
    model.y      = model_y,
    treat        = x_var,
    mediator     = m_var,
    boot         = TRUE,
    boot.ci.type = boot_type,
    sims         = sims
  )
  
  if (verbose) {
    cat("\n================================================================\n")
    cat(label, "\n")
    cat(sprintf("X = %s | M = %s | Y = %s | n = %d\n",
                x_var, m_var, y_var, nrow(data)))
    cat("================================================================\n")
    cat("\n--- Total effect (c path): X -> Y ---\n")
    print(summary(model_total)$coefficients)
    cat("\n--- Mediator model (a path): X -> M ---\n")
    print(summary(model_m)$coefficients)
    cat("\n--- Outcome model (b and c' paths): X + M -> Y ---\n")
    print(summary(model_y)$coefficients)
    cat("\n--- Mediation ---\n")
    print(summary(med_result))
  }
  
  list(label = label, family = family, n = nrow(data),
       x_var = x_var, m_var = m_var, y_var = y_var,
       model_total = model_total, model_m = model_m, model_y = model_y,
       med_result = med_result)
}

extract_mediation_summary <- function(res) {
  
  med <- res$med_result
  
  tibble(
    Family         = res$family,
    Analysis       = res$label,
    n              = res$n,
    ACME_estimate  = round(med$d0,        3),
    ACME_CI_lower  = round(med$d0.ci[1],  3),
    ACME_CI_upper  = round(med$d0.ci[2],  3),
    ACME_p         = round(med$d0.p,      3),
    ADE_estimate   = round(med$z0,        3),
    ADE_CI_lower   = round(med$z0.ci[1],  3),
    ADE_CI_upper   = round(med$z0.ci[2],  3),
    ADE_p          = round(med$z0.p,      3),
    Total_estimate = round(med$tau.coef,  3),
    Total_CI_lower = round(med$tau.ci[1], 3),
    Total_CI_upper = round(med$tau.ci[2], 3),
    Total_p        = round(med$tau.p,     3),
    Prop_mediated  = round(med$n0,        3),
    Prop_p         = round(med$n0.p,      3),
    
    Inconsistent   = sign(med$d0) != sign(med$z0)
  )
}

sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}

plot_path_text <- function(res) {
  
  a  <- coef(res$model_m)[2]
  b  <- coef(res$model_y)[3]
  cp <- coef(res$model_y)[2]
  cc <- coef(res$model_total)[2]
  
  annotation <- paste0(
    "a path (X -> M):   ", round(a,  3), "\n",
    "b path (M -> Y):   ", round(b,  3), "\n",
    "c' path (direct):  ", round(cp, 3), "\n",
    "c path (total):    ", round(cc, 3), "\n",
    "Indirect (a*b):    ", round(res$med_result$d0,   3), "\n",
    "ACME p =           ", round(res$med_result$d0.p, 3)
  )
  
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = annotation,
             size = 4, hjust = 0.5, vjust = 0.5, family = "mono") +
    labs(title = res$label) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
}

# ..............................................................................plot preparation
plot_path_diagram <- function(res, x_label, m_label, y_label,
                              pos_colour = "#B2182B", neg_colour = "#2166AC") {
  
  coef_m <- summary(res$model_m)$coefficients
  coef_y <- summary(res$model_y)$coefficients
  
  a  <- coef_m[2, 1]; a_p  <- coef_m[2, 4]
  cp <- coef_y[2, 1]; cp_p <- coef_y[2, 4]
  b  <- coef_y[3, 1]; b_p  <- coef_y[3, 4]
  ab <- res$med_result$d0; ab_p <- res$med_result$d0.p
  
  path_colour <- function(v) if (v < 0) neg_colour else pos_colour
  path_type   <- function(v) if (v < 0) "dashed" else "solid"
  
  nodes <- tibble(
    x     = c(1, 2, 3),
    y     = c(1, 2, 1),
    label = c(paste0(x_label, "\n(X)"),
              paste0(m_label, "\n(M)"),
              paste0(y_label, "\n(Y)"))
  )
  
  ggplot() +
    geom_label(data = nodes, aes(x = x, y = y, label = label),
               fill = "white", colour = "black", size = 8, fontface = "bold",
               label.padding = unit(0.6, "lines")) +
    # a path
    geom_segment(aes(x = 1.2, y = 1.2, xend = 1.8, yend = 1.8),
                 arrow = arrow(length = unit(0.3, "cm")),
                 colour = path_colour(a), linewidth = 1.5, linetype = path_type(a)) +
    # b path
    geom_segment(aes(x = 2.2, y = 1.8, xend = 2.8, yend = 1.2),
                 arrow = arrow(length = unit(0.3, "cm")),
                 colour = path_colour(b), linewidth = 1.5, linetype = path_type(b)) +
    # c' path
    geom_segment(aes(x = 1.3, y = 1.0, xend = 2.5, yend = 1.0),
                 arrow = arrow(length = unit(0.3, "cm")),
                 colour = path_colour(cp), linewidth = 1.5, linetype = path_type(cp)) +
    annotate("text", x = 1.3, y = 1.6, size = 7, fontface = "bold",
             colour = path_colour(a),
             label = sprintf("path a\nbeta = %.3f%s", a, sig_stars(a_p))) +
    annotate("text", x = 2.7, y = 1.6, size = 7, fontface = "bold",
             colour = path_colour(b),
             label = sprintf("path b\nbeta = %.3f%s", b, sig_stars(b_p))) +
    annotate("text", x = 2.0, y = 0.85, size = 7, fontface = "bold",
             colour = path_colour(cp),
             label = sprintf("direct path c'\nbeta = %.3f%s", cp, sig_stars(cp_p))) +
    annotate("label", x = 2.0, y = 1.4, size = 7, fontface = "bold",
             fill = "#F7F7F7", colour = "black", label.padding = unit(0.4, "lines"),
             label = sprintf("Indirect effect (ACME) = %.3f%s%s",
                             ab, sig_stars(ab_p),
                             if (sign(ab) != sign(cp)) "\n(inconsistent mediation)" else "")) +
    coord_cartesian(xlim = c(0.6, 3.4), ylim = c(0.5, 2.3)) +
    theme_void()
}

# ..............................................................................vlPFC - insula
data_insula <- prepare_roi_data(
  csv_name = "vlPFC_insula_T_values.csv",
  roi_map  = c(vlPFC = "vlPFC_right", Insula_R = "insula_right"),
  family   = "insula"
)

res_insula_r <- run_mediation(data_insula, "vlPFC_z", "Insula_R_z", "SI_z",
                              "vlPFC - right insula - suppression index", "insula")

summary_insula <- extract_mediation_summary(res_insula_r)
write_csv(summary_insula, file.path(out_dir, "mediation_summary_insula.csv"))

# ..............................................................................vlPFC - SMA
data_sma <- prepare_roi_data(
  csv_name = "vlPFC_SMA_T_values.csv",
  roi_map  = c(vlPFC = "vlPFC_right", SMA_R = "SMA_right"),
  family   = "SMA"
)

res_sma <- run_mediation(data_sma, "vlPFC_z", "SMA_R_z", "SI_z",
                         "vlPFC - right SMA - suppression index", "SMA")

summary_sma <- extract_mediation_summary(res_sma)
write_csv(summary_sma, file.path(out_dir, "mediation_summary_SMA.csv"))

# ..............................................................................vlPFC - PCu
data_pcu <- prepare_roi_data(
  csv_name = "vlPFC_PCu_T_values.csv",
  roi_map  = c(vlPFC = "vlPFC_right", PCu_R = "PCu_right"),
  family   = "PCu"
)

res_pcu <- run_mediation(data_pcu, "vlPFC_z", "PCu_R_z", "SI_z",
                         "vlPFC - right PCu - suppression index", "PCu")

summary_pcu <- extract_mediation_summary(res_pcu)
write_csv(summary_pcu, file.path(out_dir, "mediation_summary_PCu.csv"))

# ..............................................................................summary
all_results <- list(
  insula_right = res_insula_r,
  SMA_right    = res_sma,
  PCu_right    = res_pcu
)

summary_all <- map(all_results, extract_mediation_summary) %>% list_rbind()

cat("\n\n======================= combined summary =======================\n")
print(summary_all, width = Inf)

write_csv(summary_all, file.path(out_dir, "mediation_summary_all.csv"))

# ..............................................................................figure
# the node and arrow path diagram is only produced for the insula family.
# SMA and PCu are control regions and get the mediate() plot and the path
# coefficient text panel only. To add a diagram for another region, add an
# entry here keyed by its name in all_results.
node_labels <- list(
  insula_right = c(x = "Right vlPFC", m = "Right IC", y = "Suppression index")
)

# shared x range so the three mediate() plots are directly comparable. The
# original per-region scripts each used whatever range their own plot chose.
ci_range <- function(med) {
  range(c(med$d0.ci, med$z0.ci, med$tau.ci, med$d0, med$z0, med$tau.coef),
        na.rm = TRUE)
}
all_ranges <- map(all_results, ~ ci_range(.x$med_result))
x_min <- min(map_dbl(all_ranges, 1)) - 0.05
x_max <- max(map_dbl(all_ranges, 2)) + 0.05
cat(sprintf("\nShared x limits for mediation plots: [%.3f, %.3f]\n", x_min, x_max))

for (key in names(all_results)) {
  
  res <- all_results[[key]]
  
  png(file.path(out_dir, sprintf("mediation_plot_%s.png", key)),
      width = 800, height = 600, res = 150)
  plot(res$med_result, main = res$label, xlim = c(x_min, x_max))
  dev.off()
  
  ggsave(file.path(out_dir, sprintf("mediation_path_text_%s.tiff", key)),
         plot = plot_path_text(res),
         units = "in", width = 6, height = 4, dpi = 600, compression = "lzw")
  
  if (key %in% names(node_labels)) {
    lab <- node_labels[[key]]
    ggsave(file.path(out_dir, sprintf("mediation_path_diagram_%s.tiff", key)),
           plot = plot_path_diagram(res, x_label = lab["x"],
                                    m_label = lab["m"], y_label = lab["y"]),
           units = "in", width = 14, height = 8, dpi = 600, compression = "lzw")
  }
  
  message("Saved figures for: ", key)
}

cat("\nAll output written to:", out_dir, "\n")