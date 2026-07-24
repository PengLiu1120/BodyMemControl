# TacMem Recall phase encoded accuracy condition one-tailed t test
# ..............................................................................
# written by P. Liu
# Email: peng.liu@uni-tuebingen.de
# last updated on 23rd July 2026 by P. Liu
# ..............................................................................
# ..............................................................................load packages
library(tidyverse)

# ..............................................................................specify directory
data_dir    <- "/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/behavioural/exp3"
results_dir <- file.path(data_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ..............................................................................read data
condition_labels <- c(T = "Enhancement", B = "Baseline", NT = "Suppression")

accuracy_trial <- read_csv(file.path(data_dir, "data_learned.csv"),
                           show_col_types = FALSE) %>%
  mutate(condition = unname(condition_labels[condition]))

if (anyNA(accuracy_trial$condition)) {
  stop("Unrecognised value in `condition`; check condition_labels.")
}

cond_means <- accuracy_trial %>%
  group_by(subject, condition) %>%
  summarise(acc = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = acc)

stopifnot(all(condition_labels %in% names(cond_means)))
message(sprintf("%d participants, %d with a complete set of conditions.",
                nrow(cond_means),
                sum(complete.cases(select(cond_means, all_of(unname(condition_labels)))))))

# ..............................................................................paired t test with Cohen's dz
paired_test <- function(data, high, low, alternative = "two.sided",
                        label = NULL, conf_level = 0.95, verbose = TRUE) {
  
  label <- label %||% sprintf("%s vs %s", high, low)
  
  pair <- data %>%
    select(subject, high = all_of(high), low = all_of(low)) %>%
    drop_na(high, low)
  
  n_dropped <- nrow(data) - nrow(pair)
  if (n_dropped > 0) {
    warning(sprintf("%s: %d participant(s) dropped for incomplete data.",
                    label, n_dropped), call. = FALSE)
  }
  
  diffs <- pair$high - pair$low
  tt    <- t.test(pair$high, pair$low, paired = TRUE,
                  alternative = alternative, conf.level = conf_level)
  
  dz    <- mean(diffs) / sd(diffs)
  dz_ci <- t.test(diffs, conf.level = conf_level)$conf.int / sd(diffs)
  
  out <- tibble(
    comparison  = label,
    n           = nrow(pair),
    mean_high   = mean(pair$high),
    mean_low    = mean(pair$low),
    mean_diff   = mean(diffs),
    t           = unname(tt$statistic),
    df          = unname(tt$parameter),
    alternative = alternative,
    p           = tt$p.value,
    cohens_dz   = dz,
    dz_ci_low   = dz_ci[1],
    dz_ci_high  = dz_ci[2],
    
    shapiro_p   = shapiro.test(diffs)$p.value,
    wilcoxon_p  = suppressWarnings(
      wilcox.test(pair$high, pair$low, paired = TRUE,
                  alternative = alternative)$p.value)
  )
  
  if (verbose) {
    cat("\n----------------------------------------------------------------\n")
    cat(label, "\n")
    cat("----------------------------------------------------------------\n")
    print(tt)
    cat(sprintf("Mean difference: %.4f\n", out$mean_diff))
    cat(sprintf("Cohen's dz: %.4f [%.4f, %.4f]\n",
                out$cohens_dz, out$dz_ci_low, out$dz_ci_high))
    cat(sprintf("Shapiro-Wilk on differences: p = %.4f | Wilcoxon: p = %.4f\n",
                out$shapiro_p, out$wilcoxon_p))
  }
  
  out
}

# ..............................................................................comparison
comparisons <- tribble(
  ~high,          ~low,          ~alternative, ~label,
  "Baseline",     "Suppression", "greater",    "Suppression vs Baseline (one-tailed, B > NT)",
  "Enhancement",  "Baseline",    "greater",    "Enhancement vs Baseline (one-tailed, T > B)"
)

summary_tbl <- comparisons %>%
  pmap_dfr(function(high, low, alternative, label) {
    paired_test(cond_means, high = high, low = low,
                alternative = alternative, label = label)
  }) %>%
  
  mutate(p_holm = p.adjust(p, method = "holm"), .after = p)

cat("\n\n================ summary ================\n")
print(summary_tbl, width = Inf)

out_file <- file.path(results_dir, "recall_learned_accuracy_condition_t_test.csv")
write_csv(summary_tbl, out_file)
message("Saved: ", out_file)