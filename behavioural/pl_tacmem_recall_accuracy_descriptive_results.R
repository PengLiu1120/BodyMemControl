# TacMem Recall phase accuracy descriptive statistics
# ..............................................................................
# written by P.Liu
# Email: peng.liu@uni-tuebingen.de
# last updated on 23rd July 2026 by P. Liu
# ..............................................................................
# ..............................................................................load packages
library(tidyverse)

# ..............................................................................load data
setwd("/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/behavioural/exp3")
d <- read.csv("data_learned.csv", header = TRUE, stringsAsFactors = FALSE,
              fileEncoding = "UTF-8-BOM")

# ..............................................................................descriptives from participant-level means
describe_by <- function(data, group_vars) {
  data %>%
    group_by(across(all_of(c("subject", group_vars)))) %>%
    summarise(subj_acc = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n_subjects = n(),
      mean = mean(subj_acc, na.rm = TRUE),
      sd   = sd(subj_acc,   na.rm = TRUE),
      se   = sd(subj_acc,   na.rm = TRUE) / sqrt(n()),
      median = median(subj_acc, na.rm = TRUE),
      min  = min(subj_acc,  na.rm = TRUE),
      max  = max(subj_acc,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(c(mean, sd, se, median, min, max), ~ round(.x, 4)))
}

# ..............................................................................1. overall (per participant mean, then across)
cat("================================================================\n")
cat("1. Overall recall accuracy (participant-level)\n")
cat("================================================================\n")
overall <- d %>%
  group_by(subject) %>%
  summarise(subj_acc = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
  summarise(
    n_subjects = n(),
    mean = mean(subj_acc), sd = sd(subj_acc),
    se = sd(subj_acc)/sqrt(n()),
    median = median(subj_acc),
    min = min(subj_acc), max = max(subj_acc)
  ) %>%
  mutate(across(c(mean, sd, se, median, min, max), ~ round(.x, 4)))
print(as.data.frame(overall), row.names = FALSE)

# ..............................................................................2. by condition
cat("\n================================================================\n")
cat("2. By condition\n")
cat("================================================================\n")
by_cond <- describe_by(d, "condition")
print(as.data.frame(by_cond), row.names = FALSE)

# ..............................................................................3. by condition x body part
cat("\n================================================================\n")
cat("4. By condition x body part\n")
cat("================================================================\n")
by_cond_bp <- describe_by(d, c("condition", "bodypart"))
print(as.data.frame(by_cond_bp), row.names = FALSE)

# ..............................................................................4. by condition x stimulus texture
cat("\n================================================================\n")
cat("5. By condition x stimulus texture\n")
cat("================================================================\n")
by_cond_stim <- describe_by(d, c("condition", "stimulus"))
print(as.data.frame(by_cond_stim), row.names = FALSE)

# ..............................................................................5. by body part x stimulus texture
cat("\n================================================================\n")
cat("6. By body part x stimulus texture\n")
cat("================================================================\n")
by_bp_stim <- describe_by(d, c("bodypart", "stimulus"))
print(as.data.frame(by_bp_stim), row.names = FALSE)

# ..............................................................................save all tables to one CSV (stacked with labels)
all_desc <- bind_rows(
  by_cond      %>% mutate(grouping = "condition"),
  by_cond_bp   %>% mutate(grouping = "condition x bodypart"),
  by_cond_stim %>% mutate(grouping = "condition x stimulus"),
  by_bp_stim   %>% mutate(grouping = "bodypart x stimulus")
)