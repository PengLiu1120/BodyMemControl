# TacMem Training phase individual learning curve
# 1. Learning curves per participant per block, stopping at criterion
# 2. Group mean learning curve by attempt number
# ..............................................................................
# written by P. Liu
# Email: peng.liu@uni-tuebingen.de
# last updated on 23rd July 2026 by P. Liu
# ..............................................................................
# ..............................................................................load packages

library(tidyverse)

subjects <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 24,
              27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 39, 40, 42, 43, 44, 47)

base_dir     <- "/Users/pengliu/Documents/Work/Postdoc/TacMem/workspace/behavioural/exp3"
training_dir <- file.path(base_dir, "rawdata")
fig_dir      <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# NOTE: 0.667 is slightly ABOVE 2/3, so a block that scored exactly 2 of 3
# (0.6667 in the log) failed the `>= criterion` test and was recorded as not
# having reached criterion. Using 2/3 with a small tolerance fixes that.
# Check which convention the log uses before trusting either version.
criterion <- 2 / 3
tol       <- 1e-6

n_blocks <- 6


# ..............................................................................data extraction
# perc_correct at row N of block B = accuracy of attempt N-1 of block B
# last attempt accuracy of block B  = perc_correct at row 1 of block B+1
# for the final block, the last attempt is the last non-missing perc_correct
# attempts after the first one reaching criterion are discarded

read_training_file <- function(subj) {
  
  pattern <- sprintf("%03d_PL_TNT_Training_", subj)
  files   <- list.files(training_dir, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    warning(sprintf("No training file for subject %d.", subj), call. = FALSE)
    return(NULL)
  }
  if (length(files) > 1) {
    warning(sprintf("Subject %d has %d matching files; using %s.",
                    subj, length(files), basename(files[1])), call. = FALSE)
  }
  
  # everything read as character, then converted explicitly. PsychoPy logs have
  # long runs of empty cells, so type guessing can silently make a numeric
  # column logical and turn every value into NA
  read_csv(files[1], col_types = cols(.default = col_character()),
           show_col_types = FALSE)
}

extract_attempts <- function(subj) {
  
  d <- read_training_file(subj)
  if (is.null(d)) return(NULL)
  
  stopifnot(all(c("Block", "perc_correct") %in% names(d)))
  
  d <- d %>% mutate(perc_correct = suppressWarnings(as.numeric(perc_correct)))
  
  # last non-missing value in the file = final attempt of the final block
  final_acc <- d %>% filter(!is.na(perc_correct)) %>% pull(perc_correct) %>% last()
  
  blocks <- d %>%
    filter(!is.na(Block), Block != "") %>%
    transmute(
      # str_extract rather than a fixed "Block (\\d+)/6" pattern, which
      # returned the input unchanged and coerced to NA if the total ever differs
      block_number = as.integer(str_extract(Block, "\\d+")),
      accuracy     = perc_correct
    ) %>%
    group_by(block_number) %>%
    mutate(row_in_block = row_number()) %>%
    ungroup()
  
  # rows 2..n of a block hold attempts 1..n-1 of that block
  within_block <- blocks %>%
    filter(row_in_block >= 2) %>%
    group_by(block_number) %>%
    mutate(order_key = row_in_block) %>%
    ungroup() %>%
    select(block_number, accuracy, order_key)
  
  # row 1 of block b holds the final attempt of block b-1
  carried <- blocks %>%
    filter(row_in_block == 1, block_number > 1) %>%
    transmute(block_number = block_number - 1,
              accuracy,
              order_key = Inf)
  
  last_block <- tibble(block_number = max(blocks$block_number, na.rm = TRUE),
                       accuracy     = final_acc,
                       order_key    = Inf)
  
  out <- bind_rows(within_block, carried, last_block) %>%
    filter(!is.na(accuracy)) %>%
    arrange(block_number, order_key) %>%
    group_by(block_number) %>%
    mutate(attempt_number = row_number()) %>%
    # keep every attempt up to and including the first one at criterion
    mutate(passed = accuracy >= criterion - tol,
           prior_pass = lag(cumsum(passed), default = 0)) %>%
    filter(prior_pass == 0) %>%
    ungroup() %>%
    transmute(subject = subj, block_number, attempt_number, accuracy)
  
  out
}

attempt_data <- map(subjects, extract_attempts) %>%
  list_rbind() %>%
  mutate(block_id = paste(subject, block_number, sep = "_"),
         block_factor = factor(block_number))


# ..............................................................................quality control
# the original dropped NA attempts inside the accumulation loop, which silently
# renumbered later attempts. These checks make any such gap visible instead.

missing_subjects <- setdiff(subjects, unique(attempt_data$subject))
if (length(missing_subjects) > 0) {
  message("No data extracted for subject(s): ",
          paste(missing_subjects, collapse = ", "))
}

qc <- attempt_data %>%
  group_by(subject, block_number) %>%
  summarise(n_attempts   = n(),
            final_acc    = last(accuracy),
            reached_crit = any(accuracy >= criterion - tol),
            .groups      = "drop")

incomplete <- qc %>% filter(!reached_crit)
if (nrow(incomplete) > 0) {
  message(sprintf("%d block(s) never reached criterion:", nrow(incomplete)))
  print(incomplete)
}

short_blocks <- qc %>%
  count(subject) %>%
  filter(n != n_blocks)
if (nrow(short_blocks) > 0) {
  message("Participant(s) without all ", n_blocks, " blocks:")
  print(short_blocks)
}

max_attempts <- max(attempt_data$attempt_number)
message("Maximum attempts observed in any block: ", max_attempts)


# ..............................................................................plot 1

block_colours <- c("1" = "#e41a1c", "2" = "#ff7f00", "3" = "#daa520",
                   "4" = "#4daf4a", "5" = "#377eb8", "6" = "#984ea3")

theme_training <- function() {
  theme_minimal() +
    theme(
      panel.grid      = element_blank(),
      axis.line       = element_line(colour = "black", linewidth = 0.8),
      axis.text       = element_text(size = 8, colour = "black"),
      axis.title      = element_text(size = 10),
      strip.text      = element_text(size = 8, face = "bold"),
      legend.position = "right"
    )
}

p_curves <- ggplot(attempt_data,
                   aes(x = attempt_number, y = accuracy,
                       group = block_id, colour = block_factor)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 1.5, alpha = 0.9) +
  geom_hline(yintercept = criterion, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  facet_wrap(~ subject, ncol = 6,
             labeller = labeller(subject = function(x) paste0("Sub-", x))) +
  # breaks derived from the data, and the range set with coord_cartesian.
  # limits = c(0.5, 3.5) on scale_x_continuous DELETED any attempt beyond the
  # third rather than clipping the view, with only a silent warning
  scale_x_continuous(breaks = seq_len(max_attempts)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     labels = scales::label_percent()) +
  scale_colour_manual(values = block_colours, name = "Block") +
  coord_cartesian(xlim = c(0.8, max_attempts + 0.2), ylim = c(0, 1)) +
  labs(x = "Attempt", y = "Accuracy") +
  theme_training()

ggsave(file.path(fig_dir, "tacmem_training_learning_curves.tiff"), plot = p_curves,
       units = "in", width = 18, height = 12, dpi = 600, compression = "lzw")
message("Saved: training_learning_curves.tiff")


# ..............................................................................plot 2
# attempt_mean was computed in the original but never plotted
attempt_mean <- attempt_data %>%
  group_by(attempt_number) %>%
  summarise(n        = n(),
            mean_acc = mean(accuracy, na.rm = TRUE),
            se_acc   = sd(accuracy, na.rm = TRUE) / sqrt(n()),
            .groups  = "drop")

print(attempt_mean)

p_mean <- ggplot(attempt_mean, aes(x = attempt_number, y = mean_acc)) +
  geom_ribbon(aes(ymin = mean_acc - se_acc, ymax = mean_acc + se_acc),
              fill = "grey70", alpha = 0.4) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = criterion, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  scale_x_continuous(breaks = seq_len(max_attempts)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = scales::label_percent()) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Attempt", y = "Accuracy") +
  theme_training()

ggsave(file.path(fig_dir, "tacmem_training_mean_curve.tiff"), plot = p_mean,
       units = "in", width = 6, height = 5, dpi = 600, compression = "lzw")
message("Saved: training_mean_curve.tiff")

# NOTE: the original header promised a second plot broken down by body part and
# stimulus texture. Nothing in this script extracts those columns from the
# training logs, so it is not implemented here. Tell me the column names in the
# PsychoPy output and it is a short addition.