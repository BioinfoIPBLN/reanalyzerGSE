#!/usr/bin/env Rscript
# R_gantt_chart.R
# Usage: R_gantt_chart.R <step_times.tsv> <output_pdf>
#
# Renders the final pipeline Gantt chart from the step_times.tsv log.
# Called at the very end of the pipeline so all steps are captured.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: R_gantt_chart.R <step_times.tsv> <output_pdf>")
}

step_times_file <- args[1]
output_pdf      <- args[2]

if (!file.exists(step_times_file) || file.info(step_times_file)$size == 0) {
  cat("step_times.tsv not found or empty. Skipping Gantt chart.\n")
  quit(status = 0)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

st <- read.delim(step_times_file, header = FALSE, stringsAsFactors = FALSE)
colnames(st) <- c("step", "epoch", "event")

# Pivot start/end into wide format
st_start <- st[st$event == "start", c("step", "epoch"), drop = FALSE]
st_end   <- st[st$event == "end",   c("step", "epoch"), drop = FALSE]
colnames(st_start) <- c("step", "start_epoch")
colnames(st_end)   <- c("step", "end_epoch")
st_wide <- merge(st_start, st_end, by = "step", all = TRUE)

# Remove steps with missing start or end
st_wide <- st_wide[!is.na(st_wide$start_epoch) & !is.na(st_wide$end_epoch), ]

if (nrow(st_wide) == 0) {
  cat("No complete start/end pairs found. Skipping Gantt chart.\n")
  quit(status = 0)
}

st_wide$duration_min <- (st_wide$end_epoch - st_wide$start_epoch) / 60

# Convert epochs to POSIXct for display
st_wide$start_time <- as.POSIXct(st_wide$start_epoch, origin = "1970-01-01")
st_wide$end_time   <- as.POSIXct(st_wide$end_epoch,   origin = "1970-01-01")

# Order by start time
st_wide <- st_wide[order(st_wide$start_epoch), ]
st_wide$step <- factor(st_wide$step, levels = rev(st_wide$step))

# Total pipeline time
total_min <- round((max(st_wide$end_epoch) - min(st_wide$start_epoch)) / 60, 1)

# Dynamic height based on number of steps
n_steps <- nrow(st_wide)
plot_height <- max(5, 1.5 + n_steps * 0.6)

p_gantt <- ggplot(st_wide, aes(y = step)) +
  geom_segment(aes(x = start_time, xend = end_time, yend = step, color = step),
               linewidth = 6) +
  geom_text(aes(x = start_time + (end_time - start_time)/2,
                label = paste0(round(duration_min, 1), " min")),
            size = 2.5, color = "black") +
  scale_x_datetime(date_labels = "%H:%M") +
  labs(title = "Pipeline Step Timing (Gantt Chart)",
       subtitle = paste0("Total wall-clock: ", total_min, " min"),
       x = "Wall-clock time", y = "") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"))

pdf(output_pdf, width = 10, height = plot_height)
print(p_gantt)
dev.off()

cat(sprintf("Gantt chart saved to: %s\n", output_pdf))
