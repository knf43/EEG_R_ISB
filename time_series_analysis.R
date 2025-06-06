## Libraries
library(tidyverse)
library(janitor)
library(stringr)
library(janitor)
library(eegUtils)
library(patchwork)
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)
library(ggplot2)
library(ggpubr)
library(here)
library(fs)
library(eeguana)
library(eegUtils)
library(afex)
library(emmeans)
library(remotes)
library(akima)
library(glue)

# Import data
# ---- 1. Fix path to .vhdr files ----
vhdr_paths <- list.files("EEG_data", pattern = "\\.vhdr$", recursive = TRUE, full.names = TRUE)

# ---- 2. Updated read_ascii_erp function ----
read_ascii_erp <- function(vhdr_path) {
  lines <- readLines(vhdr_path)
  dat_file_line <- lines[grepl("^DataFile=", lines)]
  dat_filename <- str_trim(str_replace(dat_file_line, "DataFile=", ""))
  dat_path <- normalizePath(file.path(dirname(vhdr_path), dat_filename), mustWork = TRUE)
  
  dat_df <- read_table(dat_path, col_names = FALSE, skip = 0) %>%
    select(where(is.numeric))
  
  ch_lines <- lines[grepl("^Ch\\d=", lines)]
  ch_names <- str_match(ch_lines, "^Ch\\d+=([^,]+)")[, 2]
  
  dat_long <- dat_df %>%
    t() %>%
    as.data.frame() %>%
    setNames(paste0("V", seq_along(.))) %>%
    tibble::rownames_to_column("time_index") %>%
    mutate(time = seq(-200, by = 2, length.out = n())) %>%
    pivot_longer(cols = starts_with("V"), names_to = "channel_idx", values_to = "amplitude") %>%
    mutate(
      channel = ch_names[as.integer(str_remove(channel_idx, "V"))],
      amplitude = as.numeric(amplitude)
    ) %>%
    select(time, channel, amplitude)
  
  return(dat_long)
}

# ---- 3. Extract metadata from filenames ----
vhdr_df <- tibble(
  filepath = vhdr_paths,
  filename = basename(vhdr_paths),
  participant = str_extract(filename, "stress\\d+_L\\d+"),
  area = case_when(
    str_detect(filename, "F0") ~ "F0",
    str_detect(filename, "suffix") ~ "suffix",
    TRUE ~ NA_character_
  ),
  match = case_when(
    str_detect(filename, "mismatch") ~ "mismatch",
    str_detect(filename, "match") ~ "match",
    TRUE ~ NA_character_
  ),
  stress = case_when(
    str_detect(filename, "stressed") ~ "stressed",
    str_detect(filename, "unstressed") ~ "unstressed",
    TRUE ~ NA_character_
  )
)

# ---- 4. Read ERP data and merge with metadata ----
erp_all <- vhdr_df %>%
  mutate(erp_data = map(filepath, read_ascii_erp)) %>%
  select(-filepath, -filename) %>%
  unnest(erp_data) %>%
  mutate(
    L1_L2 = if_else(str_detect(participant, "_L1"), "L1", "L2"),
    Participant_ID = str_extract(participant, "\\d+_L\\d+")
  )

erp_all <- erp_all %>%
  filter(channel != "1")

# Plots
# ---- Plot ERP Data from F0 Onset ----
erp_plot_data_f0 <- erp_all %>%
  filter(
    area == "F0",
    match %in% c("match", "mismatch"),
    !is.na(channel)
  ) %>%
  mutate(
    channel = factor(channel),
    L1_L2 = factor(L1_L2, levels = c("L1", "L2")),
    amplitude = as.numeric(amplitude)
  ) %>%
  group_by(L1_L2, channel, time, match) %>%
  summarise(mean_amp = mean(amplitude, na.rm = TRUE), .groups = "drop")

ggplot(erp_plot_data_f0, aes(x = time, y = mean_amp, color = match)) +
  geom_line(size = 0.7) +
  facet_wrap(L1_L2 ~ channel, ncol = 6, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_reverse() +
  labs(
    title = "ERP Data from Onset of F0",
    x = "Time (ms)",
    y = "Amplitude (µV)",
    color = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ---- Plot ERP Data from Suffix Onset ----
erp_plot_data_suffix <- erp_all %>%
  filter(
    area == "suffix",
    match %in% c("match", "mismatch"),
    !is.na(channel)
  ) %>%
  mutate(
    channel = factor(channel),
    L1_L2 = factor(L1_L2, levels = c("L1", "L2")),
    amplitude = as.numeric(amplitude)
  ) %>%
  group_by(L1_L2, channel, time, match) %>%
  summarise(mean_amp = mean(amplitude, na.rm = TRUE), .groups = "drop")

ggplot(erp_plot_data_suffix, aes(x = time, y = mean_amp, color = match)) +
  geom_line(size = 0.7) +
  facet_wrap(L1_L2 ~ channel, ncol = 6, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_reverse() +
  labs(
    title = "ERP Data from Onset of Suffix",
    x = "Time (ms)",
    y = "Amplitude (µV)",
    color = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Topographic maps
# Electrode map
chanlocs <- tibble::tibble(
  electrode = c("Fp1", "F3", "F7", "FT9", "FC5", "FC1", "C3", "T7", "TP9",
                "CP5", "CP1", "Pz", "P3", "P7", "O1", "Oz", "O2", "P4",
                "P8", "TP10", "CP6", "CP2", "Cz", "C4", "T8", "FT10",
                "FC6", "FC2", "F4", "F8", "Fp2", "Audio", "Fz"),
  radius = rep(1, 33),
  theta = c(-90, -60, -90, -113, -69, -31, -45, -90, -113,
            -69, -31, 45, -60, -90, -90, 90, 90, 60,
            90, 113, 69, 31, 0, 45, 90, 113,
            69, 31, 60, 90, 90, 0, 45),
  phi = c(-72, -51, -36, -18, -22, -46, 0, 0, 18,
          21, 46, -90, 51, 36, 72, -90, -72, -51,
          -36, -18, -21, -46, 0, 0, 0, 18,
          21, 46, 51, 36, 72, 0, 90)
) %>%
  mutate(
    electrode = toupper(trimws(electrode)),
    theta_rad = theta * pi / 180,
    phi_rad = phi * pi / 180,
    x = radius * cos(phi_rad) * sin(theta_rad),
    y = radius * sin(phi_rad),
    z = radius * cos(phi_rad) * cos(theta_rad)
  )

# Topoplot function (without fixed scale)
plot_topos_by_condition <- function(df, area_name, time_start, time_end, 
                                    condition_var, chanlocs,
                                    palette = "plasma") {
  df <- df %>%
    mutate(channel = toupper(trimws(channel)),
           condition = .data[[condition_var]])
  
  topo_data <- df %>%
    filter(area == area_name,
           time >= time_start, time <= time_end,
           !channel %in% c("1", "AUDIO", NA),
           !is.na(condition)) %>%
    group_by(L1_L2, condition, channel) %>%
    summarise(amplitude = mean(as.numeric(amplitude), na.rm = TRUE), .groups = "drop") %>%
    rename(electrode = channel) %>%
    mutate(electrode = toupper(trimws(electrode))) %>%
    semi_join(chanlocs, by = "electrode")
  
  plots <- list()
  combos <- unique(topo_data[c("L1_L2", "condition")])
  
  for (i in seq_len(nrow(combos))) {
    grp <- combos$L1_L2[i]
    cond <- combos$condition[i]
    
    plot_df <- topo_data %>%
      filter(L1_L2 == grp, condition == cond)
    
    if (nrow(plot_df) < 3) {
      message(glue::glue("Skipping {grp} - {cond}: <3 matching electrodes"))
      next
    }
    
    message(glue::glue("Plotting {grp} – {cond}"))
    
    p <- topoplot(
      plot_df,
      chanLocs = chanlocs,
      interp_limit = "head",
      contour = FALSE,
      palette = palette
    ) +
      ggtitle(glue("{area_name} {time_start}-{time_end}ms\n{grp} – {cond}"))
    
    plots[[paste(grp, cond, sep = "_")]] <- p
  }
  
  return(plots)
}

# Generate topoplots (auto-scaled)
pran_plots <- plot_topos_by_condition(erp_all, "F0", 136, 280, "stress", chanlocs)
lan_plots  <- plot_topos_by_condition(erp_all, "suffix", 225, 400, "match", chanlocs)
p600_plots <- plot_topos_by_condition(erp_all, "suffix", 500, 800, "match", chanlocs)

# PrAN
pran_plots$L1_stressed
pran_plots$L1_unstressed
pran_plots$L2_stressed
pran_plots$L2_unstressed

# LAN
lan_plots$L1_match
lan_plots$L1_mismatch
lan_plots$L2_match
lan_plots$L2_mismatch

# P600
p600_plots$L1_match
p600_plots$L1_mismatch
p600_plots$L2_match
p600_plots$L2_mismatch


