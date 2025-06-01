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

read_lines("raw_data/Area_F0/Area136_280_F0.txt", n_max = 5)

## Load data
# F0
f0_136_280_raw <- read_fwf("raw_data/Area_F0/Area136_280_F0.txt",
                           col_positions = fwf_empty("raw_data/Area_F0/Area136_280_F0.txt"),
                           trim_ws = TRUE) %>% clean_names()

f0_280_400_raw <- read_fwf("raw_data/Area_F0/Area280_400_F0.txt",
                           col_positions = fwf_empty("raw_data/Area_F0/Area280_400_F0.txt"),
                           trim_ws = TRUE) %>% clean_names()

f0_400_600_raw <- read_fwf("raw_data/Area_F0/Area_400_600_F0.txt",
                           col_positions = fwf_empty("raw_data/Area_F0/Area_400_600_F0.txt"),
                           trim_ws = TRUE) %>% clean_names()

# Suffix 
suffix_235_415_raw <- read_fwf("raw_data/Area_suffix/Area235_415_suffix.txt",
                               col_positions = fwf_empty("raw_data/Area_suffix/Area235_415_suffix.txt"),
                               trim_ws = TRUE) %>% clean_names()

suffix_415_550_raw <- read_fwf("raw_data/Area_suffix/Area415_550_suffix.txt",
                               col_positions = fwf_empty("raw_data/Area_suffix/Area415_550_suffix.txt"),
                               trim_ws = TRUE) %>% clean_names()

suffix_550_680_raw <- read_fwf("raw_data/Area_suffix/Area_550_680_suffix.txt",
                               col_positions = fwf_empty("raw_data/Area_suffix/Area_550_680_suffix.txt"),
                               trim_ws = TRUE) %>% clean_names()

suffix_680_800_raw <- read_fwf("raw_data/Area_suffix/Area680_800_suffix.txt",
                               col_positions = fwf_empty("raw_data/Area_suffix/Area680_800_suffix.txt"),
                               trim_ws = TRUE) %>% clean_names()

## Tidy data
# Remove second row
make_tidy <- function(df) {
  df <- df[-2, ]  # remove second row
  names(df) <- as.character(unlist(df[1, ]))  # set first row as col names
  df <- df[-1, ]  # remove former header row
  return(df)
}

# F0 files
f0_136_280_tidy <- make_tidy(f0_136_280_raw)
f0_280_400_tidy <- make_tidy(f0_280_400_raw)
f0_400_600_tidy <- make_tidy(f0_400_600_raw)

# Suffix files
suffix_235_415_tidy <- make_tidy(suffix_235_415_raw)
suffix_415_550_tidy <- make_tidy(suffix_415_550_raw)
suffix_550_680_tidy <- make_tidy(suffix_550_680_raw)
suffix_680_800_tidy <- make_tidy(suffix_680_800_raw)

# Split column 1
split_participant_info <- function(df) {
  df <- df %>%
    separate(1, into = c("Participant_ID", "L1_L2"), sep = "_", remove = TRUE)
  return(df)
}

# Apply to all tidy data
f0_136_280_tidy <- split_participant_info(f0_136_280_tidy)
f0_280_400_tidy <- split_participant_info(f0_280_400_tidy)
f0_400_600_tidy <- split_participant_info(f0_400_600_tidy)

suffix_235_415_tidy <- split_participant_info(suffix_235_415_tidy)
suffix_415_550_tidy <- split_participant_info(suffix_415_550_tidy)
suffix_550_680_tidy <- split_participant_info(suffix_550_680_tidy)
suffix_680_800_tidy <- split_participant_info(suffix_680_800_tidy)

# Convert ERP value columns to numeric 
# F0 datasets
f0_136_280_tidy <- f0_136_280_tidy %>%
  mutate(across(-c(Participant_ID, L1_L2), as.numeric))

f0_280_400_tidy <- f0_280_400_tidy %>%
  mutate(across(-c(Participant_ID, L1_L2), as.numeric))

f0_400_600_tidy <- f0_400_600_tidy %>%
  mutate(across(-c(Participant_ID, L1_L2), as.numeric))

# Suffix datasets
suffix_235_415_tidy <- suffix_235_415_tidy %>%
  mutate(across(-c(Participant_ID, L1_L2), as.numeric))

suffix_415_550_tidy <- suffix_415_550_tidy %>%
  mutate(across(-c(Participant_ID, L1_L2), as.numeric))

suffix_550_680_tidy <- suffix_550_680_tidy %>%
  mutate(across(-c(Participant_ID, L1_L2), as.numeric))

suffix_680_800_tidy <- suffix_680_800_tidy %>%
  mutate(across(-c(Participant_ID, L1_L2), as.numeric))

str(suffix_680_800_tidy)

# Save tidy datasets
write.csv(suffix_235_415_tidy, "tidy_data/suffix_235_415_tidy.csv", row.names = FALSE)
write.csv(suffix_415_550_tidy, "tidy_data/suffix_415_550_tidy.csv", row.names = FALSE)
write.csv(suffix_550_680_tidy, "tidy_data/suffix_550_680_tidy.csv", row.names = FALSE)
write.csv(suffix_680_800_tidy, "tidy_data/suffix_680_800_tidy.csv", row.names = FALSE)

write.csv(f0_136_280_tidy, "tidy_data/f0_136_280_tidy.csv", row.names = FALSE)
write.csv(f0_280_400_tidy, "tidy_data/f0_280_400_tidy.csv", row.names = FALSE)
write.csv(f0_400_600_tidy, "tidy_data/f0_400_600_tidy.csv", row.names = FALSE)

## Combine datasets
insert_cols <- function(df, area_value, time_value) {
  df %>%
    mutate(area = area_value, time_period = time_value, .before = 3)
}

# Apply to each dataset
f0_136_280_tidy <- insert_cols(f0_136_280_tidy, "F0", "136_280")
f0_280_400_tidy <- insert_cols(f0_280_400_tidy, "F0", "280_400")
f0_400_600_tidy <- insert_cols(f0_400_600_tidy, "F0", "400_600")

suffix_235_415_tidy <- insert_cols(suffix_235_415_tidy, "suffix", "235_415")
suffix_415_550_tidy <- insert_cols(suffix_415_550_tidy, "suffix", "415_550")
suffix_550_680_tidy <- insert_cols(suffix_550_680_tidy, "suffix", "550_680")
suffix_680_800_tidy <- insert_cols(suffix_680_800_tidy, "suffix", "680_800")

# Change area and time period column order
f0_136_280_tidy <- f0_136_280_tidy %>% relocate(area, .after = L1_L2) %>% relocate(time_period, .after = area)
f0_280_400_tidy <- f0_280_400_tidy %>% relocate(area, .after = L1_L2) %>% relocate(time_period, .after = area)
f0_400_600_tidy <- f0_400_600_tidy %>% relocate(area, .after = L1_L2) %>% relocate(time_period, .after = area)

suffix_235_415_tidy <- suffix_235_415_tidy %>% relocate(area, .after = L1_L2) %>% relocate(time_period, .after = area)
suffix_415_550_tidy <- suffix_415_550_tidy %>% relocate(area, .after = L1_L2) %>% relocate(time_period, .after = area)
suffix_550_680_tidy <- suffix_550_680_tidy %>% relocate(area, .after = L1_L2) %>% relocate(time_period, .after = area)
suffix_680_800_tidy <- suffix_680_800_tidy %>% relocate(area, .after = L1_L2) %>% relocate(time_period, .after = area)

# Combined F0 dataset
ERP_F0_tidy <- bind_rows(
  f0_136_280_tidy,
  f0_280_400_tidy,
  f0_400_600_tidy
)

write_csv(ERP_F0_tidy, "all_erps/ERP_F0_tidy.csv")

# Combined suffix dataset
ERP_suffix_tidy <- bind_rows(
  suffix_235_415_tidy,
  suffix_415_550_tidy,
  suffix_550_680_tidy,
  suffix_680_800_tidy
)

write_csv(ERP_suffix_tidy, "all_erps/ERP_suffix_tidy.csv")

# tidy combined data
ERP_F0_long <- ERP_F0_tidy %>%
  pivot_longer(
    cols = starts_with("Fp1-"):last_col(),  # adjust as needed
    names_to = "channel_condition",
    values_to = "amplitude"
  ) %>%
  separate(channel_condition, into = c("channel", "descriptor"), sep = "-") %>%
  separate(descriptor, into = c("measure", "area", "stress", "match"), sep = "_") %>%
  select(Participant_ID, L1_L2, area, time_period, channel, stress, match, amplitude)

ERP_suffix_long <- ERP_suffix_tidy %>%
  pivot_longer(
    cols = starts_with("Fp1-"):last_col(),  # adjust as needed
    names_to = "channel_condition",
    values_to = "amplitude"
  ) %>%
  separate(channel_condition, into = c("channel", "descriptor"), sep = "-") %>%
  separate(descriptor, into = c("measure", "area", "stress", "match"), sep = "_") %>%
  select(Participant_ID, L1_L2, area, time_period, channel, stress, match, amplitude)


# Save the tidy datasets
write_csv(ERP_F0_long, "tidy_data_all_erps/ERP_F0_long.csv")
write_csv(ERP_suffix_long, "tidy_data_all_erps/ERP_suffix_long.csv")

### Stats
glimpse(ERP_F0_long)
glimpse(ERP_suffix_long)

## F0 Annovas
# 136 to 280 time window
# Filter data for the first PrAN time window and desired ROI
ERP_F0_long <- ERP_F0_long %>%
  mutate(Subject_ID = paste(Participant_ID, L1_L2, sep = "_"))

f0_136_data <- ERP_F0_long %>%
  filter(time_period == "136_280",
         channel %in% c("F3", "FZ", "FC1", "CZ"))

# Run repeated-measures ANOVA
aov_pran_136 <- aov_ez(
  id = "Subject_ID",
  dv = "amplitude",
  within = c("stress", "match"),
  between = "L1_L2",
  data = f0_136_data,
  type = 3,
  return = "afex_aov"
)

# View results
summary(aov_pran_136)

# Post-hoc comparisons
emmeans(aov_pran_136, ~ stress * match * L1_L2)

# 280-400 ms time window
# Subset the data for the 280–400 ms time window
f0_280_data <- ERP_F0_long %>%
  filter(time_period == "280_400",
         channel %in% c("F3", "FZ", "FC1", "CZ"))

# Create unique subject ID again (if not already created)
f0_280_data <- f0_280_data %>%
  mutate(Subject_ID = paste(Participant_ID, L1_L2, sep = "_"))

# Run repeated-measures ANOVA
aov_pran_280 <- aov_ez(
  id = "Subject_ID",
  dv = "amplitude",
  within = c("stress", "match"),
  between = "L1_L2",
  data = f0_280_data,
  type = 3,
  return = "afex_aov"
)

# Show ANOVA summary
summary(aov_pran_280)

# Optional: get estimated marginal means
emmeans(aov_pran_280, ~ stress * match * L1_L2)

# 400 to 600 ms 
# Subset the data for the 400–600 ms time window
f0_400_data <- ERP_F0_long %>%
  filter(time_period == "400_600",
         channel %in% c("F3", "FZ", "FC1", "CZ")) %>%
  mutate(Subject_ID = paste(Participant_ID, L1_L2, sep = "_"))

# Run repeated-measures ANOVA
aov_pran_400 <- aov_ez(
  id = "Subject_ID",
  dv = "amplitude",
  within = c("stress", "match"),
  between = "L1_L2",
  data = f0_400_data,
  type = 3,
  return = "afex_aov"
)

# View summary
summary(aov_pran_400)

# Optional: Estimated marginal means
emmeans(aov_pran_400, ~ stress * match * L1_L2)

## Suffix Annovas
# Create Subject_ID
ERP_suffix_long <- ERP_suffix_long %>%
  mutate(Subject_ID = paste(Participant_ID, L1_L2, sep = "_"))

# 235–415 ms
suffix_235_data <- ERP_suffix_long %>%
  filter(time_period == "235_415",
         channel %in% c("F7", "F3", "FC5"))

aov_suffix_235 <- aov_ez(
  id = "Subject_ID",
  dv = "amplitude",
  within = c("stress", "match"),
  between = "L1_L2",
  data = suffix_235_data,
  type = 3,
  return = "afex_aov"
)

# View summary
summary(aov_suffix_235)

# Post-hoc comparisons
emmeans(aov_suffix_235, ~ stress * match * L1_L2)

# 415–550 ms
suffix_415_data <- ERP_suffix_long %>%
  filter(time_period == "415_550",
         channel %in% c("F7", "F3", "FC5"))

aov_suffix_415 <- aov_ez(
  id = "Subject_ID",
  dv = "amplitude",
  within = c("stress", "match"),
  between = "L1_L2",
  data = suffix_415_data,
  type = 3,
  return = "afex_aov"
)

summary(aov_suffix_415)
emmeans(aov_suffix_415, ~ stress * match * L1_L2)

# 550–680 ms
suffix_550_data <- ERP_suffix_long %>%
  filter(time_period == "550_680",
         channel %in% c("F7", "F3", "FC5"))

aov_suffix_550 <- aov_ez(
  id = "Subject_ID",
  dv = "amplitude",
  within = c("stress", "match"),
  between = "L1_L2",
  data = suffix_550_data,
  type = 3,
  return = "afex_aov"
)

summary(aov_suffix_550)
emmeans(aov_suffix_550, ~ stress * match * L1_L2)

# 680-800 ms 
# Subset data for 680–800 ms window
suffix_680_data <- ERP_suffix_long %>%
  filter(time_period == "680_800",
         channel %in% c("F7", "F3", "FC5"))

aov_suffix_680 <- aov_ez(
  id = "Subject_ID",
  dv = "amplitude",
  within = c("stress", "match"),
  between = "L1_L2",
  data = suffix_680_data,
  type = 3,
  return = "afex_aov"
)

# Summary of ANOVA
summary(aov_suffix_680)

# Estimated marginal means
emmeans(aov_suffix_680, ~ stress * match * L1_L2)

## Topographic plots
# F0
names(ERP_F0_long)
unique(ERP_F0_long$channel)

# ==== F0 ANALYSIS ====
f0_diff <- ERP_F0_long %>%
  filter(channel %in% c("F3", "Fz", "FC1", "Cz"),
         time_period %in% c("136_280", "280_400", "400_600")) %>%
  group_by(channel, L1_L2, stress, time_period) %>%
  summarise(mean_amp = mean(amplitude), .groups = "drop") %>%
  pivot_wider(names_from = stress, values_from = mean_amp) %>%
  mutate(diff = stressed - unstressed)

electrode_coords <- tibble::tribble(
  ~electrode, ~x,   ~y,
  "F3",       -0.4,  0.5,
  "Fz",        0.0,  0.7,
  "FC1",       0.3,  0.4,
  "Cz",        0.0,  0.0
)

f0_diff_coords <- f0_diff %>%
  rename(electrode = channel) %>%
  left_join(electrode_coords, by = "electrode")

ggplot(f0_diff_coords, aes(x = x, y = y, fill = diff, z = diff)) +
  geom_topo(interp_limit = "head", chan_markers = "point", colour = "black") +
  scale_fill_distiller(palette = "RdBu", limits = c(-5, 5)) +
  coord_equal() +
  theme_void() +
  facet_grid(L1_L2 ~ time_period) +
  labs(
    title = "Stressed – Unstressed at F0 Onset",
    fill = expression(Delta~Amplitude~(mu*V))
  )

# Optional: View just L2 values
f0_diff_coords %>%
  filter(L1_L2 == "L2") %>%
  select(electrode, time_period, diff)

# Optional: Confirm counts
ERP_F0_long %>%
  filter(channel %in% c("F3", "Fz", "FC1", "Cz"),
         time_period %in% c("136_280", "280_400", "400_600"),
         L1_L2 == "L2") %>%
  count(channel, time_period)

# ==== SUFFIX ANALYSIS ====
suffix_diff <- ERP_suffix_long %>%
  filter(channel %in% c("F7", "F3", "FC5"),
         time_period %in% c("235_415", "415_550", "550_680", "680_800")) %>%
  group_by(channel, L1_L2, match, time_period) %>%
  summarise(mean_amp = mean(amplitude), .groups = "drop") %>%
  pivot_wider(names_from = match, values_from = mean_amp) %>%
  mutate(diff = mismatch - match)

electrode_coords_suffix <- tibble::tribble(
  ~electrode, ~x,   ~y,
  "F7",       -0.7,  0.3,
  "F3",       -0.4,  0.5,
  "FC5",      -0.5,  0.4
)

suffix_diff_coords <- suffix_diff %>%
  rename(electrode = channel) %>%
  left_join(electrode_coords_suffix, by = "electrode")

ggplot(suffix_diff_coords, aes(x = x, y = y, fill = diff, z = diff)) +
  geom_topo(interp_limit = "head", chan_markers = "point", colour = "black") +
  scale_fill_distiller(palette = "RdBu", limits = c(-5, 5)) +
  coord_equal() +
  theme_void() +
  facet_grid(L1_L2 ~ time_period) +
  labs(
    title = "Mismatch – Match at Suffix Onset",
    fill = expression(Delta~Amplitude~(mu*V))
  )

# Optional: View just L2 values
suffix_diff_coords %>%
  filter(L1_L2 == "L2") %>%
  select(electrode, time_period, diff)

# Optional: Confirm counts
ERP_suffix_long %>%
  filter(channel %in% c("F7", "F3", "FC5"),
         time_period %in% c("235_415", "415_550", "550_680", "680_800"),
         L1_L2 == "L2") %>%
  count(channel, time_period)