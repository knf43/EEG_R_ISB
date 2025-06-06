# Libraries 
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

# Load data
safe_read_csv <- purrr::possibly(read_csv, otherwise = NULL)

# Function to read and label files
read_psychopy_data <- function(version, group) {
  folder <- file.path("Psychopy", version, group)
  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  
  map_dfr(files, function(file) {
    df <- safe_read_csv(file)
    if (!is.null(df)) {
      df %>%
        mutate(
          file_name = basename(file),
          version = version,
          group = group
        )
    } else {
      message(paste("Skipping unreadable file:", file))
      NULL
    }
  })
}

# Combine all
df_all <- bind_rows(
  read_psychopy_data("V1", "L1"),
  read_psychopy_data("V1", "L2"),
  read_psychopy_data("V2", "L1"),
  read_psychopy_data("V2", "L2")
)

# Check result
df_all %>% count(version, group)

str(df_all)

df_all <- df_all %>%
  mutate(
    response_key = experimental_trials.key_resp_2.keys,
    correct = case_when(
      # Version 1
      version == "V1" & corrAns == "num_1" & response_key == "left"  ~ TRUE,
      version == "V1" & corrAns == "num_2" & response_key == "right" ~ TRUE,
      
      # Version 2
      version == "V2" & corrAns == "num_1" & response_key == "right" ~ TRUE,
      version == "V2" & corrAns == "num_2" & response_key == "left"  ~ TRUE,
      
      # Incorrect if response was made but doesn’t match the correct mapping
      corrAns %in% c("num_1", "num_2") & !is.na(response_key) ~ FALSE,
      
      # Otherwise (e.g., missing response), NA
      TRUE ~ NA
    )
  )

df_clean <- df_all %>%
  filter(
    !is.na(correct),
    !(sentence %in% c("p1", "p2", "p3", "p4"))
  )

df_clean <- df_clean %>%
  mutate(
    is_correct = case_when(
      version == "V1" & present_past == "present" & response_key == "left"  ~ TRUE,
      version == "V1" & present_past == "past"    & response_key == "right" ~ TRUE,
      version == "V2" & present_past == "present" & response_key == "right" ~ TRUE,
      version == "V2" & present_past == "past"    & response_key == "left"  ~ TRUE,
      TRUE ~ FALSE
    )
  )

glimpse(df_clean)

# Save the cleaned dataset as CSV inside that folder
write_csv(df_clean, here::here("psychopy_clean_data", "psychopy_clean.csv"))

df_clean <- readr::read_csv(here::here("psychopy_clean_data", "psychopy_clean.csv"))

# Analyze data
accuracy_summary <- df_clean %>%
  group_by(group, stress_unstress, suffix) %>%
  summarise(
    accuracy = mean(correct, na.rm = TRUE),
    sd = sd(correct, na.rm = TRUE),
    .groups = "drop"
  )

# RT
rt_summary <- df_clean %>%
  group_by(group, stress_unstress, suffix, correct) %>%
  summarise(
    mean_rt = mean(key_resp_2.rt, na.rm = TRUE),
    sd = sd(key_resp_2.rt, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
# Accuracy Plot
accuracy_by_subj <- df_clean %>%
  filter(!is.na(correct)) %>%
  group_by(participant, group, stress_unstress, suffix) %>%
  summarise(acc = mean(correct), .groups = "drop")

ggplot(accuracy_by_subj, aes(x = stress_unstress, y = acc, fill = suffix)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1) +
  facet_wrap(~ group) +
  labs(
    title = "Accuracy by Stress and Suffix",
    x = "Stress Pattern",
    y = "Accuracy (Proportion Correct)"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal()



# RT Plot
df_clean %>%
  filter(
    !is.na(experimental_trials.key_resp_2.rt),
    experimental_trials.key_resp_2.rt >= 0.3,
    experimental_trials.key_resp_2.rt <= 10
  ) %>%
  ggplot(aes(x = stress_unstress, y = experimental_trials.key_resp_2.rt, fill = suffix)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1) +
  facet_wrap(~ group) +
  labs(
    title = "Reaction Times by Stress and Suffix",
    x = "Stress Pattern",
    y = "Reaction Time (s)"
  ) +
  theme_minimal()


