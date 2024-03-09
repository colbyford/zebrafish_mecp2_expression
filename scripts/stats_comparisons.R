library(dplyr)
library(readxl)
library(stringr)

## Lightflash Data
lf_data <- read_xlsx("../data/mecp2_lightflash_kinematics_data.xlsx",
                  sheet = "summary-all turns") %>% 
  mutate(group = case_when(
      str_detect(File, "_Mu_") ~ "Mutant",
      TRUE ~ "Wild Type"))

### Latency
wilcox.test(
  (lf_data %>% filter(group == "Mutant"))$Lat,
  (lf_data %>% filter(group == "Wild Type"))$Lat,
  correct = TRUE)

### Angle
wilcox.test(
  (lf_data %>% filter(group == "Mutant"))$Ang,
  (lf_data %>% filter(group == "Wild Type"))$Ang,
  correct = TRUE)

### Duration
wilcox.test(
  (lf_data %>% filter(group == "Mutant"))$Dur,
  (lf_data %>% filter(group == "Wild Type"))$Dur,
  correct = TRUE)

### Max Angular Velocity
wilcox.test(
  (lf_data %>% filter(group == "Mutant"))$Mav,
  (lf_data %>% filter(group == "Wild Type"))$Mav,
  correct = TRUE)


## Darkflash
df_data <- read_xlsx("../data/mecp2_darkflash_kinematics_data.xlsx",
                     sheet = "Summary - O-bends") %>% 
  mutate(group = case_when(
    str_detect(File, "_Mu_") ~ "Mutant",
    TRUE ~ "Wild Type"))

### Latency
wilcox.test(
  (df_data %>% filter(group == "Mutant"))$Lat,
  (df_data %>% filter(group == "Wild Type"))$Lat,
  correct = TRUE)

### Angle
wilcox.test(
  (df_data %>% filter(group == "Mutant"))$Ang,
  (df_data %>% filter(group == "Wild Type"))$Ang,
  correct = TRUE)

### Duration
wilcox.test(
  (df_data %>% filter(group == "Mutant"))$Dur,
  (df_data %>% filter(group == "Wild Type"))$Dur,
  correct = TRUE)

### Max Angular Velocity
wilcox.test(
  (df_data %>% filter(group == "Mutant"))$Mav,
  (df_data %>% filter(group == "Wild Type"))$Mav,
  correct = TRUE)
