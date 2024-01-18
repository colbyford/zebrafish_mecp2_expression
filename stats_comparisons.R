library(dplyr)
library(readxl)
library(stringr)

data <- read_xlsx("021219_mecp2_DF Obend responsiveness-kinematics_d5.xlsx",
                  sheet = "Obend kin stim 6 to 0") %>% 
  mutate(group = case_when(
      str_detect(File, "_Mu_") ~ "Mutant",
      TRUE ~ "Wild Type"))

data_group_a <- (data %>% filter(group == "Mutant"))$Lat
data_group_b <- (data %>% filter(group == "Wild Type"))$Lat

# Perform Mann-Whitney U test with correction for ties
result <- wilcox.test(data_group_a, data_group_b, correct = TRUE)

# Display the test result
print(result)
