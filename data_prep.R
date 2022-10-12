library(readr)
library(dplyr)
library(tidyr)
library(stringr)


## Pivot RSEM data
rsem_data <- read_csv("rsem_GENE.csv")

rsem_pvt <- rsem_data %>% 
  pivot_longer(!c("gene_id",	"Symbol"),
               names_to = "sample_id",
                values_to = "expression") %>% 
  filter(!str_detect(sample_id, "ALX|alx")) %>% 
  mutate(group = substr(sample_id, 9, 9))
  
unique(rsem_pvt$sample_id)

write_csv(rsem_pvt, "rsem_GENE_pvt.csv")
