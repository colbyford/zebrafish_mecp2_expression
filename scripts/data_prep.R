library(readr)
library(dplyr)
library(tidyr)
library(stringr)


## Pivot RSEM data
rsem_data <- read_csv("rsem_GENE.csv") %>% select(c("gene_id", "Symbol", starts_with("s4d-mecp")))

rsem_pvt <- rsem_data %>% 
  pivot_longer(!c("gene_id",	"Symbol"),
               names_to = "sample_id",
                values_to = "expression") %>% 
  mutate(group = substr(sample_id, 9, 9))
  
unique(rsem_pvt$sample_id)

write_csv(rsem_pvt, "rsem_GENE_pvt.csv")



## Z transform RSEM Count
# sample_row <- rsem_data %>% filter(gene_id == "ENSDARG00000000001") %>% select(!c("gene_id", "Symbol")) %>% 
#   .[1,] %>% unlist() %>% unname()
# 
# scale(sample_row)

rsem_data_z <- apply(rsem_data %>% select(!c("gene_id", "Symbol")), 1, scale) %>% t() %>% as.data.frame()

colnames(rsem_data_z) <- colnames(rsem_data %>% select(!c("gene_id", "Symbol")))

rsem_data_z <- data.frame(gene_id = rsem_data$gene_id, 
                          Symbol = rsem_data$Symbol,
                          rsem_data_z)

write_csv(rsem_data_z, "rsem_GENE_z.csv")

## Pivot Z data

rsem_z_pvt <- rsem_data_z %>% 
  pivot_longer(!c("gene_id",	"Symbol"),
               names_to = "sample_id",
               values_to = "expression") %>% 
  mutate(group = substr(sample_id, 9, 9))

unique(rsem_z_pvt$sample_id)

write_csv(rsem_z_pvt, "rsem_GENE_z_pvt.csv")

