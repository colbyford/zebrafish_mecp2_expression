## Prepare GCT File of Expression Values
# BiocManager::install("ArrayTools")
library(ArrayTools)
library(dplyr)
library(tibble)
library(readr)

## Read in Expression Data

data <- read_csv("../rsem_GENE.csv") %>%
  select(c("Symbol", "gene_id", starts_with("s4d-mecp"))) %>% 
  rename(NAME = Symbol, Description = gene_id) # %>% 
  # column_to_rownames("NAME")
  
write_tsv(data, "rsem_GENE_shaped.gct")

## Read in Group Info

# groups <- read_tsv("groups.txt") %>% column_to_rownames("X1")
# 
# ## Create Expression Set
# 
# eSet <- createExpressionSet(pData = groups,
#                             exprs = data,
#                             annotation = "hugene10sttranscriptcluster")
# 
#   
# ## output GCT File
# 
# output.gct(data, filename = "expresion.gct")
