library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)


## Read in GO Data (From: https://current.geneontology.org/annotations/zfin.gaf.gz on https://zfin.org/downloads)
go_cols <- c("database_designation", "marker_id","gene_symbol",
             "qualifiers", "go_term_id", "reference_id", 
             "go_evidence_code", "inferred_from", "ontology",
             "marker_name", "marker_synonyms", "marker_type", "taxon", 
             "modification_date", "assigned_by", "annotation_extension", 
             "gene_product_form_id")

go_data <- read_tsv("../data/zfin.gaf", comment = "!", col_names = go_cols)

## Get Ensembl to Zfin ID lookup (from: https://zfin.org/downloads/ensembl_1_to_1.txt)
ens_to_zfin <- read_tsv("../data/ensembl_to_zfin.tsv", 
                        col_names = c("zfin_id", "so_id", 
                                      "symbol", "ensembl_id"))


## Read in GO Ontology Information (from: http://purl.obolibrary.org/obo/go.obo)
library(ontologyIndex)

go_obo <- get_OBO(
  "../data/go.obo",
  propagate_relationships = "is_a",
  extract_tags = "minimal"#,
  # merge_equivalent_terms = TRUE
)

go_obo_df <- data.frame(
  go_id = go_obo$id,
  go_name = go_obo$name
)

## Get Expression Stats
mecp2_data <- read_excel("../data/RNA Seq data edgeRglm_GENE_Mecp2M-Mecp2WT.xlsx", sheet = "edgeR GLM gene")

## Join Altogether
go_data_join <- go_data %>% 
  inner_join(ens_to_zfin, by = c("marker_id" = "zfin_id")) %>% 
  select(ensembl_id, symbol, qualifiers, go_term_id) %>% 
  inner_join(go_obo_df, by = c("go_term_id" = "go_id")) %>% 
  inner_join(mecp2_data, by = c("ensembl_id" = "Ensembl")) %>% 
  select(!c("Symbol", "Description")) %>% 
  rename(log_fc = `Mecp2M-Mecp2WT_logFC`,
         log_cpm = `Mecp2M-Mecp2WT_logCPM`,
         lr = `Mecp2M-Mecp2WT_LR`,
         p_value = `Mecp2M-Mecp2WT_PValue`,
         fdr = `Mecp2M-Mecp2WT_FDR`,
         status = `Mecp2M-Mecp2WT_Status`) %>% 
  mutate(expression_change = if_else(status == "NS", FALSE, TRUE),
         expression_up = if_else(status == "UP", TRUE, FALSE),
         expression_down = if_else(status == "DOWN", TRUE, FALSE))

## Make Summary Matrix
go_changes <- go_data_join %>% 
  group_by(go_term_id, go_name) %>% 
  summarise(n_genes = n(),
            n_changed = sum(expression_change),
            n_up = sum(expression_up),
            n_down = sum(expression_down)) %>% 
  mutate(n_same = n_genes - n_changed,
         pct_changed = n_changed / n_genes,
         pct_same = 100 - pct_changed,
         pct_up = n_up / n_genes,
         pct_down = n_down / n_genes) %>% 
  # filter(n_genes >= 300, n_genes <= 5000)
  # filter(n_genes >= 50, n_genes <= 1000)
  filter(n_genes >= 10)

write_csv(go_changes, "../data/go_changes_gte10.csv")

## Pivot Summary Data Frame for Viz
go_changes_pvt <- go_changes %>% 
  pivot_longer(!c("go_term_id", "go_name"),
               names_to = "change", values_to = "value")

## Plots
library(ggplot2)

### Stacked Counts
ggplot(data=go_changes_pvt %>% filter(change %in% c("n_up", "n_down", "n_same")), 
       aes(x=reorder(go_term_id, -value),
           y=value,
           fill=factor(change, levels=c("n_up", "n_same", "n_down")))
       ) +
  geom_bar(stat="identity", color="black") +
  scale_fill_discrete(name = "", labels = c('Up', "Not Significant", 'Down')) +
  scale_x_discrete(name = "", labels = function(x) str_wrap(x, width = 30)) +
  scale_y_continuous(name = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

