library(pathview)
library(dplyr)
library(readr)


## DEMO: http://bioconductor.org/packages/release/bioc/vignettes/pathview/inst/doc/pathview.pdf
# data(gse16873.d)
# data(demo.paths)
# 
# head(demo.paths,3)
# 
# 
# i <- 1
# pv.out <- pathview(gene.data = gse16873.d, pathway.id = demo.paths$sel.paths[i],
#                    species = "hsa", out.suffix = "gse16873", kegg.native = T)
# list.files(pattern="hsa04110", full.names=T)



## Pathway Colorization
rsem_data_z <- read_csv("rsem_GENE_z.csv")
rsem_data_z_filtered <- rsem_data_z %>% dplyr::select(!c("gene_id", "Symbol"))
rownames(rsem_data_z_filtered) <- rsem_data_z$gene_id

## Make Wild Type and Mutant Subsets
rsem_data_z_filtered_M <- rsem_data_z_filtered %>% select(contains("s4d.mecpM"))
rownames(rsem_data_z_filtered_M) <- rsem_data_z$gene_id

rsem_data_z_filtered_W <- rsem_data_z_filtered %>% select(contains("s4d.mecpW"))
rownames(rsem_data_z_filtered_W) <- rsem_data_z$gene_id

## Example: Phototransduction pathway
pt_pv_M <- pathview(gene.data = rsem_data_z_filtered_M,
                    gene.idtype="ENSEMBL",
                    pathway.id = "04744",
                    species = "dre",
                    out.suffix = "_M",
                    kegg.native = T)

pt_pv_W <- pathview(gene.data = rsem_data_z_filtered_W,
                    gene.idtype="ENSEMBL",
                    pathway.id = "04744",
                    species = "dre",
                    out.suffix = "_W",
                    kegg.native = T)

pt_pv_M_df <- pt_pv_M$plot.data.gene %>%
  select(labels, starts_with("s4d")) %>% 
  select(!ends_with(".col")) %>% 
  distinct()


pt_pv_W_df <- pt_pv_W$plot.data.gene %>%
  select(labels, starts_with("s4d")) %>% 
  select(!ends_with(".col")) %>% 
  distinct()


pt_pv_df <- pt_pv_W_df %>% 
  left_join(pt_pv_M_df) %>% 
  rename(Symbol = labels)

write_csv(pt_pv_df, "dre04744_pathview_values.csv")


rsem_data_z_filtered_dre04744 <- rsem_data_z %>% filter(Symbol %in% pt_pv_df$Symbol)

## Example: Calcium Signaling pathway

ca_pv_W <- pathview(gene.data = rsem_data_z_filtered_W,
                    gene.idtype="ENSEMBL",
                    pathway.id = "04020",
                    species = "dre",
                    out.suffix = "_W",
                    kegg.native = T)

ca_pv_M <- pathview(gene.data = rsem_data_z_filtered_M,
                    gene.idtype="ENSEMBL",
                    pathway.id = "04020",
                    species = "dre",
                    out.suffix = "_M",
                    kegg.native = T)


## Steroid hormone biosynthesis (Cholesterol Metabolism)

ch_pv_W <- pathview(gene.data = rsem_data_z_filtered_W,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00140",
                    species = "dre",
                    out.suffix = "_W",
                    low = list(gene = "magenta"),
                    high = list(gene = "green"),
                    kegg.native = T)

ch_pv_M <- pathview(gene.data = rsem_data_z_filtered_M,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00140",
                    species = "dre",
                    out.suffix = "_M",
                    low = list(gene = "magenta"),
                    high = list(gene = "green"),
                    kegg.native = T)

## Steroid biosynthesis (includes Cholesterol Synthesis)

sb_pv_W <- pathview(gene.data = rsem_data_z_filtered_W,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00100",
                    species = "dre",
                    out.suffix = "_W",
                    low = list(gene = "magenta"),
                    high = list(gene = "green"),
                    kegg.native = T)

sb_pv_M <- pathview(gene.data = rsem_data_z_filtered_M,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00100",
                    species = "dre",
                    out.suffix = "_M",
                    low = list(gene = "magenta"),
                    high = list(gene = "green"),
                    kegg.native = T)


## Mevalonate pathway (Terpenoid backbone biosynthesis)

mt_pv_W <- pathview(gene.data = rsem_data_z_filtered_W,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00900",
                    species = "dre",
                    out.suffix = "_W",
                    low = list(gene = "magenta"),
                    high = list(gene = "green"),
                    kegg.native = T)

mt_pv_M <- pathview(gene.data = rsem_data_z_filtered_M,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00900",
                    species = "dre",
                    out.suffix = "_M",
                    low = list(gene = "magenta"),
                    high = list(gene = "green"),
                    kegg.native = T)


######################
## Glycolysis --> ATP

## Glycolysis/Glucogenesis pathway

gg_pv_W <- pathview(gene.data = rsem_data_z_filtered_W,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00010",
                    species = "dre",
                    out.suffix = "_W",
                    kegg.native = T)

gg_pv_M <- pathview(gene.data = rsem_data_z_filtered_M,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00010",
                    species = "dre",
                    out.suffix = "_M",
                    kegg.native = T)

## Citrate Cycle

cc_pv_W <- pathview(gene.data = rsem_data_z_filtered_W,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00020",
                    species = "dre",
                    out.suffix = "_W",
                    kegg.native = T)

cc_pv_M <- pathview(gene.data = rsem_data_z_filtered_M,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00020",
                    species = "dre",
                    out.suffix = "_M",
                    kegg.native = T)


## Oxidative Phosphorylation

op_pv_W <- pathview(gene.data = rsem_data_z_filtered_W,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00190",
                    species = "dre",
                    out.suffix = "_W",
                    kegg.native = T)

op_pv_M <- pathview(gene.data = rsem_data_z_filtered_M,
                    gene.idtype="ENSEMBL",
                    pathway.id = "00190",
                    species = "dre",
                    out.suffix = "_M",
                    kegg.native = T)
