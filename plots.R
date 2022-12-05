library(ggplot2)
# devtools::install_github("nicolash2/ggdendroplot")
library(ggdendroplot)
library(dplyr)
library(readr)
library(readxl)

### HEAT MAP of EXPRESION w/ DENDROGRAMS

## Get Expression Stats (and filter to significant genes)
mecp2_data <- read_excel("RNA Seq data edgeRglm_GENE_Mecp2M-Mecp2WT.xlsx", sheet = "edgeR GLM gene") #%>% 
  # filter(`Mecp2M-Mecp2WT_FDR` < 0.5,
  #        `Mecp2M-Mecp2WT_Status` != "NS",
  #        `Mecp2M-Mecp2WT_PValue` <= 0.05
  #        )

## Get Entrez IDs

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         # dataset = "hsapiens_gene_ensembl",
                         dataset = "drerio_gene_ensembl",
                         host = "http://www.ensembl.org")

genes <- biomaRt::getBM(filters = "ensembl_gene_id",
                        attributes = c("ensembl_gene_id", "entrezgene_id"),
                        values = mecp2_data$Ensembl,
                        mart = mart)

mecp2_data <- mecp2_data %>% left_join(genes, by = c("Ensembl" = "ensembl_gene_id"))


## RSEM Clusters


rsem_z <- read_csv("rsem_GENE_z.csv") %>% filter(gene_id %in% mecp2_data$Ensembl)
  

rsem_z_mat <- rsem_z %>% select(!c("gene_id", "Symbol")) %>% as.matrix()
colnames(rsem_z_mat) <- colnames(rsem_z %>% select(!c("gene_id", "Symbol")))
rownames(rsem_z_mat) <- rsem_z$Symbol


## Get hierarchical clusters for rows and cols
row_clusters <- rsem_z_mat %>% dist() %>% hclust()
column_clusters <- rsem_z_mat %>% t() %>% dist() %>% hclust()

rsem_z_hm <- hmReady(rsem_z_mat, colclus=column_clusters, rowclus=row_clusters)

write_csv(rsem_z_hm, "rsem_GENE_z_clust.csv")

ggplot() + 
  geom_tile(data=rsem_z_hm, aes(x=y, y=x, fill=value)) +               ## heatmap
  scale_fill_gradientn(colors=hmGradient(), limits=c(-4,4)) +          ## options for heatmap
  geom_dendro(column_clusters, pointing="side") +
  #geom_dendro(column_clusters, ylim=c(16.5, 20)) +                     ## upper dendrogram
  #geom_dendro(row_clusters, xlim=c(8.5, 10), pointing="side") +        ## side dendrogram
  theme_hm()+                                                          ## design
  theme(axis.title=element_blank())                                    ## design




####################
## Enrichment Plots

library(clusterProfiler)
library(enrichplot)

activated_genes <- mecp2_data %>% 
  filter(`Mecp2M-Mecp2WT_Status` == "UP")

suppressed_genes <- mecp2_data %>% 
  filter(`Mecp2M-Mecp2WT_Status` == "DOWN")

activated_ego <- enrichGO(as.character(activated_genes$entrezgene_id),
                          OrgDb = "org.Dr.eg.db",
                          ont="ALL", readable=TRUE)

suppressed_ego <- enrichGO(as.character(suppressed_genes$entrezgene_id),
                           OrgDb = "org.Dr.eg.db",
                           ont="ALL", readable=TRUE)


activated_ego_plot <- goplot(activated_ego)

activated_ego_simp <- simplify(activated_ego)

activated_ego_simp_plot <- goplot(activated_ego_simp)

