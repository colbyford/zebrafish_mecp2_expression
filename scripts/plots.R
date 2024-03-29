library(ggplot2)
# devtools::install_github("nicolash2/ggdendroplot")
library(ggdendroplot)
library(ggrepel)
library(dplyr)
library(readr)
library(readxl)

## Get Expression Stats (and filter to significant genes)
# p-value < 0.05 and fold change > 2 or < -2.
mecp2_data <- read_excel("../data/RNA Seq data edgeRglm_GENE_Mecp2M-Mecp2WT.xlsx", sheet = "edgeR GLM gene") %>% 
  mutate(significance = case_when(
    (`Mecp2M-Mecp2WT_logFC` > 1 & `Mecp2M-Mecp2WT_PValue` <= 0.05) ~ "UP",
    (`Mecp2M-Mecp2WT_logFC` < -1 & `Mecp2M-Mecp2WT_PValue` <= 0.05) ~ "DOWN",
    TRUE ~ "NS"
  ))
  

mecp2_data_sig <- mecp2_data %>% 
  # filter(`Mecp2M-Mecp2WT_Status` != "NS")
  filter(significance != "NS")


## Get Entrez IDs

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         # dataset = "hsapiens_gene_ensembl",
                         dataset = "drerio_gene_ensembl",
                         host = "https://www.ensembl.org")

genes <- biomaRt::getBM(filters = "ensembl_gene_id",
                        attributes = c("ensembl_gene_id", "entrezgene_id"),
                        values = mecp2_data_sig$Ensembl,
                        useCache = FALSE,
                        mart = mart)

mecp2_data_sig <- mecp2_data_sig %>% left_join(genes, by = c("Ensembl" = "ensembl_gene_id"))

######################
### HEAT MAP of EXPRESION w/ DENDROGRAMS
## RSEM Clusters


rsem_z <- read_csv("../data/rsem_GENE_z.csv") %>% filter(gene_id %in% mecp2_data_sig$Ensembl)
  

rsem_z_mat <- rsem_z %>% dplyr::select(!c("gene_id", "Symbol")) %>% as.matrix()
colnames(rsem_z_mat) <- colnames(rsem_z %>% dplyr::select(!c("gene_id", "Symbol")))
rownames(rsem_z_mat) <- rsem_z$Symbol


## Get hierarchical clusters for rows and cols
row_clusters <- rsem_z_mat %>% dist() %>% hclust()
column_clusters <- rsem_z_mat %>% t() %>% dist() %>% hclust()

rsem_z_hm <- hmReady(rsem_z_mat, colclus=column_clusters, rowclus=row_clusters)

write_csv(rsem_z_hm, "../data/rsem_GENE_z_clust.csv")

ggplot() + 
  geom_tile(data=rsem_z_hm, aes(x=y, y=x, fill=value)) +               ## heatmap
  scale_fill_gradientn(colors=hmGradient(), limits=c(-4,4)) +          ## options for heatmap
  geom_dendro(column_clusters, pointing="side", xlim=c(0,-150)) +
  #geom_dendro(column_clusters, ylim=c(16.5, 20)) +                     ## upper dendrogram
  #geom_dendro(row_clusters, xlim=c(8.5, 10), pointing="side") +        ## side dendrogram
  theme_hm()+                                                          ## design
  theme(axis.title=element_blank())                                    ## design




####################
## Enrichment Plots

library(clusterProfiler)
library(enrichplot)

activated_genes <- mecp2_data_sig %>% 
  # filter(`Mecp2M-Mecp2WT_Status` == "UP") %>%
  filter(significance == "UP") %>%
  na.omit() %>% 
  dplyr::select(entrezgene_id)

suppressed_genes <- mecp2_data_sig %>% 
  # filter(`Mecp2M-Mecp2WT_Status` == "DOWN") %>%
  filter(significance == "DOWN") %>%
  na.omit() %>% 
  dplyr::select(entrezgene_id)

activated_ego <- enrichGO(as.character(activated_genes$entrezgene_id),
                          OrgDb = "org.Dr.eg.db",
                          ont="BP", readable=TRUE)

suppressed_ego <- enrichGO(as.character(suppressed_genes$entrezgene_id),
                           OrgDb = "org.Dr.eg.db",
                           ont="BP", readable=TRUE)


# activated_ego_plot <- goplot(activated_ego)
# activated_ego_simp <- simplify(activated_ego)
# activated_ego_simp_plot <- goplot(activated_ego_simp)


# emapplot(activated_ego %>% simplify(), showCategory = 20)
# emapplot(activated_ego, showCategory = 20)
emapplot(pairwise_termsim(activated_ego), showCategory = 20)
# barplot(activated_ego, showCategory = 20)


# suppressed_ego_plot <- goplot(suppressed_ego)
# suppressed_ego_simp <- simplify(suppressed_ego)
# suppressed_ego_simp_plot <- goplot(suppressed_ego_simp)


# emapplot(suppressed_ego %>% simplify(), showCategory = 20)
emapplot(pairwise_termsim(suppressed_ego %>% simplify()), showCategory = 20)

######################
### Volcano Plot (from: https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/)

cols <- c("UP" = "red", "DOWN" = "blue", "NS" = "grey") 
sizes <- c("UP" = 2, "DOWN" = 2, "NS" = 1) 
alphas <- c("UP" = 1, "DOWN" = 1, "NS" = 0.5)

top10up <- mecp2_data %>%
  filter(`Mecp2M-Mecp2WT_PValue` <= 0.05) %>% 
  top_n(10, `Mecp2M-Mecp2WT_logFC`)

top10down <- mecp2_data %>%
  filter(`Mecp2M-Mecp2WT_PValue` <= 0.05) %>% 
  top_n(-10, `Mecp2M-Mecp2WT_logFC`)

top10both <- rbind(
  top10up, top10down
)

mecp2_data %>%
  ggplot(aes(x = `Mecp2M-Mecp2WT_logFC`,
             y = -log10(`Mecp2M-Mecp2WT_PValue`),
             fill = significance,    
             size = significance,
             alpha = significance)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_text_repel(data = top10up, # Add labels last to appear as the top layer  
                   aes(label = Symbol),
                  size = 3,
                   force = 2,
                   nudge_y = 1) +
  geom_text_repel(data = top10down, # Add labels last to appear as the top layer  
                   aes(label = Symbol),
                  size = 3,
                   force = 2,
                   nudge_y = 1) +
  scale_fill_manual(values = cols, guide="none") + # Modify point colour
  scale_size_manual(values = sizes, guide="none") + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),       
                     limits = c(-10, 10)) +
  labs(x="log2(Fold Change)",
       y="-log10(p-Value)",
       fill="Status",
       size="Status") +
  theme(legend.position="none")


######################
## Scatter Plot of Activated and Suppressed Genes


as_results <- rbind(activated_ego@result %>% 
                      mutate(change = "activated"),
                    suppressed_ego@result %>% 
                      mutate(change = "suppressed"))

as_results$gene_ratio <- sapply(as_results$GeneRatio, function(x) eval(parse(text = x)))
as_results$bg_ratio <- sapply(as_results$BgRatio, function(x) eval(parse(text = x)))

as_results_top_n <- as_results %>% 
  # arrange(desc(gene_ratio)) %>% 
  arrange(desc(bg_ratio)) %>% 
  filter(p.adjust < 0.05) %>%
  # group_by(change) %>%
  # slice(1:20)
  slice(1:40)
  


ggplot(data = as_results_top_n,
       aes(
         # x = gene_ratio,
         x = bg_ratio,
           # y = Description,
           # y = reorder(Description, gene_ratio),
         y = reorder(Description, bg_ratio),
           color = p.adjust,
           shape = change,
           size = Count)) +
  geom_point() +
  # facet_wrap(~change,
  #            # ncol = 1, nrow = 2,
  #            ncol = 2, nrow = 1,
  #            scales = "free"
  #            ) +
  labs(y='Gene Set Description',
       # x='Gene Ratio',
       x='Ratio of Annotated Genes\n(BgRatio)',
       shape = "Change",
       color = "Adj. p-value") +
  scale_color_gradient(low = "blue",
                       high = "red") +
  scale_shape_manual(
    # values=c("▲", "▼")
    values=c(17, 25)
    )



# act_scatter <- ggplot(data = as_results_top_n %>%
#                         filter(change == "activated") %>% 
#                         slice(1:10),
#                       aes(
#                         # x = gene_ratio,
#                         x = bg_ratio,
#                         # y = Description,
#                         y = reorder(Description, bg_ratio),
#                         color = p.adjust,
#                         size = Count)) +
#   geom_point() +
#   # facet_wrap(~change) +
#   labs(y='Description',
#        # x='Gene Ratio',
#        x='Ratio of Annotated Genes (BgRatio)',
#        color = "Adj. p-value") +
#   scale_color_gradient(low = "blue",
#                        high = "red") + 
#   scale_size_discrete(range = c(1,8))
# 
# 
# sup_scatter <- ggplot(data = as_results_top_n %>%
#                         filter(change == "suppressed") %>% 
#                         slice(1:10),
#                       aes(
#                         # x = gene_ratio,
#                         x = bg_ratio,
#                         # y = Description,
#                         y = reorder(Description, bg_ratio),
#                         color = p.adjust,
#                         size = Count)) +
#   geom_point() +
#   # facet_wrap(~change) +
#   labs(y='Description',
#        # x='Gene Ratio',
#        x='Ratio of Annotated Genes (BgRatio)',
#        color = "Adj. p-value") +
#   scale_color_gradient(low = "blue",
#                        high = "red") + 
#   scale_size(breaks = c(1:8), range = c(1,8))
# 
# 
# ggarrange(act_scatter, sup_scatter, 
#           labels = c("A. Activated", "B. Suppressed"),
#           # common.legend = TRUE,
#           nrow = 1)
