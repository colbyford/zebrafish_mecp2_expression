gR1 = pathview:::parseKGML2Graph2("../data/dre00100.xml", genes = F, 
                                  expand = F, split.group = F)

node.data = node.info(gR1)
mapped.gnodes = rownames(sb_pv_W$plot.data.gene)
node.data$labels[mapped.gnodes] = sb_pv_W$plot.data.gene$labels

mapped.cnodes = rownames(sb_pv_W$cplot.data.cpd)
node.data$labels[mapped.cnodes] = sb_pv_W$cplot.data.cpd$labels


library(igraph)
library(visNetwork)
gR1_ig <- graph_from_graphnel(gR1, name = TRUE, weight = TRUE, unlist.attrs = TRUE)

bind_rows(sb_pv_W$plot.data.gene, sb_pv_W$plot.data.cpd)

gR1_vn_nodes <- toVisNetworkData(gR1_ig)$nodes %>% 
  left_join(bind_rows(sb_pv_W$plot.data.gene, sb_pv_W$plot.data.cpd) %>% 
              tibble::rownames_to_column(),
            by=c("id" = "rowname")) %>% 
  select(!c(mol.data, mol.col)) %>% 
  na.omit() %>%
  mutate(label = labels)


gR1_vn_edges <- toVisNetworkData(gR1_ig)$edges

visNetwork(nodes = gR1_vn_nodes, edges = gR1_vn_edges)



library(KEGGgraph)
dre00100_nel <- parseKGML2Graph("../data/dre00100.xml", genesOnly=FALSE)
# dre00100_df <- parseKGML2DataFrame("dre00100.xml", reactions=TRUE)

dre00100_ig <- graph_from_graphnel(dre00100_nel, name = TRUE, weight = TRUE, unlist.attrs = TRUE)

dre00100_vn_nodes <- toVisNetworkData(dre00100_ig)$nodes %>% 
  left_join(
    bind_rows(
      sb_pv_W$plot.data.gene %>% mutate(kegg.names = paste0("dre:", kegg.names)),
      sb_pv_W$plot.data.cpd %>% mutate(kegg.names = paste0("cpd:", kegg.names)) 
      ),
            by=c("id" = "kegg.names")) %>% 
  select(!c(x, y)) %>% 
  mutate(label = labels) %>% 
  distinct()

dre00100_vn_edges <- toVisNetworkData(dre00100_ig)$edges

visNetwork(nodes = dre00100_vn_nodes,
           edges = dre00100_vn_edges) %>% 
  visEdges(arrows = "to")
