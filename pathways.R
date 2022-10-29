library(pathview)
library(dplyr)
library(readr)


## DEMO: http://bioconductor.org/packages/release/bioc/vignettes/pathview/inst/doc/pathview.pdf
data(gse16873.d)
data(demo.paths)

head(demo.paths,3)


i <- 1
pv.out <- pathview(gene.data = gse16873.d, pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873", kegg.native = T)
list.files(pattern="hsa04110", full.names=T)



## Pathway Colorization
rsem_data_z <- read_csv("rsem_GENE_z.csv")
rsem_data_z_filtered <- rsem_data_z %>% dplyr::select(!c("gene_id", "Symbol"))
rownames(rsem_data_z_filtered) <- rsem_data_z$gene_id


## Example: Phototransduction pathway
pt_pv <- pathview(gene.data = rsem_data_z_filtered,
                 gene.idtype="ENSEMBL",
                 pathway.id = "04744",
                 species = "dre",
                 kegg.native = T)

## Example: Calcium Signaling pathway
ca_pv <- pathview(gene.data = rsem_data_z_filtered,
                  gene.idtype="ENSEMBL",
                  pathway.id = "04020",
                  species = "dre",
                  kegg.native = T)
