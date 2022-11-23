# Gene Set Enrichment Analysis


Download GSEA GUI app from: http://www.gsea-msigdb.org/gsea/downloads.jsp

Download `C5: ontology gene sets` GMT file from: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C5

Run `gsea.R` script to create input GCT file (also add appropriate header)

Select the following on the Run Gsea tab:
- Expression dataset: `rsem_GENE_shaped.gct`
- Gene sets database: `c5.all.v2022.1.Hs.symbols.gmt`
- Number of permutations: 1000
- Phenotype labels: `groups.cls`
- Collapse/Remap to gene symbols: Collapse
- Permutation type: phenotype
- Chip platform: `Human_Gene_Symbol_with_Remapping_MSigDB.v2022.1.Hs.chip`

