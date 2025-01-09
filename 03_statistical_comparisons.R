# ----------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------- #
# Manuscript: RBFOX2 Promotes Hepatocellular Carcinoma via regulation of cationic amino acid transport
#
# Code author: Irene Soler Sáez (isoler@cipf.es)
#
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ----------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------- #
# --- Load packages & data ---------------------------------- #
# ----------------------------------------------------------- #

library(Seurat)
library(data.table)

## Local functions

dge.local = function(cell.type, test.dge, genes, latent.dge = NULL){
  
  # Subset cell type
  so.cell = so[,so$celltype==cell.type]
  
  # Test
  ##Ident 2 as reference
  dge = FindMarkers(object = so.cell[genes,], assay = "RNA", slot = "data", ident.2 = "Normal", ident.1 = "Tumor",  test.use = test.dge, latent.vars = latent.dge)
  
  # Parse results
  dge = dge
  dge$cell.type = cell.type
  dge$Gene = rownames(dge)
  rownames(dge) = NULL
  
  # Return result
  return(dge)
  
}


## Load normalised data
so = readRDS(so)


# ----------------------------------------------------------- #
# --- Differential expression analysis ---------------------- #
# ----------------------------------------------------------- #

## Gene list
genes.core = c("SNRNP70", "SNRPA", "SNRPC", # U1
               "SF3A1", "SF3A2", "SF3A3", "SNRPA1", "SNRPB2", # U2
               "PPIH", "PRPF3", "PRPF4", "PRPF31", "SNU13", # U4/U6
               "CD2BP2", "DDX23", "EFTUD2", "PRPF6", "PRPF8", "SNRNP40", "SNRNP200", "TSSC4", "TXNL4A") # U5
table(rownames(so) %in% genes.core)

alt.splicing = c("PTBP2", "SCAF8", "MBNL2", "RBPMS", "BICC1", "RBFOX2", "NOVA1", "RBM6", "RBM47", "RALYL", "RBMS3", "RBMS1", "RBM33", "MBNL1", "A1CF", "SREK1", "NOVA2", "TRA2B", "SF3B1", "SRSF1", "SRSF2", "SRSF3", "SRSF5", "SRSF7", "SLU7", "NONO", "HNRNPA1", "SRSF10", "ESRP2", "PRPF8")
table(rownames(so) %in% alt.splicing)

targets = c("SLC7A2", "SCARB1", "NUMB", "OGDH")
table(rownames(so) %in% targets)

all.genes = c(genes.core, alt.splicing, targets)


## Establish as index variable to test
Idents(so) = "site"

# Test linear regression by cell type adding as latent variable the donor
b = dge.local(cell.type = "B", test.dge = "LR", latent.dge = "patient", genes = genes2explore)
endo = dge.local(cell.type = "Endothelial", test.dge = "LR", latent.dge = "patient", genes = genes2explore)
fibro = dge.local(cell.type = "Fibroblast", test.dge = "LR", latent.dge = "patient", genes = genes2explore)
hepa = dge.local(cell.type = "Hepatocyte", test.dge = "LR", latent.dge = "patient", genes = genes2explore)
mye = dge.local(cell.type = "Myeloid", test.dge = "LR", latent.dge = "patient", genes = genes2explore)
n = dge.local(cell.type = "T/NK", test.dge = "LR", latent.dge = "patient", genes = genes2explore)

# Join results
all = rbind(b, endo, fibro, hepa, mye, n)

# Parse results
colnames(all) = c("p.value", "logFC", "percent.Tumor", "percent.Normal", "p.adj", "cell.type", "Gene")
all$p.adj = p.adjust(all$p.value, method = "BH")

all = all[,c("Gene", "cell.type", "percent.Normal", "percent.Tumor", "logFC", "p.value", "p.adj")]

# Generate table
DT::datatable(all, rownames = FALSE) %>% DT::formatSignif(c(5:7), digits = 3) %>% 
  DT::formatStyle("p.adj", backgroundColor = styleInterval(c(0.05), c("#cce3de", "NA"))) %>% 
  DT::formatStyle("logFC", backgroundColor = styleInterval(c(0), c("#ffccd5", "#faedcd")))

write.table(x = all, file = "dge.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# ----------------------------------------------------------- #
# --- Cell proportion analysis ------------------------------ #
# ----------------------------------------------------------- #

## Total cells
cells = data.frame(table(so$celltype, so$site))
colnames(cells) = c("Cells", "Site", "Freq")

cells.normal = cells[cells$Site=="Normal",]
cells.tumor = cells[cells$Site=="Tumor",]

## Porportions
all = all[all$Gene=="RBFOX2",]
abundances = all[,c("cell.type", "percent.Normal", "percent.Tumor")]
abundances.all = merge(abundances, cells.normal, by.x = "cell.type", by.y = "Cells")
colnames(abundances.all) = gsub(pattern = "Freq", replacement = "N.normal", x = colnames(abundances.all))

abundances.all = merge(abundances.all, cells.tumor, by.x = "cell.type", by.y = "Cells")
colnames(abundances.all) = gsub(pattern = "Freq", replacement = "N.tumor", x = colnames(abundances.all))

## Number of cells that express the gene
abundances.all$N.normal.test = abundances.all$percent.Normal * abundances.all$N.normal
abundances.all$N.tumor.test = abundances.all$percent.Tumor * abundances.all$N.tumor

## Prop test
p = prop.test(x = c(abundances.all$N.normal.test, abundances.all$N.tumor.test),
              n = c(abundances.all$N.normal, abundances.all$N.tumor), 
              conf.level = 0.95)

## Prop test pairwise
df = data.frame(matrix(0, nrow = 0, ncol = 5))
colnames(df) = c("Cells", "percent.Normal", "percent.Tumor", "X-squared", "p.value")

for (i in 1:nrow(abundances.all)) {
  
  ## Select cell type
  cell = abundances.all[i, ]
  
  p = prop.test(x = c(cell$N.normal.test, cell$N.tumor.test),
                n = c(cell$N.normal, cell$N.tumor), 
                conf.level = 0.95)
  
  stat = p$statistic
  p.val = p$p.value
  
  df.cell = data.frame("Cells" = cell$cell.type, 
                       "percent.Normal" = cell$percent.Normal,
                       "percent.Tumor" = cell$percent.Tumor, 
                       "X-squared" = stat, 
                       "p.value" = p.val)
  
  df = rbind(df, df.cell)
}

df$p.adjusted = p.adjust(p = df$p.value, method = "BH")

DT::datatable(df, rownames = FALSE) %>% DT::formatSignif(c(4:6), digits = 3) %>%
  DT::formatStyle("p.adjusted", backgroundColor = styleInterval(c(0.05), c("#cce3de", "NA"))) 
