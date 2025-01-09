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
# --- Load packages ----------------------------------------- #
# ----------------------------------------------------------- #

library(Seurat)
library(data.table)

# ----------------------------------------------------------- #
# --- Load data & parsing ----------------------------------- #
# ----------------------------------------------------------- #

path = "./"

## Read file
raw.counts = data.table::fread(paste0(path, "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"))
dim(raw.counts)
raw.counts[1:10, 1:10]

## Convert into a data.frame
raw.counts = as.data.frame(raw.counts)
class(raw.counts)
raw.counts[1:10, 1:10]

## Parsing gene annotation
rownames(raw.counts) = raw.counts[,1]
raw.counts = raw.counts[,-1]
dim(raw.counts)
raw.counts[1:10, 1:10]

## Check order
all.equal(colnames(raw.counts), rownames(metadata))

## Is present RBFOX2? YES
which(rownames(raw.counts)=="RBFOX2")
raw.RBFOX2 = raw.counts["RBFOX2",]

## Construct Seurat object
so = CreateSeuratObject(counts = as.matrix(raw.counts))
so = AddMetaData(so, metadata = metadata)
rm(raw.counts)

# ----------------------------------------------------------- #
# --- Quality control --------------------------------------- #
# ----------------------------------------------------------- #

## Percentage mito genes
so[["ratio.mito"]] <- PercentageFeatureSet(so, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "ratio.mito"), ncol = 3)

## Quality control filtering
data = subset(so, subset = nFeature_RNA >= 200 & nFeature_RNA <= 8000 & nCount_RNA >= 200 & ratio.mito < 10)


# ----------------------------------------------------------- #
# --- Normalised data --------------------------------------- #
# ----------------------------------------------------------- #

## Normalize data
so = NormalizeData(object = so)
p = so

## Filter data
so = subset(x = so, subset = site %in% c("Normal", "Tumor"))
dim(so)


