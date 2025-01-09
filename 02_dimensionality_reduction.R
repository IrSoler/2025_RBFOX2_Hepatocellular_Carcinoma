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
library(ggplot2)

so = readRDS(so)


# ----------------------------------------------------------- #
# --- High Variable Genes ----------------------------------- #
# ----------------------------------------------------------- #

Idents(so) = "sample"

# Split object by samples
so[["RNA"]] = split(so[["RNA"]], f = so$sample)

so = FindVariableFeatures(so, selection.method = "vst", mean.cutoff = c(0.125, 3), dispersion.cutoff = c(0.5, Inf))


# ----------------------------------------------------------- #
# --- Principal Component Analysis -------------------------- #
# ----------------------------------------------------------- #

all.genes = rownames(so)
so = ScaleData(so, features = all.genes)

so = RunPCA(so)


# ----------------------------------------------------------- #
# --- Canonical Correlation Analysis ------------------------ #
# ----------------------------------------------------------- #

# Integration
so = IntegrateLayers(object = so, method = CCAIntegration, assay = "RNA", orig.reduction = "pca", new.reduction = "integrated.cca",
                     verbose = FALSE)

# Re-join labels after integration
so[["RNA"]] = JoinLayers(so[["RNA"]])


# ----------------------------------------------------------- #
# --- tSNE -------------------------------------------------- #
# ----------------------------------------------------------- #

# tSNE
so = RunTSNE(so, dims = 1:30, reduction = "integrated.cca")

# Plot
DimPlot(so, reduction = "tsne", dims = c(2,1))

