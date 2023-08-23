#!/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)
library(stringr)
print(paste0("Seurat Version: ", packageVersion("Seurat")))

library(argparser, quietly = TRUE)

parser <- arg_parser("Deconvolve visium spots in an RDS file and output cell fractions")

# Add command line arguments
parser <- add_argument(parser, "inputRDSFile", help = "RDS File to deconvolve and label")
parser <- add_argument(parser, "--atlas", help = "Atlas file in h5Seurat format", default = "KidneyAtlas_snCV3_20percent.h5Seurat")
parser <- add_argument(parser, "--rds", help = "Adjusted RDS file", default = "")
parser <- add_argument(parser, "--cell", help = "Calculated Cell Fractions", default = "")
parser <- add_argument(parser, "--spot", help = "Calculated Spot Coordinates", default = "")

# Parse the command line arguments
argv <- parse_args(parser)

spatial <- readRDS(argv$inputRDSFile)
print("Read RDS:")
print(names(spatial@assays))
DefaultAssay(spatial) <- "SCT"

kbrname <- argv$atlas
print(paste0("KBR Atlas Name: ", kbrname))
KBR <- LoadH5Seurat(kbrname, assays = c("counts", "scale.data"), tools = TRUE, images = FALSE)
print("Loaded H5Seurat")

Idents(KBR) <- KBR@meta.data$subclass.l2
print("Set Idents")
KBR <- subset(KBR, idents = "NA", invert = T)
print("Created Subset")

# There is a bug in UpdateSCTAssays. Currently, you can use as to convert integrated assay from Assay to SCTAssay
KBR <- UpdateSeuratObject(KBR)
print("Updated Seurat Object")
KBR[["RNA"]] <- as(object = KBR[["RNA"]], Class = "SCTAssay")
print("Retyped RNA")
DefaultAssay(KBR) <- "RNA"
print("Set DefaultAssay")

Idents(KBR) <- KBR@meta.data[["subclass.l2"]]
print("Set Idents")
anchors <- FindTransferAnchors(
    reference = KBR, query = spatial, normalization.method = "SCT",
    query.assay = "SCT", recompute.residuals = FALSE)
print("FindTransferAnchors")
predictions.assay <- TransferData(
    anchorset = anchors, refdata = KBR@meta.data[["subclass.l2"]],
    prediction.assay = TRUE,
    weight.reduction = spatial[["pca"]], dims = 1:30)
print("TransferData:")
print(head(predictions.assay[, 1:5]))
spatial[["predsubclassl2"]] <- predictions.assay

print("Set spatial")
df_pred <- predictions.assay@data
print("Set df_pred")
max_pred <- apply(df_pred, 2, function(x) max.col(t(x), "first"))
max_pred_val <- apply(df_pred, 2, function(x) max(t(x)))
max_pred <- as.data.frame(max_pred)
max_pred$Seurat_subset <- rownames(df_pred)[max_pred$max_pred]
max_pred$score <- max_pred_val
max_pred$Barcode <- rownames(max_pred)

spatial@meta.data$subclass.l2 <- max_pred$Seurat_subset
spatial@meta.data$subclass.l2_score <- max_pred$score
print("Update spatial")

Idents(KBR) <- KBR@meta.data[["subclass.l1"]]
print("Set Idents")
anchors <- FindTransferAnchors(
    reference = KBR, query = spatial, normalization.method = "SCT",
    query.assay = "SCT", recompute.residuals = FALSE)
print("FindTransferAnchors")
predictions.assay <- TransferData(
    anchorset = anchors, refdata = KBR@meta.data[["subclass.l1"]],
    prediction.assay = TRUE,
    weight.reduction = spatial[["pca"]], dims = 1:30)
print("TransferData")
print(head(predictions.assay[, 1:5]))
spatial[["predsubclassl1"]] <- predictions.assay
print(names(spatial@assays))

print("Set spatial")
df_pred <- predictions.assay@data
print("Set ds_pred")
max_pred <- apply(df_pred, 2, function(x) max.col(t(x), "first"))
max_pred_val <- apply(df_pred, 2, function(x) max(t(x)))
max_pred <- as.data.frame(max_pred)
max_pred$Seurat_subset <- rownames(df_pred)[max_pred$max_pred]
max_pred$score <- max_pred_val
max_pred$Barcode <- rownames(max_pred)

spatial@meta.data$subclass.l1 <- max_pred$Seurat_subset
spatial@meta.data$subclass.l1_score <- max_pred$score

print("Update spatial")
if (argv$rds != "") {
    saveRDS(spatial, argv$rds)

    print(paste0("Saved rds: ", argv$rds))
}

cell_type_fract <- GetAssayData(spatial@assays[["predsubclassl2"]])

# Normalizing so that columns sum to 1
cell_type_fract <- cell_type_fract[1:nrow(cell_type_fract) - 1, ]
cell_type_norm <- cell_type_fract / colSums(cell_type_fract)
cell_type_norm[is.na(cell_type_norm)] = 0

# Getting spot barcodes and coordinates
spot_coords <- spatial@images[["slice1"]]@coordinates

if (argv$cell != "") {
    write.csv(cell_type_norm, argv$cell)
    print(paste0("Saved cell csv: ", argv$cell))
}
if (argv$spot != "") {
    write.csv(spot_coords, argv$spot)
    print(paste0("Saved spot coordinates csv: ", argv$spot))
}
