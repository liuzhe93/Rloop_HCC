library(infercnv)
library(Matrix)
library(dplyr)

OUTDIR <- "/mnt/zhangzheng_group/liuz-52/Test_R/01_Figure1_python"
IN_DIR <- file.path(OUTDIR, "inferCNV_input")
INFER_OUT <- file.path(OUTDIR, "inferCNV_output")
RES_OUT <- file.path(OUTDIR, "inferCNV_results")

dir.create(INFER_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(RES_OUT, showWarnings = FALSE, recursive = TRUE)

# 1. 读取 sparse counts
counts <- Matrix::readMM(file.path(IN_DIR, "inferCNV_counts_matrix.mtx"))

genes <- read.table(file.path(IN_DIR, "inferCNV_genes.txt"), stringsAsFactors = FALSE)[,1]
cells <- read.table(file.path(IN_DIR, "inferCNV_cells.txt"), stringsAsFactors = FALSE)[,1]

rownames(counts) <- genes
colnames(counts) <- cells

# 2. annotation
anno_file <- file.path(IN_DIR, "inferCNV_cell_annotation.txt")
anno <- read.table(anno_file, sep = "\t", stringsAsFactors = FALSE)
colnames(anno) <- c("cell", "celltype")

print(table(anno$celltype))

# 3. reference groups
ref_groups <- c("T cell", "B cell", "TAM", "CAF", "TEC")
ref_groups <- intersect(ref_groups, unique(anno$celltype))

cat("Reference groups used:\n")
print(ref_groups)

if (length(ref_groups) < 2) {
  stop("Reference groups are too few. Please check celltype names.")
}

# 4. gene order
gene_order_file <- file.path(OUTDIR, "hg38_gene_order_file.txt")

if (!file.exists(gene_order_file)) {
  stop("hg38_gene_order_file.txt not found.")
}

# 5. Create inferCNV object
counts <- as.matrix(counts)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = anno_file,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = ref_groups
)

# 6. Run inferCNV
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = INFER_OUT,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  num_threads = 8
)

# 7. CNV score
cnv_mat <- infercnv_obj@expr.data

cnv_score <- Matrix::colMeans(abs(cnv_mat - 1))
cnv_score <- as.numeric(cnv_score)
names(cnv_score) <- colnames(cnv_mat)

meta <- anno
rownames(meta) <- meta$cell

meta$inferCNV_score <- NA
common_cells <- intersect(rownames(meta), names(cnv_score))
meta[common_cells, "inferCNV_score"] <- cnv_score[common_cells]

# 8. threshold based on reference cells
ref_cells <- rownames(meta)[meta$celltype %in% ref_groups]
ref_cells <- intersect(ref_cells, names(cnv_score))

cnv_cutoff <- median(cnv_score[ref_cells], na.rm = TRUE) +
  3 * mad(cnv_score[ref_cells], na.rm = TRUE)

cat("CNV cutoff:\n")
print(cnv_cutoff)

meta$marker_malignant <- meta$celltype == "Malignant cell"
meta$inferCNV_positive <- meta$inferCNV_score > cnv_cutoff

meta$malignant_inferCNV <- "Non-malignant"

meta$malignant_inferCNV[
  meta$marker_malignant & !meta$inferCNV_positive
] <- "Marker-defined malignant only"

meta$malignant_inferCNV[
  meta$marker_malignant & meta$inferCNV_positive
] <- "High-confidence malignant"

meta$malignant_inferCNV[
  !meta$marker_malignant & meta$inferCNV_positive
] <- "inferCNV-positive non-marker cell"

meta$malignant_inferCNV <- factor(
  meta$malignant_inferCNV,
  levels = c(
    "Non-malignant",
    "Marker-defined malignant only",
    "inferCNV-positive non-marker cell",
    "High-confidence malignant"
  )
)

print(table(meta$malignant_inferCNV, useNA = "ifany"))

# 9. 保存结果
write.csv(
  meta,
  file = file.path(RES_OUT, "Figure1_marker_malignant_inferCNV_intersection_metadata.csv"),
  quote = FALSE
)

high_conf <- meta[meta$malignant_inferCNV == "High-confidence malignant", ]

write.table(
  rownames(high_conf),
  file = file.path(RES_OUT, "High_confidence_malignant_cells_marker_and_inferCNV.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

saveRDS(
  infercnv_obj,
  file = file.path(RES_OUT, "infercnv_obj_final.rds")
)

cat("Done.\n")
cat("Results saved to:\n")
cat(RES_OUT, "\n")
