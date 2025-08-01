options(timeout = 300)  # Set timeout to 5 minutes
# Vectors
patient_ids <- c("GSM6592048", "GSM6592054", "GSM6592056", "GSM6592057", "GSM6592062",
                            "GSM6592049", "GSM6592050", "GSM6592051", "GSM6592052", "GSM6592053",
                            "GSM6592055", "GSM6592058")


sample_ids <- c("M1_CL-like3", "M7_CL-like4", "M9_CL-like3", "M10_CL-like4", "M16_MBC_NST",
                "M2_CL-like1", "M3_Unstable1", "M4_Unstable2", "M5_Unstable3", "M6_Unstable4",
                "M8_CL-like2", "M11_CL-like3")

suffixes <- c("M1", "M7", "M9", "M10", "M16", "M2", "M3", "M4", "M5", "M6", "M8", "M11")

new_patient_ids <- c("p1_normal", 
                 "p2_normal", 
                 "p3_normal", 
                 "p4_normal", 
                 "p5_cancer_untreated", 
                 "p6_cancer_untreated", 
                 "p7_cancer_untreated", 
                 "p8_cancer_untreated", 
                 "p9_cancer_untreated", 
                 "p10_cancer_untreated", 
                 "p11_cancer_untreated", 
                 "p12_cancer_untreated")


patient_condition <- c(
  "normal",   # M1
  "normal",   # M2
  "normal",   # M3
  "normal",   # M4
  "cancer_untreated",   # M5
  "cancer_untreated",   # M6
  "cancer_untreated",   # M7
  "cancer_untreated",   # M8
  "cancer_untreated",   # M9
  "cancer_untreated",   # M10
  "cancer_untreated",   # M11
  "cancer_untreated"
)


# Make sure sample_ids length matches folder count
stopifnot(length(patient_ids) == 12)

# Create folders
sapply(patient_ids, function(x) {
  dir.create(x, showWarnings = TRUE, recursive = TRUE)
})
# Vectors

file_types <- list(
  filtered_matrix = "filtered_feature_bc_matrix.h5",
  annotations     = "pathologist_annotations.csv.gz",
  scalefactors    = "scalefactors_json.json.gz",
  hires_image     = "tissue_hires_image.png.gz",
  positions_list  = "tissue_positions_list.csv.gz"
)

# Loop through each patient and download files into their folder
for (i in seq_along(patient_ids)) {
  pid <- patient_ids[i]
  sid <- suffixes[i]
  
  # Create the directory if it doesn't exist
  if (!dir.exists(pid)) {
    dir.create(pid)
  }
  
  # Loop through each file type
  for (ftype in names(file_types)) {
    fname <- file_types[[ftype]]
    encoded_fname <- URLencode(paste0(pid, "_", sid, "_", fname))
    
    url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", pid, "&format=file&file=", encoded_fname)
    destfile <- file.path(pid, paste0(pid, "_", sid, "_", fname))  # Save inside the folder
    
    # Download the file
    download.file(url = url, destfile = destfile, mode = "wb")
  }
}
# Load required libraries
library(R.utils)
library(Seurat)
library(dplyr)
# Working directory
wd <- getwd()
for (i in seq_along(patient_ids)) {
  pid    <- patient_ids[i]
  suf    <- suffixes[i]
  sampid <- sample_ids[i]
  combo <- paste0(pid, "_", suf)
  patient_dir <- file.path(wd, pid)
  spatial_dir <- file.path(patient_dir, "spatial")
  
  # Create spatial folder if it doesn't exist
  if (!dir.exists(spatial_dir)) dir.create(spatial_dir)
  
  # Step 1: Move and decompress spatial files
  spatial_files <- c("tissue_positions_list.csv.gz",
                     "tissue_hires_image.png.gz",
                     "scalefactors_json.json.gz")
  
  for (file in spatial_files) {
    src <- file.path(patient_dir, paste0(combo, "_", file))
    dst <- file.path(spatial_dir, paste0(combo, "_", file))
    if (file.exists(src)) {
      file.rename(src, dst)
    } else {
      message("Missing file: ", src)
    }
  }
  
  # Step 2: Decompress files
  gz_files <- list.files(spatial_dir, pattern="\\.gz$", full.names=TRUE)
  for (gzf in gz_files) {
    gunzip(gzf, destname=sub("\\.gz$", "", gzf), remove=TRUE)
  }
  
  # Step 3: Rename files to match standard names
  rename_map <- c("tissue_positions_list.csv",
                  "tissue_hires_image.png",
                  "scalefactors_json.json")
  
  for (file in rename_map) {
    oldname <- file.path(spatial_dir, paste0(combo, "_", file))
    newname <- file.path(spatial_dir, file)
    if (file.exists(oldname)) {
      file.rename(oldname, newname)
    } else {
      message("File not found for renaming: ", oldname)
    }
  }
  
  # Step 4: Load matrix and image
  h5f <- file.path(patient_dir, paste0(combo, "_filtered_feature_bc_matrix.h5"))
  if (!file.exists(h5f)) {
    message("Missing H5: ", h5f, " — skipping.")
    next
  }
  matrix <- Read10X_h5(h5f)
  
  imgf <- file.path(spatial_dir, "tissue_hires_image.png")
  if (!file.exists(imgf)) {
    message("Missing image: ", imgf, " — skipping.")
    next
  }
  img <- Read10X_Image(spatial_dir, image.name="tissue_hires_image.png")
  
  # Step 5: Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = matrix, assay = "Spatial")
  seurat_obj@images$spatial <- img
  
  # Step 6: Add annotation if available
  annf <- file.path(patient_dir, paste0(combo, "_pathologist_annotations.csv.gz"))
  if (file.exists(annf)) {
    gunzip(annf, remove=FALSE)
    ann_txt <- sub("\\.gz$", "", annf)
    if (file.exists(ann_txt)) {
      annotations <- read.csv(ann_txt, stringsAsFactors = FALSE)
      colnames(annotations) <- c("barcode", "annotations")
      
      seurat_obj$barcode <- rownames(seurat_obj@meta.data)
      metadata <- left_join(seurat_obj@meta.data, annotations, by = "barcode")
      rownames(metadata) <- metadata$barcode
      metadata$barcode <- NULL
      seurat_obj@meta.data <- metadata
    }
  }
  
  # Step 7: Add metadata columns
  seurat_obj[["sampleid"]]   <- sampid
  seurat_obj[["patient_id"]] <- pid
  seurat_obj[["suffix"]]     <- suf
  seurat_obj[["new_patient_id"]] <- new_patient_ids[i]
  seurat_obj[["condition"]] <- patient_condition[i]
  
  # Step 8: QC & processing
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^MT-")
  seurat_obj <- subset(seurat_obj, subset = nCount_Spatial > 200 & nFeature_Spatial < 7500 & percent.mt < 10)
  seurat_obj <- SCTransform(seurat_obj, assay="Spatial", verbose=FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method="vst", nfeatures=2000)
  seurat_obj <- RunPCA(seurat_obj, verbose=FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction="pca", dims=1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution=0.5)
  seurat_obj <- RunUMAP(seurat_obj, reduction="pca", dims=1:30)
  
  # Step 9: Save
  saveRDS(seurat_obj, file=file.path(patient_dir, paste0(combo, "_seurat.rds")))
}
unique(GSM6592050_M3_seurat @meta.data)
#####################################integration######################
GSM6592048_M1_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592048/GSM6592048_M1_seurat.rds")
GSM6592049_M2_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592049/GSM6592049_M2_seurat.rds")
GSM6592050_M3_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592050/GSM6592050_M3_seurat.rds")
GSM6592051_M4_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592051/GSM6592051_M4_seurat.rds")
GSM6592052_M5_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592052/GSM6592052_M5_seurat.rds")
GSM6592053_M6_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592053/GSM6592053_M6_seurat.rds")
GSM6592054_M7_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592054/GSM6592054_M7_seurat.rds")
GSM6592055_M8_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592055/GSM6592055_M8_seurat.rds")
GSM6592056_M9_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592056/GSM6592056_M9_seurat.rds")
GSM6592057_M10_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592057/GSM6592057_M10_seurat.rds")
GSM6592058_M11_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592058/GSM6592058_M11_seurat.rds")
GSM6592062_M16_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592062/GSM6592062_M16_seurat.rds")
library(Seurat)
library(dplyr)
library(harmony)
library(patchwork)
# 1. Read in all Seurat objects into a list
seurat_list <- list(
  GSM6592048_M1 = GSM6592048_M1_seurat,
  GSM6592049_M2 = GSM6592049_M2_seurat,
  GSM6592050_M3 = GSM6592050_M3_seurat,
  GSM6592051_M4 = GSM6592051_M4_seurat,
  GSM6592052_M5 = GSM6592052_M5_seurat,
  GSM6592053_M6 = GSM6592053_M6_seurat,
  GSM6592054_M7 = GSM6592054_M7_seurat,
  GSM6592055_M8 = GSM6592055_M8_seurat,
  GSM6592056_M9 = GSM6592056_M9_seurat,
  GSM6592057_M10 = GSM6592057_M10_seurat,
  GSM6592058_M11 = GSM6592058_M11_seurat,
  GSM6592062_M16 = GSM6592062_M16_seurat)



# Merge all Seurat objects while ensuring unique cell names
merged <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = sample_ids,  # exclude first ID because it's already in x
  project = "HarmonyIntegration"
)

library(stringr)

library(dplyr)
library(stringr)

merged@meta.data <- merged@meta.data %>%
  mutate(
    annotations = annotations %>%
      str_replace_all("Artefacts|Artifacts", "Artifacts") %>%
      str_replace_all("Tumor Cells|Tumour Cells|Tumor cells|Tumour cells|\\?", "Tumor cells") %>%
      str_replace_all("In [Ss]itu [Cc]arcinoma\\*?", "In situ carcinoma") %>%
      str_replace_all("Tumor Stroma|Tumour Stroma", "Tumor Stroma") %>%
      str_replace_all("Tumor cells Tumor cells", "Tumor cells") %>%
      str_trim()
  ) 





unique(merged@meta.data$annotations)
unique(merged@meta.data$new_patient_id)

# Basic normalization
merged <- NormalizeData(merged, verbose = FALSE)

# Find highly variable genes
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

# Scale data
merged <- ScaleData(merged, verbose = FALSE)

# Run PCA
merged <- RunPCA(merged, verbose = FALSE)

pca_coords <- Embeddings(merged, "pca")


# 5. Run Harmony integration on patients condition 
combined1 <- RunHarmony(merged, group.by.vars = "condition", verbose = TRUE)
combined11 <- RunPCA(combined1, verbose = FALSE)

combined11 <- subset(combined11, 
                     subset = !is.na(annotations) & annotations != "")
unique(combined11@meta.data$annotations)

harmony_coords <- Embeddings(combined11, "harmony")

###################PCA before and after harmoney integration by sample id as a legend ####################
library(Seurat)
library(harmony)
library(irlba)
cols_sampleid <- setNames(rainbow(length(sample_ids)), sample_ids)
# Split sample IDs into 3 groups for 3 lines
sample_ids <- names(cols_sampleid)
n <- length(sample_ids)
split_ids <- split(sample_ids, ceiling(seq_along(sample_ids) / ceiling(n / 3)))

# Setup layout: 2 plots on top row, legend box full width below 
#that is code with the plots with legends below 
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(4, 1))  # 4 units tall plots, 1 unit tall legend

# Plot margins for PCA plots
par(mar = c(5, 5, 4, 2))  # leave space for axis labels and titles

# Plot 1 this one + names from above 
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_sampleid[combined11$sampleid],
  main = "PCA (Before Harmony)\nColored by SampleID",
  pch = 20, xlab = "PC1", ylab = "PC2", bty = "o",xaxt = "s",    # show x-axis
  yaxt = "s"     # show y-axis
)

# Plot 2
plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_sampleid[combined11$sampleid],
  main = "PCA (After Harmony)\nColored by SampleID",
  pch = 20, xlab = "PC1", ylab = "PC2", bty = "o",xaxt = "s",    # show x-axis
  yaxt = "s"     # show y-axis
)

# Now plot legend in bottom row
par(mar = c(0, 0, 0, 0), xpd = TRUE, xaxt = 'n', yaxt = 'n', bty = 'n')
plot.new()

# Legend plotting area dimensions [0,1] in both x and y

# Number of lines
lines_count <- length(split_ids)

# Vertical positions for lines (spread nicely)
y_pos <- seq(0.75, 0.25, length.out = lines_count)

# Horizontal spacing - spread across full width, adjust multiplier for spacing
max_per_line <- max(sapply(split_ids, length))
x_spacing <- 1 / (max_per_line + 1)

# Plot colored points and labels line by line
for (i in seq_along(split_ids)) {
  ids <- split_ids[[i]]
  n_ids <- length(ids)
  
  # Horizontal positions spaced evenly
  x_positions <- seq(from = x_spacing, by = x_spacing, length.out = n_ids)
  y <- y_pos[i]
  
  # Plot colored points
  points(x = x_positions, y = rep(y, n_ids),
         col = cols_sampleid[ids],
         pch = 20, cex = 2)
  
  # Add sample ID labels slightly above points
  text(x = x_positions, y = rep(y + 0.08, n_ids),
       labels = ids, cex = 1.4, adj = 0.5)
}

###################PCA before and after harmoney integration by annotations  as a legend ####################
# Load libraries
library(colorspace)
pca_coords <- Embeddings(merged, "pca")

# Step 3: Harmony integration
harmony_coords <- Embeddings(combined11, "harmony")

# Step 4: Color settings
annotation_groups <- unique(combined11$annotations)
# Define 19 visually distinct colors manually
distinct_colors_20 <- c(
  "#E6194B",  # red
  "#3CB44B",  # green
  "#FFE119",  # yellow
  "#0082C8",  # blue
  "#F58231",  # orange
  "#911EB4",  # purple
  "#46F0F0",  # cyan
  "#F032E6",  # magenta
  "#D2F53C",  # lime
  "#FABEBE",  # pink
  "#008080",  # teal
  "#E6BEFF",  # lavender
  "#AA6E28",  # brown
  "#FFFAC8",  # cream
  "#800000",  # maroon
  "#Aaffc3",  # mint
  "#FFD8B1",  # apricot
  "#000080",  # navy
  "#808080"   # grey
)

# Assign colors for annotations
annotation_groups <- annotation_groups[annotation_groups != ""]
cols_annotations <- setNames(distinct_colors_20[seq_along(annotation_groups)], annotation_groups)

# Split annotation IDs for legend
n <- length(annotation_groups)
split_ids <- split(annotation_groups, ceiling(seq_along(annotation_groups) / ceiling(n / 3)))

# Close open graphic devices
dev.off()

# Layout: PCA Before & After Harmony + Annotation Legend
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(4, 1))

par(mar = c(5, 5, 4, 2))

# Plot 1: PCA Before Harmony
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_annotations[combined11$annotations],
  main = "PCA (Before Harmony)\nColored by Annotations",
  pch = 20, xlab = "PC1", ylab = "PC2"
)

# Plot 2: PCA After Harmony
plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_annotations[combined11$annotations],
  main = "PCA (After Harmony)\nColored by Annotations",
  pch = 20, xlab = "PC1", ylab = "PC2"
)

# Legend area
par(mar = c(0, 0, 0, 0), xpd = TRUE, xaxt = 'n', yaxt = 'n', bty = 'n')
plot.new()

lines_count <- length(split_ids)
y_pos <- seq(0.75, 0.1, length.out = lines_count)
max_per_line <- max(sapply(split_ids, length))

for (i in seq_along(split_ids)) {
  ids <- split_ids[[i]]
  n_ids <- length(ids)
  base_spacing <- 1 / (max_per_line + 1)
  adjusted_spacing <- base_spacing * 3
  start_x <- 0.1
  x_positions <- start_x + (0:(n_ids - 1)) * adjusted_spacing
  
  if (max(x_positions) > 1) {
    x_positions <- seq(from = 0.05, to = 0.95, length.out = n_ids)
  }
  
  y <- y_pos[i]
  points(x = x_positions, y = rep(y, n_ids),
         col = cols_annotations[ids], pch = 20, cex = 2)
  text(x = x_positions, y = rep(y + 0.08, n_ids),
       labels = ids, cex = 1.4, adj = 0.5)
}


###################PCA before and after harmoney integration by patient id as a legend ####################


# Assign colors for new_patient_id
unique_patients <- unique(combined11@meta.data$new_patient_id)
cols_patient <- setNames(rainbow(length(unique_patients)), unique_patients)

### ---- PCA Before & After Harmony Colored by New Patient ID ----

# Layout: 2 plots (top) + 1 legend (bottom)
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(4, 1))

# Margins for PCA plots
par(mar = c(5, 5, 4, 2))

# Plot 1: PCA Before Harmony colored by new_patient_id
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_patient[combined11$new_patient_id],
  main = "PCA (Before Harmony)\nColored by New Patient ID",
  pch = 20, xlab = "PC1", ylab = "PC2",
  bty = "o", xaxt = "s",   # show x-axis
  yaxt = "s"    # show y-axis
)

# Plot 2: PCA After Harmony colored by new_patient_id
plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_patient[combined11$new_patient_id],
  main = "PCA (After Harmony)\nColored by New Patient ID",
  pch = 20, xlab = "PC1", ylab = "PC2",
  bty = "o", xaxt = "s",   # show x-axis
  yaxt = "s"    # show y-axis
)

# Legend for new_patient_id
par(mar = c(0, 0, 0, 0), xpd = TRUE, xaxt = 'n', yaxt = 'n', bty = 'n')
plot.new()

annotation_groups <- unique_patients
n <- length(annotation_groups)
split_ids <- split(annotation_groups, ceiling(seq_along(annotation_groups) / ceiling(n / 3)))

lines_count <- length(split_ids)
y_pos <- seq(0.75, 0.1, length.out = lines_count)
max_per_line <- max(sapply(split_ids, length))

for (i in seq_along(split_ids)) {
  ids <- split_ids[[i]]
  n_ids <- length(ids)
  base_spacing <- 1 / (max_per_line + 1)
  spacing_multiplier <- 3
  adjusted_spacing <- base_spacing * spacing_multiplier
  start_x <- 0.1
  x_positions <- start_x + (0:(n_ids - 1)) * adjusted_spacing
  
  if (max(x_positions) > 1) {
    x_positions <- seq(from = 0.05, to = 0.95, length.out = n_ids)
  }
  
  y <- y_pos[i]
  points(x = x_positions, y = rep(y, n_ids),
         col = cols_patient[ids],
         pch = 20, cex = 2)
  text(x = x_positions, y = rep(y + 0.08, n_ids),
       labels = ids, cex = 1.4, adj = 0.5)
}

dev.off()

############################pca before and after harmoney colored by condition and treatment ####################
condition_groups <- unique(combined11$condition)
cols_condition <- setNames(rainbow(length(condition_groups)), condition_groups)

dev.off()  # Close any open graphics devices

# Layout: two PCA plots top, one legend row bottom
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(4, 1))

# Set margins for PCA plots
par(mar = c(5, 5, 4, 2))

# PCA Plot 1: Before Harmony
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_condition[combined11$condition],
  main = "PCA (Before Harmony)\nColored by Annotations",
  pch = 20, xlab = "PC1", ylab = "PC2",
  bty = "o", xaxt = "s", yaxt = "s"
)

# PCA Plot 2: After Harmony
plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_condition[combined11$condition],
  main = "PCA (After Harmony)\nColored by Annotations",
  pch = 20, xlab = "PC1", ylab = "PC2",
  bty = "o", xaxt = "s", yaxt = "s"
)

# One-line Legend
par(mar = c(0, 0, 0, 0), xpd = TRUE, xaxt = 'n', yaxt = 'n', bty = 'n')
plot.new()

ids <- condition_groups
n_ids <- length(ids)

# Compute horizontal spacing
x_positions <- seq(from = 0.05, to = 0.95, length.out = n_ids)
y <- 0.5  # Vertical center

# Plot colored dots
points(x = x_positions, y = rep(y, n_ids),
       col = cols_condition[ids],
       pch = 20, cex = 2)

# Add labels above the points
text(x = x_positions, y = rep(y + 0.08, n_ids),
     labels = ids, cex = 1.2, adj = 0.5)

########################extract metadata from the serat objects######################
# Loop over each Seurat object in the list
for (sample_name in names(seurat_list)) {
  # Get the Seurat object
  seurat_obj <- seurat_list[[sample_name]]
  
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Add CellID column
  metadata$CellID <- rownames(metadata)
  
  # Reorder columns to have CellID first
  metadata <- metadata[, c("CellID", setdiff(colnames(metadata), "CellID"))]
  
  # Define filename (e.g., GSM6592048_M1_metadata.tsv)
  filename <- paste0(sample_name, "_metadata.tsv")
  
  # Write to TSV file
  write.table(metadata, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

p1_normal_tissue_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592048/p1_normal_tissue_seurat.rds")
p2_normal_tissue_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592054/p2_normal_tissue_seurat.rds")
p3_normal_tissue_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592056/p3_normal_tissue_seurat.rds")
p4_normal_tissue_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592057/p4_normal_tissue_seurat.rds")
p5_cancer_untreated_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592062/p5_cancer_untreated_seurat.rds")
p6_cancer_untreated_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592049/p6_cancer_untreated_seurat.rds")
p7_cancer_untreated_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592050/p7_cancer_untreated_seurat.rds")
p8_cancer_untreated_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592051/p8_cancer_untreated_seurat.rds")
p9_cancer_untreated_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592052/p9_cancer_untreated_seurat.rds")
p10_cancer_untreated_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592053/p10_cancer_untreated_seurat.rds")
p11_cancer_untreated_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592055/p11_cancer_untreated_seurat.rds")
p12_cancer_untreated_seurat <- readRDS("~/integrated by condition normal and untreated only/GSM6592058/p12_cancer_untreated_seurat.rds")
library(tidyverse)
new_seurat_list <- list(
  p1_normal_tissue    = p1_normal_tissue_seurat,
  p2_normal_tissue = p2_normal_tissue_seurat,
  p3_normal_tissue = p3_normal_tissue_seurat,
  p4_normal_tissue = p4_normal_tissue_seurat,
  p5_cancer_untreated = p5_cancer_untreated_seurat,
  p6_cancer_untreated = p6_cancer_untreated_seurat,
  p7_cancer_untreated = p7_cancer_untreated_seurat,
  p8_cancer_untreated = p8_cancer_untreated_seurat,
  p9_cancer_untreated = p9_cancer_untreated_seurat,
  p10_cancer_untreated = p10_cancer_untreated_seurat,
  p11_cancer_untreated = p11_cancer_untreated_seurat,
  p12_cancer_untreated = p12_cancer_untreated_seurat)

# Create a data frame to collect annotation counts
annotation_matrix <- list()

# Loop through each sample
for (sample_name in names(new_seurat_list)) {
  meta <- new_seurat_list[[sample_name]]@meta.data
  annotation_counts <- table(meta$annotations)
  annotation_df <- as.data.frame(annotation_counts)
  colnames(annotation_df) <- c("Annotation", "Count")
  annotation_df$Sample <- sample_name
  annotation_matrix[[sample_name]] <- annotation_df
}

# Combine all into one data frame
annotation_all <- bind_rows(annotation_matrix)
annotation_all <- subset(annotation_all, 
                                   subset = !is.na(Annotation ) & Annotation  != "")
library(dplyr)
library(stringr)

annotation_all <- annotation_all %>%
  mutate(
    Annotation = Annotation %>%
      str_replace_all("Artefacts|Artifacts", "Artifacts") %>%
      str_replace_all("Tumor Cells|Tumour Cells|Tumor cells|Tumour cells|\\?", "Tumor cells") %>%
      str_replace_all("In [Ss]itu [Cc]arcinoma\\*?", "In situ carcinoma") %>%
      str_replace_all("Tumor Stroma|Tumour Stroma", "Tumor Stroma") %>%
      str_replace_all("Tumor cells Tumor cells", "Tumor cells") %>%
      str_trim()
  )

# Convert to wide format for heatmap input
annotation_wide <- annotation_all %>%
  pivot_wider(names_from = Annotation, values_from = Count, values_fill = 0) %>%
  column_to_rownames("Sample")

library(pheatmap)

# Plot
pheatmap(as.matrix(annotation_wide),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Annotation Frequency per Tissue Sample")

################################################################merge all metadata files in one merged file################
library(readxl)
library(dplyr)
library(writexl)
install.packages("writexl")
# Step 1: List all Excel files in your folder
files <- list.files(path = "~/integrated by condition normal and untreated only/metadat excel files", pattern = "\\.xlsx$", full.names = TRUE)

# Step 2: Read and combine all Excel files into one data frame
combined_data <- bind_rows(lapply(files, function(file) {
  df <- read_excel(file)
  df$Source <- tools::file_path_sans_ext(basename(file))  # optional: keep track of original file
  return(df)
}))

# Step 3: Write the combined data into a single-sheet Excel file
write_xlsx(combined_data, path = "All_Metadata_Combined.xlsx")
