library(Seurat)
library(dplyr)

# import the integrated object
seu <- readRDS("data/integrated_prdm1_180924.rds")

# how the atac fragments look like (we have suffix "-1")
frag_head <- read.table(
  gzfile("data/atac_fragments.tsv.gz"),
  sep = "\t",
  nrows = 10,
  stringsAsFactors = FALSE
)
frag_head

# check barcode format in the seurat object (contains also prefix)
head(colnames(seu))
tail(colnames(seu))
unique(seu$predicted.id)

# we only want to keep the ones with prefix "ctl_". Subset CTL cells
ctl_cells <- colnames(seu)[grepl("^ctl_", colnames(seu))]
length(ctl_cells)
ctl <- seu[, ctl_cells]

# we extract barcodes correctly
clean_barcodes <- sub("^ctl_", "", colnames(ctl))
# check 
head(clean_barcodes)

# build the correct cell_data.tsv
cell_data <- data.frame(
  barcode   = clean_barcodes,
  cell_type = ctl$predicted.id,
  sample_id = "ctl",
  stringsAsFactors = FALSE
)

colnames(cell_data) 

# final checks
frag_barcodes <- read.table(
  gzfile("data/atac_fragments.tsv.gz"),
  sep = "\t",
  stringsAsFactors = FALSE
)[,4] |> unique()

cat("Unique fragment barcodes:", length(unique(frag_barcodes)), "\n")
cat("Cell data barcodes:", nrow(cell_data), "\n")
cat("Overlap:", sum(cell_data$barcode %in% frag_barcodes), "\n")

# save dataframe as tsv
write.table(
  cell_data,
  file = "data/cell_data.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

