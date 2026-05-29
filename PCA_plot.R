suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# -----------------------------
# Input files
# -----------------------------
geno_file <- "pca_input.raw"
race_file <- "racefile.txt"

# -----------------------------
# Read genotype data
# -----------------------------
geno <- read.table(geno_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Keep IDs
ids <- geno[, c("FID", "IID")]

# Remove non-genotype columns
# PLINK .raw format first 6 columns are:
# FID, IID, PAT, MAT, SEX, PHENOTYPE
X <- geno[, -(1:6), drop = FALSE]

# Convert to numeric matrix
X <- as.matrix(X)
mode(X) <- "numeric"

# -----------------------------
# Handle missing values
# Replace NA with SNP mean
# -----------------------------
for (j in seq_len(ncol(X))) {
  if (anyNA(X[, j])) {
    X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
  }
}

# -----------------------------
# Remove zero-variance SNPs
# -----------------------------
snp_sd <- apply(X, 2, sd, na.rm = TRUE)
keep <- !is.na(snp_sd) & snp_sd > 0
X <- X[, keep, drop = FALSE]

# Store kept SNP names
kept_snps <- colnames(X)

# -----------------------------
# Standardize SNPs
# -----------------------------
X_scaled <- scale(X, center = TRUE, scale = TRUE)

# -----------------------------
# Run PCA
# -----------------------------
pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)

# -----------------------------
# Save sample scores
# -----------------------------
scores <- as.data.frame(pca$x, stringsAsFactors = FALSE)
scores <- cbind(ids, scores)
write.table(scores, "pca_scores.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Save SNP loadings
# -----------------------------
loadings <- as.data.frame(pca$rotation, stringsAsFactors = FALSE)
loadings$SNP <- kept_snps
loadings <- loadings[, c("SNP", setdiff(colnames(loadings), "SNP"))]
write.table(loadings, "pca_loadings.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Save variance explained
# -----------------------------
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
var_table <- data.frame(
  PC = paste0("PC", seq_along(var_exp)),
  Variance = var_exp,
  Cumulative = cumsum(var_exp),
  stringsAsFactors = FALSE
)
write.table(var_table, "pca_variance.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Read race file and merge
# racefile.txt expected columns: FID IID race
# -----------------------------
race <- read.table(race_file, header = TRUE, stringsAsFactors = FALSE)
scores2 <- merge(scores, race, by = c("FID", "IID"))

# -----------------------------
# Plot settings to mimic your MDS plot style
# -----------------------------
shape_vals <- c(
  "EUR" = 1,
  "ASN" = 1,
  "AMR" = 1,
  "AFR" = 1,
  "OWN" = 3
)

color_vals <- c(
  "EUR" = "green",
  "ASN" = "red",
  "AMR" = "orange",
  "AFR" = "blue",
  "OWN" = "black"
)

# Make sure race is ordered similarly
scores2$race <- factor(scores2$race, levels = c("EUR", "ASN", "AMR", "AFR", "OWN"))

# Only label OWN samples that are not HG or NA, matching your MDS logic
label_data <- subset(scores2, race == "OWN" & !grepl("HG", IID) & !grepl("NA", IID))

# Axis labels with variance explained
pc1_lab <- paste0("PC1 (", round(100 * var_exp[1], 2), "%)")
pc2_lab <- paste0("PC2 (", round(100 * var_exp[2], 2), "%)")
pc3_lab <- paste0("PC3 (", round(100 * var_exp[3], 2), "%)")

# Use the same label rule as MDS_plot3.png
label_data <- subset(scores2, race == "OWN" & !grepl("HG", IID) & !grepl("NA", IID))

# -----------------------------
# PCA plot 1: PC1 vs PC2
# -----------------------------
p1 <- scores2 %>%
  arrange(race) %>%
  ggplot(aes(x = PC1, y = PC2, color = race, shape = race)) +
  geom_point(size = 2) +
  scale_shape_manual(values = c(1, 1, 1, 1, 3), drop = FALSE) +
  scale_color_manual(values = color_vals, drop = FALSE) +
  geom_text_repel(
    data = label_data,
    aes(x = PC1, y = PC2, label = IID),
    color = "black",
    segment.color = "gray",
    segment.size = 0.5,
    nudge_x = 0.01,
    nudge_y = 0.05
  ) +
  theme(text = element_text(size = 20)) +
  labs(x = pc1_lab, y = pc2_lab, color = "race", shape = "race")

ggsave("PCA_plot.png", plot = p1, height = 8, width = 10, units = "in", dpi = 300)

# -----------------------------
# PCA plot 2: PC3 vs PC2
# -----------------------------
p2 <- scores2 %>%
  arrange(race) %>%
  ggplot(aes(x = PC3, y = PC2, color = race, shape = race)) +
  geom_point(size = 2) +
  scale_shape_manual(values = c(1, 1, 1, 1, 3), drop = FALSE) +
  scale_color_manual(values = color_vals, drop = FALSE) +
  geom_text_repel(
    data = label_data,
    aes(x = PC3, y = PC2, label = IID),
    color = "black",
    segment.color = "gray",
    segment.size = 0.5,
    nudge_x = 0.01,
    nudge_y = 0.05
  ) +
  theme(text = element_text(size = 20)) +
  labs(x = pc3_lab, y = pc2_lab, color = "race", shape = "race")

ggsave("PCA_plot2.png", plot = p2, height = 8, width = 10, units = "in", dpi = 300)