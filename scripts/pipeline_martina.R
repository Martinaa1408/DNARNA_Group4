# Load the required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")
library(minfi)

# Clear the R environment
rm(list=ls())

# Set the working directory
setwd("C:/Users/marty/OneDrive/Desktop/github/DNARNA_project/Final_Report-20250524")

# Load the required library again (in case it was unloaded)
library(minfi)

# Define the base directory containing the raw data and load the sample sheet
baseDir <- "C:/Users/marty/OneDrive/Desktop/github/DNARNA_project/Final_Report-20250524/Input_Data"
targets <- read.metharray.sheet(baseDir)

# Create an RGChannelSet object from the raw data
RGset <- read.metharray.exp(targets = targets)
save(RGset, file = "RGset.RData")

# Extract Green and Red channels
Red <- data.frame(getRed(RGset))
dim(Red)  # Shows number of probes (rows) and samples (columns)
head(Red)
Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)

# Extract fluorescence values at a specific address
address <- "44666390"
probe_red <- Red[address, ]
probe_green <- Green[address, ]

# Load the cleaned manifest file to check probe design type
load('Illumina450Manifest_clean.RData')
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID == "44666390", 'Infinium_Design_Type']
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID == "44666390", 'Infinium_Design_Type']

# Create a summary dataframe for the selected address
df_address <- data.frame(
  Sample = colnames(probe_green),
  Red_fluor = unlist(probe_red, use.names = FALSE),
  Green_fluor = unlist(probe_green, use.names = FALSE),
  Type = "I"
)
df_address

# Convert the RGChannelSet to an MSet object using preprocessRaw
BiocManager::install("IlluminaHumanMethylation450kmanifest")
MSet.raw <- preprocessRaw(RGset)
MSet.raw
save(MSet.raw, file = "MSet_raw.RData")

# Quality control
qc <- getQC(MSet.raw)
qc
plotQC(qc)
controlStripPlot(RGset, controls = "NEGATIVE")

# Calculate detection p-values and flag unreliable probes
threshold <- 0.01
detP <- detectionP(RGset)
failed <- detP > threshold
num_failed <- colSums(failed)
df_failed <- data.frame(
  Sample = colnames(probe_green),
  Failed_Position = num_failed,
  perc_failed_probes = colMeans(failed) * 100
)
df_failed

# Select groups (CTRL and DIS) from targets
ctrl <- basename(targets[targets$Group == "CTRL", "Basename"])
dis <- basename(targets[targets$Group == "DIS", "Basename"])

# Remove any suffixes (e.g., .idat)
ctrl <- gsub("\\.idat$", "", ctrl)
dis <- gsub("\\.idat$", "", dis)

# Print the group samples for verification
cat("CTRL samples:", ctrl, "\n")
cat("DIS samples:", dis, "\n")

# Subset the MSet object by groups
ctrlSet <- MSet.raw[, colnames(MSet.raw) %in% ctrl]
disSet <- MSet.raw[, colnames(MSet.raw) %in% dis]

# Calculate Beta and M values
ctrlBeta <- getBeta(ctrlSet)
ctrlM <- getM(ctrlSet)
disBeta <- getBeta(disSet)
disM <- getM(disSet)

# Calculate mean values per probe across samples
mean_ctrlBeta <- apply(ctrlBeta, 1, mean, na.rm = TRUE)
mean_disBeta <- apply(disBeta, 1, mean, na.rm = TRUE)
mean_ctrlM <- apply(ctrlM, 1, mean, na.rm = TRUE)
mean_disM <- apply(disM, 1, mean, na.rm = TRUE)

# Calculate density estimates, checking for enough data points
if (sum(!is.na(mean_ctrlBeta)) >= 2) {
  d_mean_ctrlBeta <- density(na.omit(mean_ctrlBeta))
} else {
  warning("Not enough points in mean_ctrlBeta for density plot.")
}

if (sum(!is.na(mean_disBeta)) >= 2) {
  d_mean_disBeta <- density(na.omit(mean_disBeta))
} else {
  warning("Not enough points in mean_disBeta for density plot.")
}

if (sum(!is.na(mean_ctrlM)) >= 2) {
  d_mean_ctrlM <- density(na.omit(mean_ctrlM))
} else {
  warning("Not enough points in mean_ctrlM for density plot.")
}

if (sum(!is.na(mean_disM)) >= 2) {
  d_mean_disM <- density(na.omit(mean_disM))
} else {
  warning("Not enough points in mean_disM for density plot.")
}

# Plot density of mean Beta and M values
par(mfrow = c(1, 2))  # 1 row, 2 columns

# Density of mean Beta values
plot(d_mean_ctrlBeta,
     main = "Density of Beta Values",
     col = "#0067E6",
     lwd = 2.5,
     xlab = "Mean Beta Values",
     ylab = "Density")
lines(d_mean_disBeta, col = "#E50068", lwd = 2.5)
legend("topright", legend = c("CTRL", "DIS"), fill = c("#0067E6", "#E50068"), cex = 1)

# Density of mean M values
plot(d_mean_ctrlM,
     main = "Density of M Values",
     col = "#0067E6",
     lwd = 2.5,
     xlab = "Mean M Values",
     ylab = "Density")
lines(d_mean_disM, col = "#E50068", lwd = 2.5)
legend("topright", legend = c("CTRL", "DIS"), fill = c("#0067E6", "#E50068"), cex = 1)

# Reset plotting layout
par(mfrow = c(1, 1))

# Subset the manifest by probe chemistry (Type I and Type II)
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "I", ]
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "II", ]
dfI <- droplevels(dfI)
dfII <- droplevels(dfII)

# Subset raw beta values by chemistry
beta_raw <- getBeta(MSet.raw)
beta_I_raw <- beta_raw[rownames(beta_raw) %in% dfI$IlmnID, ]
beta_II_raw <- beta_raw[rownames(beta_raw) %in% dfII$IlmnID, ]

# Calculate mean and sd for raw Beta
mean_beta_I_raw <- apply(beta_I_raw, 1, mean, na.rm = TRUE)
mean_beta_II_raw <- apply(beta_II_raw, 1, mean, na.rm = TRUE)
sd_beta_I_raw <- apply(beta_I_raw, 1, sd, na.rm = TRUE)
sd_beta_II_raw <- apply(beta_II_raw, 1, sd, na.rm = TRUE)

# Density plots for raw Beta values
density_mean_beta_I_raw <- density(na.omit(mean_beta_I_raw))
density_mean_beta_II_raw <- density(na.omit(mean_beta_II_raw))
density_sd_beta_I_raw <- density(na.omit(sd_beta_I_raw))
density_sd_beta_II_raw <- density(na.omit(sd_beta_II_raw))

# Normalize using preprocessFunnorm
MSet.norm <- preprocessFunnorm(RGset)
beta_norm <- getBeta(MSet.norm)

# Subset normalized Beta by chemistry
beta_I_norm <- beta_norm[rownames(beta_norm) %in% dfI$IlmnID, ]
beta_II_norm <- beta_norm[rownames(beta_norm) %in% dfII$IlmnID, ]

# Calculate mean and sd for normalized Beta
mean_beta_I_norm <- apply(beta_I_norm, 1, mean, na.rm = TRUE)
mean_beta_II_norm <- apply(beta_II_norm, 1, mean, na.rm = TRUE)
sd_beta_I_norm <- apply(beta_I_norm, 1, sd, na.rm = TRUE)
sd_beta_II_norm <- apply(beta_II_norm, 1, sd, na.rm = TRUE)

# Density plots for normalized Beta values
density_mean_beta_I_norm <- density(na.omit(mean_beta_I_norm))
density_mean_beta_II_norm <- density(na.omit(mean_beta_II_norm))
density_sd_beta_I_norm <- density(na.omit(sd_beta_I_norm))
density_sd_beta_II_norm <- density(na.omit(sd_beta_II_norm))

# Plotting: 2 rows x 3 columns
par(mfrow = c(2, 3))

# Raw mean Beta densities
plot(density_mean_beta_I_raw, main = "Raw Beta Mean (Type I)", col = "blue", lwd = 2, xlab = "Mean Beta")
lines(density_mean_beta_II_raw, col = "red", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), fill = c("blue", "red"))

# Raw sd Beta densities
plot(density_sd_beta_I_raw, main = "Raw Beta SD (Type I)", col = "blue", lwd = 2, xlab = "Beta SD")
lines(density_sd_beta_II_raw, col = "red", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), fill = c("blue", "red"))

# Raw Beta boxplot
boxplot(beta_raw, outline = FALSE,
        col = rep(c("#0067E6", "#E50068"), length.out = ncol(beta_raw)),
        main = "Raw Beta Boxplot", ylab = "Beta Values")

# Normalized mean Beta densities
plot(density_mean_beta_I_norm, main = "Norm. Beta Mean (Type I)", col = "blue", lwd = 2, xlab = "Mean Beta")
lines(density_mean_beta_II_norm, col = "red", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), fill = c("blue", "red"))

# Normalized sd Beta densities
plot(density_sd_beta_I_norm, main = "Norm. Beta SD (Type I)", col = "blue", lwd = 2, xlab = "Beta SD")
lines(density_sd_beta_II_norm, col = "red", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), fill = c("blue", "red"))

# Normalized Beta boxplot
boxplot(beta_norm, outline = FALSE,
        col = rep(c("#0067E6", "#E50068"), length.out = ncol(beta_norm)),
        main = "Norm. Beta Boxplot", ylab = "Beta Values")

# Reset plotting layout
par(mfrow = c(1, 1))

# Assuming beta_norm is your matrix of normalized beta values
pca_results <- prcomp(t(beta_norm), scale = TRUE)
library(factoextra)
fviz_eig(pca_results, addlabels = TRUE, xlab = 'PC number', 
ylab = '% of variance', barfill = "#0063A6", barcolor = "black")

# PCA plot by Group
targets$Group <- as.factor(targets$Group)
palette(c("#E50068", "#0067E6"))
xlim <- range(pca_results$x[, 1])
ylim <- range(pca_results$x[, 2])

# Increase right margin, reduce top/bottom margins
par(mar = c(5, 4, 4, 6), xpd = TRUE)
plot(pca_results$x[, 1], pca_results$x[, 2],
     cex = 1.2, pch = 19, col = targets$Group,
     xlab = "PC1", ylab = "PC2",
     xlim = xlim, ylim = ylim,
     main = "PCA (Groups)")
text(pca_results$x[, 1], pca_results$x[, 2],
     labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("right", inset = c(-0.15, 0),
       legend = levels(targets$Group),
       col = c(1, 2), pch = 19,
       cex = 0.9, pt.cex = 1.2, x.intersp = 1.2,
       box.lty = 0, bg = "white")

# PCA plot by Sex
targets$Sex <- as.factor(targets$Sex)
palette(c("#C364CA", "#8BC4F9"))
par(mar = c(5, 4, 4, 6), xpd = TRUE)
plot(pca_results$x[, 1], pca_results$x[, 2],
     cex = 1.2, pch = 19, col = targets$Sex,
     xlab = "PC1", ylab = "PC2",
     xlim = xlim, ylim = ylim,
     main = "PCA (Sex)")
text(pca_results$x[, 1], pca_results$x[, 2],
     labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("right", inset = c(-0.15, 0),
       legend = levels(targets$Sex),
       col = 1:nlevels(targets$Sex), pch = 19,
       cex = 0.9, pt.cex = 1.2, x.intersp = 1.2,
       box.lty = 0, bg = "white")

# PCA plot by Slide (Batch)
targets$Slide <- as.factor(targets$Slide)
palette(c("#9e0059", "green"))
par(mar = c(5, 4, 4, 6), xpd = TRUE)
plot(pca_results$x[, 1], pca_results$x[, 2],
     cex = 1.2, pch = 19, col = targets$Slide,
     xlab = "PC1", ylab = "PC2",
     xlim = xlim, ylim = ylim,
     main = "PCA (Batch)")
text(pca_results$x[, 1], pca_results$x[, 2],
     labels = rownames(pca_results$x), cex = 0.4, pos = 3)
legend("right", inset = c(-0.15, 0),
       legend = levels(targets$Slide),
       col = c(1, 2), pch = 19,
       cex = 0.9, pt.cex = 1.2, x.intersp = 1.2,
       box.lty = 0, bg = "white")

# Define the t-test function: compares beta values between groups
My_t_test <- function(x) {
  t_test <- t.test(x ~ targets$Group)
  return(t_test$p.value)
}

# Apply the t-test function to each probe (row)
p_values <- apply(beta_norm, 1, My_t_test)

# Create a data frame combining beta values with the raw p-values
final_ttest <- data.frame(beta_norm, p_values_ttest = p_values)

# Sort the probes by p-value in ascending order
final_ttest <- final_ttest[order(final_ttest$p_values_ttest), ]

# Filter significant probes with p-value <= 0.01
final_ttest_th <- final_ttest[final_ttest$p_values_ttest <= 0.01, ]

# Plot histogram of p-value distribution
hist(final_ttest$p_values_ttest, main = "P-value distribution (t-test)", 
     xlab = 'p-value')
abline(v = 0.01, col = "#ef233c")

# Check the number of differentially methylated probes before correction
dim(final_ttest_th)

# Apply multiple testing corrections
corr_pval_BH <- p.adjust(final_ttest$p_values_ttest, method = "BH")  # Benjamini-Hochberg
corr_pval_bonferroni <- p.adjust(final_ttest$p_values_ttest, method = "bonferroni")

# Combine the corrected p-values into the final data frame
final_ttest_corr <- data.frame(final_ttest, corr_pval_BH, corr_pval_bonferroni)

# Count the number of differentially methylated probes before and after correction
before_correction <- nrow(final_ttest_corr[final_ttest_corr$p_values_ttest <= 0.01, ])
after_Bonferroni <- nrow(final_ttest_corr[final_ttest_corr$corr_pval_bonferroni <= 0.01, ])
after_BH <- nrow(final_ttest_corr[final_ttest_corr$corr_pval_BH <= 0.01, ])

# Create a summary data frame with the results
diff_meth_df <- data.frame(before_correction, after_Bonferroni, after_BH)
rownames(diff_meth_df) <- c("# differentially methylated probes")

# Print the summary data frame
print(diff_meth_df)

# Set plotting area to a single plot
par(mfrow = c(1,1))

# Plot boxplots of raw and adjusted p-values
boxplot(final_ttest_corr[,9:11], 
        col = c("#7038FF", "#C7FF38", "seashell2"),
        names = c("Raw", "BH", "Bonferroni"),
        main = 'P-value distribution (raw and adjusted)',
        ylab = "P-value")

# Add legend to the plot
legend("topright", 
       legend = c("raw", "BH", "Bonferroni"),
       col = c("#7038FF", "#C7FF38", "seashell2"),
       pch = 19, cex = 0.5, xpd = TRUE)

# Extract beta values (columns 1 to 8)
beta <- final_ttest_corr[, 1:8]

# Calculate mean beta-values for each group
beta_groupCTRL <- beta[, targets$Group == "CTRL"]
mean_beta_groupCTRL <- apply(beta_groupCTRL, 1, mean)

beta_groupDIS <- beta[, targets$Group == "DIS"]
mean_beta_groupDIS <- apply(beta_groupDIS, 1, mean)

# Calculate delta (DIS - CTRL)
delta_average <- mean_beta_groupDIS - mean_beta_groupCTRL

# Identify probes significant after BH correction
BH_sig <- final_ttest_corr$corr_pval_BH < 0.01

# Create a data frame for the volcano plot
toVolcPlot <- data.frame(
  delta_average,
  minus_log10_p_val = -log10(final_ttest_corr$p_values_ttest),
  BH_sig
)

# Plot the volcano plot
plot(toVolcPlot[,1], toVolcPlot[,2], 
     pch=16, cex=0.4, col='grey16',
     xlab='Δμ_beta(DIS-CTRL)', ylab='-log10(p)',
     xlim=c(-0.8,0.8))
abline(h = -log10(0.01))

# Add points for nominally significant probes
nominal_sig <- toVolcPlot[
  abs(toVolcPlot[,1]) > 0.1 & 
    toVolcPlot[,2] > -log10(0.01), ]
points(nominal_sig[,1], nominal_sig[,2], pch=16, cex=0.4, col="snow4")

# Add points for BH-significant probes
BH_sig_points <- toVolcPlot[
  abs(toVolcPlot[,1]) > 0.1 & 
    toVolcPlot$BH_sig, ]
points(BH_sig_points[,1], BH_sig_points[,2], pch=19, cex=0.6, col="hotpink3")

# Add legend
legend("topright", 
       legend=c("not significant", "p < 0.01", "BH p < 0.01"),
       col=c("grey16", "snow4", "hotpink3"), pch = 19)

# Manhattan Plot

# Step 1: Add probe IDs as a column
final_ttest_corr_anno <- data.frame(rownames(final_ttest_corr), final_ttest_corr)
colnames(final_ttest_corr_anno)[1] <- "IlmnID"

# Step 2: Merge with Illumina manifest file (must contain IlmnID, CHR, MAPINFO)
final_ttest_corr_anno <- merge(final_ttest_corr_anno, Illumina450Manifest_clean, by = "IlmnID")

# Step 3: Prepare the data frame for qqman
input_Manhattan <- data.frame(
  ID = final_ttest_corr_anno$IlmnID,
  CHR = final_ttest_corr_anno$CHR,
  MAPINFO = final_ttest_corr_anno$MAPINFO,
  PVAL = final_ttest_corr_anno$p_values_ttest
)

# Step 4: Convert chromosomes X and Y to numeric (23 and 24)
levels(input_Manhattan$CHR)[levels(input_Manhattan$CHR) == "X"] <- "23"
levels(input_Manhattan$CHR)[levels(input_Manhattan$CHR) == "Y"] <- "24"
input_Manhattan$CHR <- as.numeric(as.character(input_Manhattan$CHR))

# Step 5: Load the qqman library and define a color palette
library(qqman)
color_palette <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2",
  "#7F7F7F", "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896",
  "#C5B0D5", "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"
)

# Step 6: Plot the Manhattan plot
manhattan(
  input_Manhattan,
  snp = "ID",
  chr = "CHR",
  bp = "MAPINFO",
  p = "PVAL",
  annotatePval = 0.00001,
  col = color_palette,
  suggestiveline = FALSE,
  genomewideline = -log10(0.00001)
)

# Load required library
if (!require("gplots")) install.packages("gplots")
library(gplots)

# Prepare the input matrix for the heatmap
input_heatmap <- as.matrix(final_ttest[1:100, 1:8])

# Assign colors to groups
group_color <- c()
for (name in colnames(beta)) {
  if (targets$Group[which(colnames(beta) == name)] == "CTRL") {
    group_color <- c(group_color, "#CBCBCB")  # Light gray
  } else {
    group_color <- c(group_color, "#393939")  # Dark gray
  }
}

# Save Average linkage heatmap
png("Average_linkage_heatmap.png", width = 1200, height = 1000)
heatmap.2(input_heatmap,
          col = terrain.colors(100),
          Rowv = TRUE,
          Colv = TRUE,
          hclustfun = function(x) hclust(x, method = 'average'),
          dendrogram = "both",
          key = TRUE,
          ColSideColors = group_color,
          density.info = "none",
          trace = "none",
          scale = "none",
          symm = FALSE,
          main = "Average linkage",
          key.xlab = 'beta-val',
          key.title = NA,
          keysize = 1,
          labRow = NA)
legend("topright",
       legend = levels(targets$Group),
       col = c('#CBCBCB', '#393939'),
       pch = 19,
       cex = 0.7)
dev.off()

# Save Complete linkage heatmap
png("Complete_linkage_heatmap.png", width = 1200, height = 1000)
heatmap.2(input_heatmap,
          col = terrain.colors(100),
          Rowv = TRUE,
          Colv = TRUE,
          dendrogram = "both",
          key = TRUE,
          ColSideColors = group_color,
          density.info = "none",
          trace = "none",
          scale = "none",
          symm = FALSE,
          main = "Complete linkage",
          key.xlab = 'beta-val',
          key.title = NA,
          keysize = 1,
          labRow = NA)
legend("topright",
       legend = levels(targets$Group),
       col = c('#CBCBCB', '#393939'),
       pch = 19,
       cex = 0.7)
dev.off()

# Save Single linkage heatmap
png("Single_linkage_heatmap.png", width = 1200, height = 1000)
heatmap.2(input_heatmap,
          col = terrain.colors(100),
          Rowv = TRUE,
          Colv = TRUE,
          hclustfun = function(x) hclust(x, method = 'single'),
          dendrogram = "both",
          key = TRUE,
          ColSideColors = group_color,
          density.info = "none",
          trace = "none",
          scale = "none",
          symm = FALSE,
          main = "Single linkage",
          key.xlab = 'beta-val',
          key.title = NA,
          keysize = 1,
          labRow = NA)
legend("topright",
       legend = levels(targets$Group),
       col = c('#CBCBCB', '#393939'),
       pch = 19,
       cex = 0.7)
dev.off()
