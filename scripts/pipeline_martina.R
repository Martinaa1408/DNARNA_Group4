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

# Group: Convert to factor
targets$Group <- as.factor(targets$Group)  

# Set color palette for Group plot
palette(c("#E50068", "#0067E6"))  

# Plot PCA with points colored by Group
plot(pca_results$x[, 1], pca_results$x[, 2], 
     cex = 1, pch = 19, col = targets$Group,
     xlab = "PC1", ylab = "PC2", 
     xlim = c(-700, 700), ylim = c(-700, 700),
     main = 'PCA (Groups)')  
# Add sample labels
text(pca_results$x[, 1], pca_results$x[, 2], 
     labels = rownames(pca_results$x), cex = 0.4, pos = 3)  
# Add legend
legend("bottomright", legend = levels(targets$Group), 
       col = c(1, 2), pch = 19, cex = 1.0)  

# Sex: Convert to factor
targets$Sex <- as.factor(targets$Sex)  

# Set color palette for Sex plot
palette(c("#C364CA", "#8BC4F9"))  

# Plot PCA with points colored by Sex
plot(pca_results$x[, 1], pca_results$x[, 2], 
     cex = 1, pch = 19, col = targets$Sex,
     xlab = "PC1", ylab = "PC2", 
     xlim = c(-700, 700), ylim = c(-700, 700),
     main = 'PCA (Sex)')  
# Add sample labels
text(pca_results$x[, 1], pca_results$x[, 2], 
     labels = rownames(pca_results$x), cex = 0.4, pos = 3)  
# Add legend
legend("bottomright", legend = levels(targets$Sex), 
       col = c(1:nlevels(targets$Sex)), pch = 19, cex = 1.0)  

# Slide: Convert to factor (for batch plotting)
targets$Slide <- as.factor(targets$Slide)  

# Set color palette for Slide plot
palette(c("#9e0059", "green"))  

# Plot PCA with points colored by Slide
plot(pca_results$x[, 1], pca_results$x[, 2], 
     cex = 1, pch = 19, col = targets$Slide,
     xlab = "PC1", ylab = "PC2", 
     xlim = c(-700, 700), ylim = c(-700, 700),
     main = 'PCA (Batch)')  
# Add sample labels
text(pca_results$x[, 1], pca_results$x[, 2], 
     labels = rownames(pca_results$x), cex = 0.4, pos = 3)  
# Add legend
legend("bottomright", legend = levels(targets$Slide), 
       col = c(1, 2), pch = 19, cex = 1.0)  
