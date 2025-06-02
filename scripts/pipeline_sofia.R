if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")
library(minfi)

# Set the directory in which the raw data are stored
rm(list=ls())
setwd("C:/Users/sofia/OneDrive/Documenti/DRD - R/Final_Report-20250526/Input_Data")
library(minfi)

# Load the samplesheet using the function read,metharray.sheet
baseDir <- ("C:/Users/sofia/OneDrive/Documenti/DRD - R/Final_Report-20250526/Input_Data")
targets <- read.metharray.sheet(baseDir)

# Create an object of class RGChannelSet using the function read.metharray.exp
RGset <- read.metharray.exp(targets = targets)
save(RGset,file="RGset.RData")

# We extract the Green and Red Channels using the function getGreen and getRed
Red <- data.frame(getRed(RGset))
dim(Red) #rows: probes; columns: samples
# [1] 622399      8
head(Red)
Green <- data.frame(getGreen(RGset))
dim(Green)
# [1] 622399      8
head(Green)

# Address of interest
address <- "44666390"

#Get Red and Green fluorescence values for the specified address
probe_red <- Red[address, ]
probe_green <- Green[address, ]

# Load the cleaned manifest file to check the type of the probe
load("C:/Users/sofia/OneDrive/Documenti/DRD - R/DRD/Illumina450Manifest_clean.RData")

# check whether it is of type I or II from column 'Infinium_Design_Type'
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="44666390", 'Infinium_Design_Type']
# factor()
# Levels: I II
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="44666390", 'Infinium_Design_Type']
# [1] I
# Levels: I II
# type I

#Extract all data related to sample, red fluorescence, green fluorescence and type and color to fill the datable
df_address <- data.frame(Sample=colnames(probe_green), Red_fluor=unlist(probe_red, use.names = FALSE ), Green_fluor=unlist(probe_green, use.names = FALSE ), Type = "I")
df_address

# Convert the RGChannelSet object to an MSet object:
BiocManager::install("IlluminaHumanMethylation450kmanifest")
MSet.raw <- preprocessRaw(RGset) ## Use the preprocessRaw function to create the MSet.raw object from the RGChannelSet.
MSet.raw
save(MSet.raw,file="MSet_raw.RData")

# Quality control
qc <- getQC(MSet.raw)
qc
plotQC(qc)
controlStripPlot(RGset, controls="NEGATIVE")

# Set the detection p-value threshold to determine which probes are considered reliably detected.
threshold <- 0.01
detP <- detectionP(RGset)
failed <- detP>threshold
num_failed = colSums(failed)
df_failed <- data.frame(Sample=colnames(probe_green), Failed_Position=num_failed, 
                        perc_failed_probes=(means_of_columns <- colMeans(failed))*100, 
                        row.names = NULL)
df_failed

#Select groups (CTRL and DIS) from targets:
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

