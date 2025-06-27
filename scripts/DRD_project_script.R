# 0. Load required libraries, clear the environment and set the working directory
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")
library(minfi)

rm(list = ls())
directory = "C:/Users/aurim/Documents/Università/4anno/DNA_RNA_dynamics/Final_Report-20250529" #insert here your directory
setwd(directory)

#---------------------
# 1. Load raw data with minfi and create an object called RGset storing the RGChannelSet object

# Define the base directory containing the raw data and load the sample sheet
baseDir <- (directory)
targets <- read.metharray.sheet(baseDir)
# verify and print the number of imported samples
print(paste("Number of found samples: ", nrow(targets))) 
#note: for every sample we have: SampleID, Sex, Group, Sentrix_ID, Sentrix_Position

# Create an RGChannelSet object from the raw data
RGset <- read.metharray.exp(targets = targets)
save(RGset, file = "RGset.RData")



#---------------------
# 2.Create the dataframes Red and Green to store the red and green fluorescences respectively

# Extract Green and Red channels
Red <- data.frame(getRed(RGset))
dim(Red)  # Shows number of probes (rows) and samples (columns)
head(Red)
Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)


#---------------------
# 3. What are the Red and Green fluorescences for the address assigned to your group? 

# Extract fluorescence values at a specific address
address <- "44666390"
probe_red <- Red[address, ]
probe_green <- Green[address, ]
print(probe_red)
print(probe_green)

# 3. OPTIONAL: check in the manifest file if the address corresponds to a Type I or a 
#     Type II probe and, in case of Type I probe, report its color.

# Load the cleaned manifest file to check probe design type
load(file.path(directory, "Illumina450Manifest_clean.RData"))

# Check if the address is in AddressA_ID or AddressB_ID
type_A <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID == "44666390", 'Infinium_Design_Type']
type_B <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID == "44666390", 'Infinium_Design_Type']
paste("Type from AddressA_ID:",type_A)
paste("Type from AddressB_ID:",type_B)

# Determine the actual type 
actual_type <- ifelse(length(type_A) > 0 && !is.na(type_A), type_A, 
                      ifelse(length(type_B) > 0 && !is.na(type_B), type_B, "Not found"))
paste("Type of address:",actual_type)

#Determine the probe color 
probe_row <- Illumina450Manifest_clean[
  Illumina450Manifest_clean$AddressA_ID == "44666390" | 
    Illumina450Manifest_clean$AddressB_ID == "44666390", ]

actual_color_text <- as.character(probe_row$Color_Channel[1])
if (actual_color_text == "Grn") {
  actual_color_text <- "Green"
}
paste("Color:",actual_color_text)

#Create a summary dataframe for the selected address, 
#combining the Red and Green fluorescence values for each sample 
#and including the probe type and color.
df_address <- data.frame(
  Sample = colnames(probe_green),
  Red_fluor = unlist(probe_red, use.names = FALSE),
  Green_fluor = unlist(probe_green, use.names = FALSE),
  Type = actual_type,
  Color = actual_color_text
)
df_address

#---------------------
# 4. Create the object MSet.raw

# Convert the RGChannelSet to an MSet object using preprocessRaw
BiocManager::install("IlluminaHumanMethylation450kmanifest")
MSet.raw <- preprocessRaw(RGset)
MSet.raw
save(MSet.raw,file="MSet_raw.RData")


#---------------------
# 5. Perform the following quality checks and provide a brief comment to each step:
#   5a. QCplot

# Quality control
qc <- getQC(MSet.raw)
qc
png("qc_plot.png", width = 1200, height = 1000)
plotQC(qc)
dev.off()


#   5b. Check the intensity of negative controls using minfi

# Plot for the controls
png("negative_controls.png", width = 1200, height = 1000)
controlStripPlot(RGset, controls = "NEGATIVE")
dev.off()

#   5c. Calculate detection pValues; for each sample, how many probes have a detection p-value higher than the threshold assigned to each group?

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


#---------------------
# 6. Calculate raw beta and M values and plot the densities of mean methylation values, 
#     dividing the samples in CTRL and DIS (suggestion: subset the beta and M values matrixes 
#     in order to retain CTRL or DIS subjects and apply the function mean to the 2 subsets). 
#     Do you see any difference between the two groups?

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
ctrlSet <- MSet.raw[, ctrl]
disSet <- MSet.raw[, dis]

# Calculate Beta and M values
ctrlBeta <- getBeta(ctrlSet)
ctrlM <- getM(ctrlSet)
disBeta <- getBeta(disSet)
disM <- getM(disSet)

# Calculate mean values per probe across samples
mean_ctrlBeta <- rowMeans(ctrlBeta, na.rm = TRUE)
mean_disBeta <- rowMeans(disBeta, na.rm = TRUE)
mean_ctrlM <- rowMeans(ctrlM, na.rm = TRUE)
mean_disM <- rowMeans(disM, na.rm = TRUE)

# Plot density of mean Beta and M values
png("beta_m_values.png", width = 1200, height = 1000)
par(mfrow = c(1, 2))

# Density of mean Beta values
plot(density(mean_ctrlBeta, na.rm = TRUE),
     main = "Density of Beta Values",
     col = "#0067E6", lwd = 2.5,
     xlab = "Mean Beta Values", ylab = "Density")
lines(density(mean_disBeta, na.rm = TRUE), col = "#E50068", lwd = 2.5)
legend("topright", legend = c("CTRL", "DIS"), col = c("#0067E6", "#E50068"), lwd = 2)

# Density of mean M values
plot(density(mean_ctrlM, na.rm = TRUE),
     main = "Density of M Values",
     col = "#0067E6", lwd = 2.5,
     xlab = "Mean M Values", ylab = "Density")
lines(density(mean_disM, na.rm = TRUE), col = "#E50068", lwd = 2.5)
legend("topright", legend = c("CTRL", "DIS"), col = c("#0067E6", "#E50068"), lwd = 2)
dev.off()

par(mfrow = c(1, 1))


#---------------------
# 7. Normalize the data using the function assigned to each group and compare raw data 
#     and normalized data. Produce a plot with 6 panels in which, for both raw and 
#     normalized data, you show the density plots of beta mean values according to the 
#     chemistry of the probes, the density plot of beta standard deviation values 
#     according to the chemistry of the probes and the boxplot of beta values. 
#     Provide a short comment about the changes you observe. 

#     OPTIONAL: do you think that the normalization approach that you used is appropriate 
#     considering this specific dataset? Try to color the boxplots according to the group 
#     (CTRL and DIS) and check whether the distribution of methylation values is different 
#     between the two groups, before and after normalization ...

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
mean_beta_I_raw <- rowMeans(beta_I_raw, na.rm = TRUE)
mean_beta_II_raw <- rowMeans(beta_II_raw, na.rm = TRUE)
sd_beta_I_raw <- apply(beta_I_raw, 1, sd, na.rm = TRUE)
sd_beta_II_raw <- apply(beta_II_raw, 1, sd, na.rm = TRUE)

# Normalize using preprocessFunnorm
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
MSet.norm <- preprocessFunnorm(RGset)
beta_norm <- getBeta(MSet.norm)

# Subset normalized Beta by chemistry
beta_I_norm <- beta_norm[rownames(beta_norm) %in% dfI$IlmnID, ]
beta_II_norm <- beta_norm[rownames(beta_norm) %in% dfII$IlmnID, ]

# Calculate mean and sd for normalized Beta
mean_beta_I_norm <- rowMeans(beta_I_norm, na.rm = TRUE)
mean_beta_II_norm <- rowMeans(beta_II_norm, na.rm = TRUE)
sd_beta_I_norm <- apply(beta_I_norm, 1, sd, na.rm = TRUE)
sd_beta_II_norm <- apply(beta_II_norm, 1, sd, na.rm = TRUE)

#Density plots for normalized Beta values.

density_mean_beta_I_norm <- density(na.omit(mean_beta_I_norm))
density_mean_beta_II_norm <- density(na.omit(mean_beta_II_norm))
density_sd_beta_I_norm <- density(na.omit(sd_beta_I_norm))
density_sd_beta_II_norm <- density(na.omit(sd_beta_II_norm))


# Prepare colors for groups
sample_prefixes <- sapply(strsplit(colnames(beta_raw), "_"), `[`, 1)
group_colors <- ifelse(sample_prefixes %in% targets[targets$Group == "CTRL", "SampleID"], 
                       "#0067E6", "#E50068")

# Plotting: 2 rows x 3 columns
png("Raw_normalised_beta.png", width = 1200, height = 1000)
par(mfrow = c(2, 3))

# RAW DATA (first row)
# 1. Raw mean Beta densities
plot(density(mean_beta_I_raw, na.rm = TRUE), main = "Raw Beta Mean", col = "blue", lwd = 2, xlab = "Mean Beta")
lines(density(mean_beta_II_raw, na.rm = TRUE), col = "red", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), col = c("blue", "red"), lwd = 2)

# 2. Raw sd Beta densities
plot(density(sd_beta_I_raw, na.rm = TRUE), main = "Raw Beta SD", col = "blue", lwd = 2, xlab = "Beta SD")
lines(density(sd_beta_II_raw, na.rm = TRUE), col = "red", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), col = c("blue", "red"), lwd = 2)

# 3. Raw Beta boxplot (colored by group)
boxplot(beta_raw, outline = FALSE, col = group_colors,
        main = "Raw Beta Boxplot", ylab = "Beta Values", names = rep("", ncol(beta_raw)))
legend("topright", legend = c("CTRL", "DIS"), fill = c("#0067E6", "#E50068"))

# NORMALIZED DATA (second row)
# 4. Normalized mean Beta densities
plot(density(mean_beta_I_norm, na.rm = TRUE), main = "Norm. Beta Mean", col = "blue", lwd = 2, xlab = "Mean Beta")
lines(density(mean_beta_II_norm, na.rm = TRUE), col = "red", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), col = c("blue", "red"), lwd = 2)

# 5. Normalized sd Beta densities
plot(density(sd_beta_I_norm, na.rm = TRUE), main = "Norm. Beta SD", col = "blue", lwd = 2, xlab = "Beta SD")
lines(density(sd_beta_II_norm, na.rm = TRUE), col = "red", lwd = 2)
legend("topright", legend = c("Type I", "Type II"), col = c("blue", "red"), lwd = 2)

# 6. Normalized Beta boxplot (colored by group)
boxplot(beta_norm, outline = FALSE, col = group_colors,
        main = "Norm. Beta Boxplot", ylab = "Beta Values", names = rep("", ncol(beta_norm)))
legend("topright", legend = c("CTRL", "DIS"), fill = c("#0067E6", "#E50068"))
dev.off()

par(mfrow = c(1, 1))



#---------------------
# 8. Perform a PCA on the matrix of normalized beta values generated in step 7, 
#     after normalization. Comment the plot (Do you see any outlier? Do the samples 
#     divide according to the group? Do they divide according to the sex of the samples? 
#     Do they divide according to the batch, that is the column Sentrix_ID?).

# Assuming beta_norm is your matrix of normalized beta values
pca_results <- prcomp(t(beta_norm), scale = TRUE)
install.packages(c("ggplot2", "ggpubr", "cluster", "factoMineR"))
install.packages("factoextra")
library(factoextra)

#Scree plot
png("scree_plot.png", width = 1200, height = 1000)
fviz_eig(pca_results, addlabels = TRUE, xlab = 'PC number', 
         ylab = '% of variance', barfill = "#0063A6", barcolor = "black")
dev.off()

# PCA plot by Group
png("PCA_group.png", width = 1200, height = 1000)
targets$Group <- as.factor(targets$Group)
palette(c("#E50068", "#0067E6"))
xlim <- range(pca_results$x[, 1])
ylim <- range(pca_results$x[, 2])

# Increase right margin, reduce top/bottom margins
par(mar = c(5, 4, 4, 10), xpd = TRUE)
plot(pca_results$x[, 1], pca_results$x[, 2],
     cex = 1.2, pch = 19, col = targets$Group,
     xlab = "PC1", ylab = "PC2",
     xlim = xlim, ylim = ylim,
     main = "PCA (Groups)")
text(pca_results$x[, 1], pca_results$x[, 2],
     labels = rownames(pca_results$x), 
     cex = 0.7, pos = 1, offset = 0.3, col = "black")
legend("right", inset = c(-0.05, 0),
       legend = levels(targets$Group),
       col = c(1, 2), pch = 19,
       cex = 0.9, pt.cex = 1.2, x.intersp = 1.2,
       box.lty = 0, bg = "white")
dev.off()

# PCA plot by Sex
png("PCA_sex.png", width = 1200, height = 1000)
targets$Sex <- as.factor(targets$Sex)
palette(c("#C364CA", "#8BC4F9"))
par(mar = c(5, 4, 4, 10), xpd = TRUE)
plot(pca_results$x[, 1], pca_results$x[, 2],
     cex = 1.2, pch = 19, col = targets$Sex,
     xlab = "PC1", ylab = "PC2",
     xlim = xlim, ylim = ylim,
     main = "PCA (Sex)")
text(pca_results$x[, 1], pca_results$x[, 2],
     labels = rownames(pca_results$x), 
     cex = 0.7, pos = 1, offset = 0.3, col = "black")
legend("right", inset = c(-0.05, 0),
       legend = levels(targets$Sex),
       col = 1:nlevels(targets$Sex), pch = 19,
       cex = 0.9, pt.cex = 1.2, x.intersp = 1.2,
       box.lty = 0, bg = "white")
dev.off()

# PCA plot by Slide (Batch)
png("PCA_batch.png", width = 1200, height = 1000)
targets$Slide <- as.factor(targets$Slide)
palette(c("#9e0059", "green"))
par(mar = c(5, 4, 4, 10), xpd = TRUE)
plot(pca_results$x[, 1], pca_results$x[, 2],
     cex = 1.2, pch = 19, col = targets$Slide,
     xlab = "PC1", ylab = "PC2",
     xlim = xlim, ylim = ylim,
     main = "PCA (Batch)")
text(pca_results$x[, 1], pca_results$x[, 2],
     labels = rownames(pca_results$x), 
     cex = 0.7, pos = 1, offset = 0.3, col = "black")
legend("right", inset = c(-0.05, 0),
       legend = levels(targets$Slide),
       col = c(1, 2), pch = 19,
       cex = 0.9, pt.cex = 1.2, x.intersp = 1.2,
       box.lty = 0, bg = "white")
dev.off()

#---------------------
# 9. Using the matrix of normalized beta values generated in step 7, identify differentially 
#     methylated probes between group CTRL and group DIS using the function assigned to each group.

# Define the t-test function: compares beta values between groups
My_t_test <- function(x) {
  t_test <- t.test(x ~ targets$Group)
  return(t_test$p.value)
}

# Apply the t-test function to each probe (row)
p_values <- apply(beta_norm, 1, My_t_test)

# Alternative faster method using rowttests from genefilter package
if (!requireNamespace("genefilter", quietly = TRUE)) {
  BiocManager::install("genefilter")
}
library(genefilter)
# Create factor for the groups 
group_factor <- factor(targets$Group)
# Apply vectorized t-test  
t_test_results <- rowttests(beta_norm, group_factor)
# Extract only the p-values
p_values <- t_test_results$p.value

# Create a data frame combining beta values with the raw p-values
final_ttest <- data.frame(beta_norm, p_values_ttest = p_values)

# Sort the probes by p-value in ascending order
final_ttest <- final_ttest[order(final_ttest$p_values_ttest), ]

# Filter significant probes with p-value <= 0.05
final_ttest_th <- final_ttest[final_ttest$p_values_ttest <= 0.05, ]

# Plot histogram of p-value distribution
png("Histogram_pvalues.png", width = 1200, height = 1000)
hist(final_ttest$p_values_ttest, main = "P-value distribution (t-test)", 
     xlab = 'p-value')
abline(v = 0.05, col = "#ef233c")
dev.off()

# Check the number of differentially methylated probes before correction
dim(final_ttest_th)

#---------------------
# 10. Apply multiple test correction and set a significant threshold of 0.05.
#     How many probes do you identify as differentially methylated considering nominal pValues? 
#     How many after Bonferroni correction? How many afteR BH correction?

# Apply multiple testing corrections
corr_pval_BH <- p.adjust(final_ttest$p_values_ttest, method = "BH")  # Benjamini-Hochberg
corr_pval_bonferroni <- p.adjust(final_ttest$p_values_ttest, method = "bonferroni")

# Combine the corrected p-values into the final data frame
final_ttest_corr <- data.frame(final_ttest, corr_pval_BH, corr_pval_bonferroni)

# Count the number of differentially methylated probes before and after correction
before_correction <- nrow(final_ttest_corr[final_ttest_corr$p_values_ttest <= 0.05, ])
after_Bonferroni <- nrow(final_ttest_corr[final_ttest_corr$corr_pval_bonferroni <= 0.05, ])
after_BH <- nrow(final_ttest_corr[final_ttest_corr$corr_pval_BH <= 0.05, ])

# Create a summary data frame with the results
diff_meth_df <- data.frame(before_correction, after_Bonferroni, after_BH)
rownames(diff_meth_df) <- c("# differentially methylated probes")

# Print the summary data frame
print(diff_meth_df)

# Set plotting area to a single plot
par(mfrow = c(1,1))

# Plot boxplots of raw and adjusted p-values
png("Boxplot_corrections.png", width = 1200, height = 1000)
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
dev.off()

#---------------------
# 11. Produce a volcano plot and a Manhattan plot of the results of differential methylation analysis

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
png("Volcano_plot.png", width = 1200, height = 1000)
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
dev.off()

# Manhattan Plot

# Add probe IDs as a column
final_ttest_corr_anno <- data.frame(rownames(final_ttest_corr), final_ttest_corr)
colnames(final_ttest_corr_anno)[1] <- "IlmnID"

# Merge with Illumina manifest file (must contain IlmnID, CHR, MAPINFO)
final_ttest_corr_anno <- merge(final_ttest_corr_anno, Illumina450Manifest_clean, by = "IlmnID")

# Prepare the data frame for qqman
input_Manhattan <- data.frame(
  ID = final_ttest_corr_anno$IlmnID,
  CHR = final_ttest_corr_anno$CHR,
  MAPINFO = final_ttest_corr_anno$MAPINFO,
  PVAL = final_ttest_corr_anno$p_values_ttest
)

# Convert chromosomes X and Y to numeric (23 and 24)
levels(input_Manhattan$CHR)[levels(input_Manhattan$CHR) == "X"] <- "23"
levels(input_Manhattan$CHR)[levels(input_Manhattan$CHR) == "Y"] <- "24"
input_Manhattan$CHR <- as.numeric(as.character(input_Manhattan$CHR))

# Load the qqman library and define a color palette
install.packages("qqman")
library(qqman)
color_palette <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2",
  "#7F7F7F", "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896",
  "#C5B0D5", "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"
)

# Plot the Manhattan plot
png("Manhattan_plot.png", width = 600, height = 500)
par(cex = 1.3)
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
dev.off()


#---------------------
# 12. Produce an heatmap of the top 100 significant, differentially methylated probes.

# Load required library
if (!require("gplots")) install.packages("gplots")
library(gplots)

# Prepare the input matrix for the heatmap
input_heatmap <- as.matrix(final_ttest[1:100, 1:8])

# Assign colors to groups
# RESET group_color
group_color <- c()
for (i in 1:ncol(input_heatmap)) {
  name <- colnames(input_heatmap)[i]
  sample_id <- strsplit(name, "_")[[1]][1]
  matching_row <- which(targets$SampleID == sample_id)
  if (length(matching_row) > 0) {
    if (targets$Group[matching_row] == "CTRL") {
      group_color <- c(group_color, "#CBCBCB")
    } else {
      group_color <- c(group_color, "#393939")
    }
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
