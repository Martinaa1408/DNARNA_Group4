library(minfi)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

rm(list=ls())
setwd("C:/Users/coric/Desktop/1 YEAR/2 SEMESTER/DNA-RNA/PROGETTO/Final_Report-20250526/")
library(minfi)
# Set the directory in which the raw data are stored and load the samplesheet 
# using the function read.metharray.sheet
baseDir <- ("C:/Users/coric/Desktop/1 YEAR/2 SEMESTER/DNA-RNA/PROGETTO/Final_Report-20250526/")
targets <- read.metharray.sheet(baseDir)
# Create an object of class RGChannelSet using the function read.metharray.exp
RGset <- read.metharray.exp(targets = targets)
save(RGset,file="RGset.RData")

# We extract the Green and Red Channels using the functions getGreen and getRed
Red <- data.frame(getRed(RGset))
dim(Red) # rows: probes, columns: samples
head(Red)
Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)

# Address of interest
address <- "44666390"

# Get Red and Green fluorescence values for the specified address
probe_red <- Red[address, ]
probe_green <- Green[address, ]

# Load the cleaned manifest file to check the type of the probe
load('Illumina450Manifest_clean.RData')
# check whether it is of type I or II from column 'Infinium_Design_Type'
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="44666390",'Infinium_Design_Type']

# extract all data related to sample, red fluorescence, green fluorescence and type and color to fill the datable
df_address <- data.frame(Sample=colnames(probe_green), Red_fluor=unlist(probe_red, use.names = FALSE ), Green_fluor=unlist(probe_green, use.names = FALSE), Type = "II")
df_address

Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="44666390",'Infinium_Design_Type']
df_address <- data.frame(Sample=colnames(probe_green), Red_fluor=unlist(probe_red, use.names = FALSE ), Green_fluor=unlist(probe_green, use.names = FALSE), Type = "I")
df_address
