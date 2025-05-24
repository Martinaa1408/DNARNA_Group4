# Group 4 - DNA Methylation Analysis Project

## 📑 Table of Contents
- [Project Overview](#-project-overview)
- [Assigned Parameters](#-assigned-parameters)
- [Tools and Technologies](#-tools-and-technologies)
- [Workflow Summary](#-workflow-summary)
- [Motivation and Rationale](#-motivation-and-rationale)
- [Outputs and Deliverables](#-outputs-and-deliverables)
- [Academic Context](#-academic-context)

---

## Project Overview
This project involves the analysis of DNA methylation data using the Illumina Infinium BeadChip microarray platform. Our aim is to investigate potential epigenetic differences between wild-type (WT) and mutant (MUT) samples through a complete preprocessing, quality control, normalization, and statistical analysis pipeline.

**Group Members:**
- Andrea Pusiol  
- Aurora Mazzoni  
- Martina Castellucci  
- Alessia Corica  
- Sofia Natale  
- Bianca Mastroddi  
- Perla Lucaboni  

---

## Assigned Parameters

| Parameter                  | Value              |
|---------------------------|--------------------|
| Group ID                  | 4                  |
| Assigned Probe Address    | 44666390           |
| Detection p-value cutoff  | 0.01               |
| Normalization Method      | preprocessFunnorm  |

---

## Tools and Technologies

- **Language**: R  
- **Key packages**: `minfi`, `BiocManager`, `ggplot2`, `factoextra`, `gplots`, `qqman`  
- **Platform**: Illumina Infinium 450K BeadChip  

---

## Workflow Summary

### 1. Data Import
- Load IDAT files via `read.metharray.exp()` and create `RGset`

### 2. Fluorescence Extraction
- Extract red and green signals for address `44666390` and determine probe type

### 3. Create MSet.raw
- Use `preprocessRaw()` to convert RGset into MSet.raw

### 4. Quality Control
- `plotQC()`, `controlStripPlot()`  
- Detection p-values and table of failed probes

### 5. Beta and M Values
- Calculate and compare beta and M values for WT and MUT  
- Density plots for raw values

### 6. Normalization
- Use `preprocessFunnorm()`  
- Visualize changes via 6-panel plots and boxplots

### 7. PCA Analysis
- Perform PCA on normalized beta matrix  
- Plot by group, sex, and batch

### 8. Differential Methylation Analysis
- t-test per probe + BH/Bonferroni correction  
- Histogram, volcano and Manhattan plots

### 9. Heatmap
- Top 100 most significant probes

---

## Motivation and Rationale

This pipeline follows best practices for methylation microarray analysis:
- **QC** ensures sample quality
- **Funnorm** adjusts for probe-type and technical biases
- **Differential analysis** uncovers methylation changes linked to biological conditions

The choice of `preprocessFunnorm` is optimal for this dataset as it corrects known biases while preserving biological signal.

---

## Outputs and Deliverables

- Annotated `.R` script  
- Red/Green intensity table for probe `44666390`  
- All plots (6-panel, PCA, volcano, Manhattan, heatmap)  
- Summary of differentially methylated probes before and after correction  

---

## Academic Context

This project was developed as part of the **DNA/RNA course – Module 2** (**Prof. Francesco Ravaioli**) of the MSc in Bioinformatics at the **University of Bologna**.  
It applies statistical and computational methods to analyze DNA methylation microarray data using real experimental data.

---
> This README explains our methodology and decisions, guiding the project development and report writing.
