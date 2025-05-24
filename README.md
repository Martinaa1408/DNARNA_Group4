# Group 4 - DNA Methylation Analysis Project


## Table of Contents
- [Project Overview](#project-overview)
- [Assigned Parameters](#assigned-parameters)
- [Tools and Technologies](#tools-and-technologies)
- [Workflow Summary](#workflow-summary)
- [Motivation and Rationale](#motivation-and-rationale)
- [Outputs and Deliverables](#outputs-and-deliverables)
- [Academic Context](#academic-context)

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
- Use `plotQC()`, `controlStripPlot()`  
- Calculate detection p-values and count failed probes

### 5. Beta and M Values
- Calculate and compare raw beta and M values for WT and MUT  
- Plot density distributions

### 6. Normalization
- Apply `preprocessFunnorm()`  
- Compare raw vs normalized data using 6-panel plot and boxplots

### 7. PCA Analysis
- Perform PCA on normalized data  
- Plot clustering by group, sex, and batch

### 8. Differential Methylation Analysis
- t-test per probe, p-value correction (BH and Bonferroni)  
- Visualize with histogram, volcano plot, and Manhattan plot

### 9. Heatmap
- Plot top 100 most significant differentially methylated probes

---

## Motivation and Rationale

This pipeline follows best practices in methylation data analysis:
- Ensures **quality control** through detection p-values and QC plots
- Applies **Funnorm normalization**, which adjusts for Type I/II probe differences and technical variation
- Performs **statistical testing** to detect biologically relevant methylation changes

The `preprocessFunnorm` function is particularly well suited for correcting bias without eliminating biological signal in heterogeneous datasets like ours.

---

## Outputs and Deliverables

- Commented R script (`pipeline_group4.R`)
- Table with fluorescence intensities for probe `44666390`
- Quality control metrics and visualizations
- 6-panel comparison plots (raw vs normalized)
- PCA plots grouped by biological variables
- Volcano and Manhattan plots of differential analysis
- Heatmap of top 100 probes
- Summary tables for differentially methylated probes (raw, BH, Bonferroni)

---

## Academic Context

This project was developed within the course **DNA/RNA dynamics (Module 2, Prof Ravaioli)** of the MSc in **Bioinformatics** at the **University of Bologna**.

The aim of the course is to provide hands-on experience with the analysis of DNA methylation and gene expression using real experimental datasets and R/Bioconductor tools.

---

> This README provides the methodological foundation for our group project and supports the creation of the final report and codebase.
