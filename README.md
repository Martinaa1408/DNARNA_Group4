# Group 4 – DNA Methylation Analysis Project

## Table of Contents
- [Project Overview](#project-overview)
- [Group Members](#group-members)
- [Assigned Parameters](#assigned-parameters)
- [Tools and Technologies](#tools-and-technologies)
- [Data Structure](#data-structure)
- [Workflow Summary](#workflow-summary)
- [Motivation and Rationale](#motivation-and-rationale)
- [Outputs and Deliverables](#outputs-and-deliverables)
- [Academic Context](#academic-context)

---

## Project Overview

This project involves the analysis of DNA methylation data generated using the **Illumina Infinium HumanMethylation450K BeadChip** platform.  
Our goal is to investigate potential **epigenetic differences** between **wild-type (WT)** and **mutant (MUT)** samples through a complete pipeline including preprocessing, quality control, normalization, and statistical analysis.

---

## Group Members

- Andrea Pusiol  
- Aurora Mazzoni  
- Martina Castellucci  
- Alessia Corica  
- Sofia Natale  
- Bianca Mastroddi  
- Perla Lucaboni

---

## Assigned Parameters

| Parameter                 | Value          |
|---------------------------|----------------|
| **Group ID**              | 4              |
| **Assigned Probe**        | 44666390       |
| **Detection p-value**     | 0.01           |
| **Normalization method**  | `preprocessFunnorm` |

---

## Tools and Technologies

- **Language**: R  
- **Platform**: Illumina HumanMethylation450K BeadChip  
- **Key Packages**:  
  - `minfi`, `BiocManager`, `ggplot2`  
  - `gplots`, `factoextra`, `qqman`

---

## Data Structure

Processed object (`RGset.RData`) available for direct download:

📁 [Download RGset.RData from Google Drive](https://drive.google.com/uc?export=download&id=1eIU1pHnwIDmMTmn73Zu3RdZgdcb_ZFux)

> All data paths are referenced programmatically in the pipeline to ensure reproducibility.
> The `.RData` file is not stored in this repository due to GitHub file size limitations.


> All paths are referenced programmatically in the pipeline to ensure reproducibility.

---

## 🔬 Workflow Summary

1. **Data Import**  
   Load `.idat` files from `Input_data/idat/` using `read.metharray.exp()` and create an `RGset` object.

2. **Fluorescence Extraction**  
   Extract red and green signals for probe `44666390` and determine probe type (Infinium I/II).

3. **Create MSet.raw**  
   Convert RGset to MSet using `preprocessRaw()` to compute raw β-values and M-values.

4. **Quality Control**  
   Visualize QC metrics (`plotQC()`, `controlStripPlot()`), compute detection p-values, and flag failed probes.

5. **Beta and M Values**  
   Compare raw methylation levels across WT and MUT; visualize with density plots.

6. **Normalization**  
   Apply `preprocessFunnorm()` and compare pre- and post-normalization with 6-panel plots and boxplots.

7. **PCA Analysis**  
   Run PCA on normalized β-values and visualize sample clustering by group, sex, and batch.

8. **Differential Methylation Analysis**  
   Perform t-tests per probe, apply BH and Bonferroni correction, visualize results with volcano and Manhattan plots.

9. **Heatmap**  
   Plot top 100 differentially methylated probes using hierarchical clustering.

---

## Motivation and Rationale

This pipeline follows best practices in methylation analysis:

- Ensures rigorous **quality control** through detection p-values and control plots  
- Applies **functional normalization** to correct probe type bias and reduce technical variability  
- Enables robust **statistical detection** of biologically meaningful methylation changes  

`preprocessFunnorm` is particularly suitable for heterogeneous datasets like ours, balancing bias correction and signal preservation.

---

## 📤 Outputs and Deliverables

- ✅ Full analysis script: `pipeline_group4.R`  
- ✅ Fluorescence intensity table for probe `44666390`  
- ✅ Quality control metrics and detection statistics  
- ✅ Normalized vs raw beta-value comparisons  
- ✅ PCA plots grouped by sample metadata  
- ✅ Volcano and Manhattan plots for DMPs  
- ✅ Heatmap of top 100 DMPs  
- ✅ Summary tables: raw p-values, BH, Bonferroni

---

## Academic Context

This project was developed within the course **DNA/RNA Dynamics (Module 2, Prof. Ravaioli)** of the **MSc in Bioinformatics** at the **University of Bologna**.

> _This README serves as the methodological reference for our group project and supports the creation of the final report and reproducible codebase._

