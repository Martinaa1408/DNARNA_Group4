# Group 4 – DNA Methylation Analysis Project

## Table of Contents
- [Project Overview](#project-overview)
- [Group Members](#group-members)
- [Theoretical Background](#theoretical-background)
- [Assigned Parameters](#assigned-parameters)
- [Tools and Technologies](#tools-and-technologies)
- [Input Data Structure](#input-data-structure)
- [Workflow Summary](#workflow-summary)
- [Motivation and Rationale](#motivation-and-rationale)
- [Outputs and Deliverables](#outputs-and-deliverables)
- [Academic Context](#academic-context)
- [Resources and References](#resources-and-references)

---

## Project Overview
This repository contains the DNA methylation analysis performed by **Group 4** for the *DNA/RNA Dynamics* course (Module 2, Prof. Francesco Ravaioli), MSc in Bioinformatics – University of Bologna. The objective is to identify differentially methylated positions (DMPs) between **wild-type (WT)** and **mutant (MUT)** samples using the Illumina HumanMethylation450K array and the `minfi` R package.

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

## Theoretical Background

**DNA methylation** is an epigenetic modification where a methyl group is added to cytosine bases, primarily at CpG dinucleotides. This mechanism regulates gene expression without altering the DNA sequence, often repressing transcription when methylation occurs in gene promoters or CpG islands.

**Illumina BeadChip arrays** quantify methylation at >450,000 CpG sites using bisulfite-converted DNA. Each CpG site is interrogated using probes that detect either the methylated or unmethylated version of the sequence:

- **Infinium I**: two probes per CpG site, single-color detection  
- **Infinium II**: one probe per site, dual-color detection

After bisulfite treatment:

- Unmethylated cytosines become uracils (then thymines after PCR)
- Methylated cytosines remain unchanged

Fluorescence detection:
- Green fluorescence → indicates methylated cytosines
- Red fluorescence → indicates unmethylated cytosines

![Infinium I and II chemistries]((https://github.com/Martinaa1408/DNARNA_Group4/blob/main/figures/infinium_scheme.png))

The methylation level is calculated using:

- **Beta value (β)**: proportion of methylated signal, ranging from 0 (unmethylated) to 1 (fully methylated)
- **M-value**: log2 ratio of methylated/unmethylated intensity

**Normalization** is essential to correct for technical and probe-design biases. We used `preprocessFunnorm`, a method optimized for heterogeneous samples.

---

## Assigned Parameters

| Parameter               | Value          |
|------------------------|----------------|
| Group ID               | 4              |
| Probe Address          | 44666390       |
| Detection p-value cut-off | 0.01         |
| Normalization method   | preprocessFunnorm |

---

## Tools and Technologies

- **Language**: R  
- **Platform**: Illumina HumanMethylation450K  
- **Packages**: `minfi`, `BiocManager`, `ggplot2`, `gplots`, `factoextra`, `qqman`

---

## Input Data Structure

/Input_data/ → .idat files (Red and Green channels) and SampleSheet.csv
/scripts/ → pipeline_group4.R
/output/ → figures, tables, PCA, volcano, heatmap

**Download RGset.RData** (processed methylation object):  
https://drive.google.com/uc?export=download&id=1eIU1pHnwIDmMTmn73Zu3RdZgdcb_ZFux

---

## Workflow Summary

- Load IDAT files using `read.metharray.exp()`
- Extract red/green signal for probe 44666390
- Create `MSet.raw` with `preprocessRaw()`
- Perform quality control (detection p-values, QC plots)
- Normalize data using `preprocessFunnorm`
- Conduct PCA for sample stratification
- Perform differential methylation analysis (t-test + BH and Bonferroni)
- Visualize results: PCA plots, volcano, Manhattan, heatmap

---

## Motivation and Rationale

This workflow ensures:
- Reliable probe filtering via detection p-values
- Correction for type I/II probe bias through Funnorm
- Reproducible identification of biologically relevant DMPs
- Clear data visualization for interpretation and reporting

---

## Outputs and Deliverables

- Annotated R script (`pipeline_group4.R`)
- Fluorescence table for probe 44666390
- Quality control metrics (plotQC, detection p)
- Raw vs normalized beta value plots
- PCA plots by group, sex, batch
- Volcano & Manhattan plots of DMPs
- Heatmap of top 100 differentially methylated probes
- Summary tables: raw p, BH-adjusted, Bonferroni

---

## Academic Context

This project was developed for the *DNA/RNA Dynamics* course (Module 2, Prof. Ravaioli) within the **MSc in Bioinformatics** program at the **University of Bologna**.

---

## Resources and References

- [`minfi` – Bioconductor package](https://bioconductor.org/packages/release/bioc/html/minfi.html)  
- [Illumina 450K Product Files](https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html)  
- [GEO Datasets – Illumina 450K](https://www.ncbi.nlm.nih.gov/geo/query/?search=450k)  
- [BeadArray Technology Overview](https://www.illumina.com/techniques/microarrays/array-data-analysis-experimental-design/beadarrays.html)

> _The repository documents a reproducible methylation analysis workflow combining theoretical background and practical skills._
