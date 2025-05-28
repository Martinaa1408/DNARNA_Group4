# Group 4 – DNA Methylation Analysis Project

## Table of Contents
- [Project Overview](#project-overview)
- [Group Members](#group-members)
- [Theoretical Background](#theoretical-background)
- [Assigned Parameters](#assigned-parameters)
- [Tools and Technologies](#tools-and-technologies)
- [Repository Structure](#repository-structure)
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

**DNA methylation** is an epigenetic modification in which a methyl group (–CH₃) is covalently added to the 5-carbon position of cytosine residues, primarily within CpG dinucleotides. This biochemical mark modulates gene expression without altering the nucleotide sequence, typically resulting in transcriptional repression when located in gene promoters or CpG islands.

**Illumina BeadChip arrays** (e.g. HumanMethylation450K) assess methylation at over 450,000 CpG sites by analyzing bisulfite-treated genomic DNA. After bisulfite conversion, each CpG is interrogated using probes designed to discriminate between methylated and unmethylated sequences:

-**Infinium I**: employs two separate probes per CpG site, each specific to either the methylated or 
 unmethylated sequence, and uses single-color fluorescence detection (either red or green depending on 
 the base).

-**Infinium II**: uses a single probe per CpG site, relying on two-color detection to distinguish between 
 methylated and unmethylated states via differential base incorporation at the single-base extension site.



After **bisulfite treatment**:

-**Unmethylated cytosines** are deaminated to **uracils**, which are then amplified as thymines during PCR.

-**Methylated cytosines** remain **unchanged**, preserving the original cytosine signal.

![Infinium I and II](https://github.com/Martinaa1408/DNARNA_Group4/blob/main/figures/infinium_scheme.png)

**Fluorescence detection** is based on labeled nucleotide incorporation at the single-base extension site:

-**Green fluorescence** (e.g., Cy3) → indicates incorporation of a base matching a methylated cytosine.

-**Red fluorescence** (e.g., Cy5) → indicates incorporation corresponding to an unmethylated cytosine.


The **methylation level** is calculated using:

- **Beta value (β)**:
  (β) = M / (M + U + 100)
  
  Proportion of methylated signal (0 = unmethylated, 1 = fully methylated). Easy to interpret but 
  compressed at extremes.

- **M-value**:
  log2((M + 1) / (U + 1))
  
  Log2 ratio of methylated vs unmethylated intensity. Preferred for statistical modeling due to better 
  distribution properties.

**Note**: M-values are approximately linear around β = 0.5, but expand better dynamic range at the extremes (β ≈ 0 or 1).


**Normalization** is crucial to remove technical biases (batch effects, probe design, dye bias). Several methods exist:

- `preprocessFunnorm`: **Functional normalization** using control probes; ideal for datasets with global 
   methylation differences (e.g. cancer).
- `preprocessQuantile`: Subset-quantile normalization for comparability across arrays.
- `preprocessNoob`: Background correction using out-of-band probes.
- `preprocessSWAN`: Adjusts for probe-type bias (Infinium I vs II).

In this project, we applied **`preprocessFunnorm`**, which is robust against heterogeneity and preserves biological differences.

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

## Repository Structure

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
A powerful R/Bioconductor package for analyzing Illumina Infinium methylation arrays. It supports     
preprocessing, normalization, quality control, and downstream statistical analysis for 450K and EPIC 
arrays.

- [Illumina 450K Product Files](https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html)  
Official page from Illumina providing downloadable content related to the 450K array: annotation files, manifest files, sample sheets, and reference genome mappings.

- [Infinium HumanMethylation450 BeadChip – Product Datasheet (PDF)](https://www.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/datasheet_humanmethylation450.pdf)  
Technical datasheet describing the design, probe chemistry, and applications of the 450K array. Includes performance metrics and hybridization details.


> _The repository documents a reproducible methylation analysis workflow combining theoretical background and practical skills._
