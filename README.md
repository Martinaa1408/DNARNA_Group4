# Group 4 – DNA Methylation Analysis Project

## Table of Contents
- [Project Overview](#project-overview)
- [Group Members](#group-members)
- [Assigned Parameters](#assigned-parameters)
- [Tools and Technologies](#tools-and-technologies)
- [Repository Structure](#repository-structure)
- [Workflow Summary](#workflow-summary)
- [Motivation and Rationale](#motivation-and-rationale)
- [Resources and References](#resources-and-references)

---

## Project Overview
This repository contains the DNA methylation analysis performed by **Group 4** for the *DNA/RNA Dynamics* course (Module 2, Prof. Francesco Ravaioli), MSc in Bioinformatics – University of Bologna. The project aims to identify differentially methylated positions (DMPs) between **control (CTRL)** and **disease (DIS)** samples using the Illumina HumanMethylation450K array and the `minfi` R package.

---

## Group Members

- Andrea Pusiol — andrea.pusiol@studenti.unibo.it
- Aurora Mazzoni — aurora.mazzoni2@studenti.unibo.it
- Martina Castellucci — martina.castellucci@studenti.unibo.it
- Alessia Corica — alessia.corica@studenti.unibo.it
- Sofia Natale — sofia.natale@studenti.unibo.it
- Bianca Mastroddi — bianca.mastroddi@studenti.unibo.it
- Perla Lucaboni — perla.lucaboni@studenti.unibo.it

---

## Assigned Parameters

| Parameter                 | Value          |
|---------------------------|----------------|
| Group ID                  | 4              |
| Probe Address             | 44666390       |
| Detection p-value cut-off | 0.01           |
| Normalization method      | preprocessFunnorm |

---

## Tools and Technologies

- **Language**: R  
- **Platform**: Illumina HumanMethylation450K  
- **Packages**: `minfi`, `BiocManager`, `gplots`, `factoextra`, `qqman`

---

## Repository Structure

### `/Input_data/` — Raw Data & Metadata
- `*.idat` (16 files): Red and green channel raw data files containing probe intensities per sample.
- `SampleSheet_Report_II.csv`: Sample metadata including IDs, experimental groups, batch information, and corresponding .idat file references.


### `/scripts/` — Main Pipeline & Report
- `pipeline_group4.R`: Main R script executing the entire workflow.
- `report.html`: Auto-generated interactive HTML report summarizing results with plots and tables.


### `/outputs/` — Results & Visualizations

#### Quality Control & Intensity
- `Density_plot.png`: Density distribution of beta values or intensities.
- `qc_plot_msetraw.png`: Raw MSet data distribution before normalization.
- `raw_vs_normalized.plot.png`: Scatter plot comparing raw vs normalized values.
- `negative_control_intensity_check.png`: Background control plot using negative probes.
- `df_address.pdf`: Sample plate addresses for QC.
- `df_failed.pdf`: Summary of failed or excluded samples.

#### PCA 
- `PCA_batch_plot.png`: PCA plot colored by batch.
- `PCA_groups_plot.png`: PCA plot colored by experimental group (e.g., CTRL/DIS).
- `PCA_sex_plot.png`: PCA plot colored by sex (Female/Male).
- `scree_plot.png`: Scree plot showing explained variance per principal component.

#### Clustering
- `Average_linkage_heatmap.png`: Heatmap with average linkage clustering.
- `Complete_linkage_heatmap.png`: Heatmap with complete linkage clustering.
- `Single_linkage_heatmap.png`: Heatmap with single linkage clustering.

#### Statistical Analysis
- `p-value_distribution_plot.png`: Histogram of raw p-values (t-tests).
- `p-value_distribution_raw_adjusted_plot.png`: Comparison of raw and adjusted p-values (BH, Bonferroni).
- `manhattan_plot.png`: Manhattan plot of –log₁₀ p-values across genomic positions.
- `volcano_plot.png`: Volcano plot of ΔBeta vs –log₁₀ p-value (effect size vs significance).

#### Fluorescence Tables
- `Green_Fluorescences_Table.yaml`: Green channel intensities by sample.
- `Red_Fluorescences_Table.yaml`: Red channel intensities by sample.


### `/diagram_workflow/` — Workflow Overview
- `workflow.png`: Diagram illustrating the main steps of the DNA methylation analysis pipeline, from raw data input to final output and visualization.


### `/report_pipeline/` — Report Guidelines
- `20250605_Report_pipeline_FINAL.pdf`: Document containing the professor’s official instructions and structure for writing the final project report.


`DNAmethylation_analysis_manual.pdf`: Manual (LaTeX) including R usage guide, function descriptions, package references, and pipeline instructions.

**Download processed RGset object**:  
[RGset.RData](https://drive.google.com/uc?export=download&id=1eIU1pHnwIDmMTmn73Zu3RdZgdcb_ZFux)

---

## Workflow Summary

image..

---

## Motivation and Rationale

This workflow ensures:
- Rigorous probe filtering using detection p-values.
- Correction for probe-type bias with Funnorm.
- Reproducible identification of biologically meaningful DMPs.
- Comprehensive visualization of results to facilitate interpretation.

---

## Resources and References

- #%% md
#### Installing packages – brief description and references

- **BiocManager**  
  Used to install and manage packages from the Bioconductor project.  
  📖 Vignette: [BiocManager Vignette](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)  
  📚 Citation:  
  Shepherd L. (2024). *BiocManager: Access the Bioconductor Project Package Repository*. R package version 1.30.22.

- **minfi**  
  Core package for analyzing Illumina 450K/EPIC methylation arrays. Includes preprocessing, QC, DMP analysis, and visualization tools.  
  Vignette: [minfi Vignette](https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html)  
  Citation: Aryee, M. J., et al. (2014). *Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation              microarrays*. Bioinformatics, 30(10), 1363–1369. DOI: [10.1093/bioinformatics/btu049](https://doi.org/10.1093/bioinformatics/btu049)

- **Illumina 450K Manifest**  
  CSV annotation file for probe IDs, positions, and types on the 450K array.  
  Documentation: [Illumina Manifest Column Headings](https://support.illumina.com/bulletins/2016/05/infinium-methylationk-manifest-column-headings.html)

- **factoextra**  
  Simplifies extraction and visualization of multivariate analyses (e.g., PCA).  
  Documentation: [factoextra Website](https://rpkgs.datanovia.com/factoextra/index.html)  
  Citation:  
  Kassambara, A. (2020). *factoextra: Extract and Visualize the Results of Multivariate Data Analyses*. R package version 1.0.7.

- **qqman**  
  Produces Manhattan and Q-Q plots, primarily used in GWAS and EWAS visualizations.  
  Vignette: [qqman Vignette](https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html)  
  Citation:  
  Turner, S. D. (2014). *qqman: an R package for visualizing GWAS results using Q-Q and Manhattan plots*. bioRxiv. DOI: [10.1101/005165]                     (https://doi.org/10.1101/005165)

- **gplots**  
  Offers plotting tools including `heatmap.2()` for hierarchical clustering and visualizations.  
  CRAN Page: [gplots on CRAN](https://cran.r-project.org/web/packages/gplots/index.html)  
  Citation:  
  Warnes, G. R., et al. (2022). *gplots: Various R Programming Tools for Plotting Data*. R package version 3.1.3.

- [Illumina 450K Product Files](https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html)  
  Official documentation and downloads: manifest files, annotation files, and sample sheets.

- [Infinium HumanMethylation450 BeadChip – Datasheet (PDF)](https://www.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/datasheet_humanmethylation450.pdf)  
  Detailed technical overview: probe design, detection chemistry, and array performance.

- [minfi::preprocessFunnorm() function](https://www.rdocumentation.org/packages/minfi/versions/1.18.4/topics/preprocessFunnorm)
  Official documentation describing how to apply functional normalization to Illumina HumanMethylation arrays.

---

> _This repository documents a reproducible methylation analysis workflow combining theoretical insights and practical bioinformatics skills, tailored for CTRL vs DIS comparisons._

