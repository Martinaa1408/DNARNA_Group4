# Group 4 – DNA Methylation Analysis Project

## Table of Contents
- [Project Overview](#project-overview)
- [Assigned Parameters](#assigned-parameters)
- [Tools and Technologies](#tools-and-technologies)
- [Repository Structure](#repository-structure)
- [Workflow Summary](#workflow-summary)
- [Resources and References](#resources-and-references)
- [License](#license)
- [Contact](#contact)
  
---

## Project Overview
This repository contains the DNA methylation analysis developed by **Group 4** for the *DNA/RNA Dynamics* course (Module 2, Prof. Francesco Ravaioli), MSc in Bioinformatics – University of Bologna.

The project investigates genome-wide CpG methylation patterns using data from the Illumina HumanMethylation450K BeadChip, with the goal of identifying methylation changes associated with disease. The analysis, entirely performed in R using Bioconductor packages, integrates statistical testing, quality control, and biological interpretation.

Designed as both a scientific case study and an educational exercise, the repository provides a reproducible framework for exploring epigenetic variation and its potential role in disease mechanisms.

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
- **Packages**: `minfi`, `BiocManager`, `gplots`, `factoextra`, `qqman`, `genefilter`, `ggplot2`, `ggpubr`, `cluster`, `factoMineR`


---

## Repository Structure

### `/Input_data/` — Raw Data & Metadata
- `*.idat` (16 files): Red and green channel raw data files containing probe intensities per sample.
- `SampleSheet_Report_II.csv`: Sample metadata including IDs, experimental groups, batch information, and corresponding .idat file references.


### `/scripts/` — Main Pipeline & Report
- `DRD_project_script.R`: Standalone R script with the core analysis code, separated from the RMarkdown report. Useful for re-running the analysis or integrating it into other workflows.
- `DRD_project_final.html`: HTML version of the report generated from the .Rmd file. Allows for interactive viewing in a web browser.
- `DRD_project_final.Rmd`: RMarkdown source file containing the full project report, including code, results, and explanations. Can be compiled into HTML or PDF format.
- `DRD_project_final.pdf`: PDF version of the final report. Useful for printing or sharing as a static document.
  

### `/outputs/` — Results & Visualizations

#### QC & Intensity
- `beta_m_values.png`: Distribution of Beta and M values (CTRL vs DIS).
- `qc_plot.png`: Raw MSet data distribution before normalization.
- `Raw_normalised_beta.png`: Comparison of raw vs normalized values (mean, SD, boxplot).
- `controlStripPlot.png`: Background control plot using negative probes.
- `df_address.pdf`: Sample plate addresses for QC.
- `df_failed.pdf`: Summary of failed or excluded samples.

#### PCA
- `PCA_batch.png`: PCA plot colored by batch.
- `PCA_groups.png`: PCA plot colored by experimental group (e.g., CTRL/DIS).
- `PCA_sex.png`: PCA plot colored by sex (Female/Male).
- `scree_plot.png`: Scree plot showing explained variance per principal component.

#### Clustering
- `Average_linkage_heatmap.png`: Heatmap with average linkage clustering.
- `Complete_linkage_heatmap.png`: Heatmap with complete linkage clustering.
- `Single_linkage_heatmap.png`: Heatmap with single linkage clustering.

#### Statistics
- `Histogram_pvalues.png`: Histogram of raw p-values (t-tests).
- `Boxplot_corrections.png`: Comparison of raw and adjusted p-values (BH, Bonferroni).
- `manhattan_plot.png`: Manhattan plot of –log₁₀ p-values across genomic positions.
- `volcano_plot.png`: Volcano plot of ΔBeta vs –log₁₀ p-value (effect size vs significance).
  

### `/diagram_workflow/` — Workflow Overview
- `workflow.png`: Diagram illustrating the main steps of the DNA methylation analysis pipeline, from raw data input to final output and visualization.


### `/report_pipeline/` — Report Guidelines
- `20250605_Report_pipeline_FINAL.pdf`: Document containing the professor’s official instructions and structure for writing the final project report.


`DNAmethylation_analysis_manual.pdf`: Manual (LaTeX) including R usage guide, function descriptions, package references, and pipeline instructions.

---

## Workflow Summary

![Data Import](https://github.com/user-attachments/assets/b2a8ca31-ca1f-413a-baed-f4d2e99f264e)


---

## Resources and References

- **BiocManager**  
  Used to install and manage packages from the Bioconductor project.  
  Vignette: [BiocManager Vignette](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)  
  Citation:  
  Shepherd L. (2024). *BiocManager: Access the Bioconductor Project Package Repository*. R package version 1.30.22.
 
- **minfi**  
  Core package for analyzing Illumina 450K/EPIC methylation arrays. Includes preprocessing, QC, DMP analysis, and visualization tools.  
  Vignette: [minfi Vignette](https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html)  
  Citation: Aryee, M. J., et al. (2014). *Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation              microarrays*. Bioinformatics, 30(10), 1363–1369. DOI: [10.1093/bioinformatics/btu049](https://doi.org/10.1093/bioinformatics/btu049)

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

- **genefilter**
  Bioconductor package for high-throughput filtering and statistical testing.
  [genefilter](https://bioconductor.riken.jp/packages/3.8/bioc/html/genefilter.html)

- **Illumina 450K Product Files**  
  Official documentation and downloads: manifest files, annotation files, control probe info, and sample sheets.  
  [support.illumina.com – Product Files](https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html)

- **Infinium HumanMethylation450 BeadChip – Datasheet (PDF)**  
  Technical summary of the 450K array platform, including probe design, detection chemistry, and assay performance metrics.  
  [Download Datasheet](https://www.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/datasheet_humanmethylation450.pdf)

- **`minfi::preprocessFunnorm()`**  
  Official function documentation describing how to apply **functional normalization** to HumanMethylation arrays using internal control probes.  
  Reduces unwanted technical variation while preserving biological signals.  
  [Function Reference](https://www.rdocumentation.org/packages/minfi/versions/1.18.4/topics/preprocessFunnorm)

- **ggplot2**  
  Grammar of graphics implementation for elegant and layered data visualization.  
  [ggplot2 – tidyverse.org](https://ggplot2.tidyverse.org)

- **ggpubr**  
  Publication-ready plots built on top of `ggplot2`, with simplified syntax.  
  [ggpubr – datanovia.com](https://rpkgs.datanovia.com/ggpubr/)

- **cluster**  
  Core clustering algorithms and validation methods for statistical computing.  
  [cluster – CRAN](https://cran.r-project.org/package=cluster)

- **FactoMineR**  
  Multivariate exploratory data analysis including PCA, MCA, and CA.  
  [FactoMineR – Bioconductor](https://bioconductor.org/packages/FactoMineR)

---

## License

This project is licensed under the MIT License.  
See the [LICENSE](LICENSE) file for more details.

---

## Contact

For questions, feedback, or reproducibility concerns, feel free to reach out to the project members:

- Martina Castellucci — [martina.castellucci@studenti.unibo.it](mailto:martina.castellucci@studenti.unibo.it)  
- Alessia Corica — [alessia.corica@studenti.unibo.it](mailto:alessia.corica@studenti.unibo.it)  
- Sofia Natale — [sofia.natale@studenti.unibo.it](mailto:sofia.natale@studenti.unibo.it)  
- Andrea Pusiol — [andrea.pusiol@studenti.unibo.it](mailto:andrea.pusiol@studenti.unibo.it)  
- Perla Lucaboni — [perla.lucaboni@studenti.unibo.it](mailto:perla.lucaboni@studenti.unibo.it)  
- Aurora Mazzoni — [aurora.mazzoni2@studenti.unibo.it](mailto:aurora.mazzoni2@studenti.unibo.it)  
- Bianca Mastroddi — [bianca.mastroddi@studenti.unibo.it](mailto:bianca.mastroddi@studenti.unibo.it)

---

> _This repository documents a reproducible methylation analysis workflow combining theoretical insights and practical bioinformatics skills._

