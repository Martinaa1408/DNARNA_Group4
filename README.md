# 🧬 Group 4 – DNA Methylation Analysis Project

## 📑 Table of Contents
- [Project Overview](#project-overview)
- [Group Members](#group-members)
- [Assigned Parameters](#assigned-parameters)
- [Tools and Technologies](#tools-and-technologies)
- [Workflow Summary](#workflow-summary)
- [Motivation and Rationale](#motivation-and-rationale)
- [Outputs and Deliverables](#outputs-and-deliverables)
- [Academic Context](#academic-context)

---

## 📌 Project Overview

This project involves the analysis of DNA methylation data generated using the **Illumina Infinium HumanMethylation450K BeadChip** platform.  
Our goal is to investigate potential **epigenetic differences** between **wild-type (WT)** and **mutant (MUT)** samples through a complete pipeline including:

- Data preprocessing  
- Quality control  
- Functional normalization  
- Statistical analysis and visualization  

---

## 👥 Group Members

- Andrea Pusiol  
- Aurora Mazzoni  
- Martina Castellucci  
- Alessia Corica  
- Sofia Natale  
- Bianca Mastroddi  
- Perla Lucaboni

---

## 📋 Assigned Parameters

| Parameter                 | Value          |
|---------------------------|----------------|
| Group ID                  | 4              |
| Assigned Probe Address    | 44666390       |
| Detection p-value cutoff  | 0.01           |
| Normalization Method      | `preprocessFunnorm` |

---

## 🛠️ Tools and Technologies

- **Language**: R  
- **Platform**: Illumina HumanMethylation450K BeadChip  
- **Key packages**:  
  - `minfi`  
  - `BiocManager`  
  - `ggplot2`  
  - `factoextra`  
  - `gplots`  
  - `qqman`  

> 📁 **Input Data**: IDAT files and `SampleSheet.csv` are stored in the `input_data/` directory.

---

## 🔬 Workflow Summary

1. **Data Import**  
   Load IDAT files using `read.metharray.exp()` and generate the `RGset` object.

2. **Fluorescence Extraction**  
   Extract red and green channel signals for probe `44666390` and determine the probe type.

3. **Create MSet.raw**  
   Use `preprocessRaw()` to compute raw β-values and M-values.

4. **Quality Control**  
   Visualize QC metrics with `plotQC()` and `controlStripPlot()`, calculate detection p-values, and count failed probes.

5. **Beta and M Values**  
   Compare raw β-values and M-values between WT and MUT groups; plot density distributions.

6. **Normalization**  
   Apply `preprocessFunnorm()` and compare raw vs normalized values using 6-panel plots and boxplots.

7. **PCA Analysis**  
   Perform PCA on normalized data and plot sample clustering by group, sex, and batch.

8. **Differential Methylation Analysis**  
   Apply a t-test per probe, correct p-values using Benjamini-Hochberg and Bonferroni methods, and visualize results using histograms, volcano plots, and Manhattan plots.

9. **Heatmap**  
   Visualize the top 100 most significant differentially methylated probes with hierarchical clustering.

---

## 🧠 Motivation and Rationale

This pipeline follows best practices in methylation data analysis:

- Ensures data reliability through detection p-values and QC plots  
- Applies functional normalization (`preprocessFunnorm`) to correct for probe design bias and technical variation  
- Performs rigorous statistical testing to identify biologically meaningful methylation changes  

The use of `preprocessFunnorm` is particularly suited for heterogeneous datasets, as it corrects technical biases without removing biological signal.

---

## 📤 Outputs and Deliverables

- ✅ Commented R script (`pipeline_group4.R`)  
- ✅ Fluorescence intensity table for probe `44666390`  
- ✅ Quality control plots and detection p-value statistics  
- ✅ Raw vs normalized β-values comparison plots  
- ✅ PCA plots grouped by biological variables  
- ✅ Volcano and Manhattan plots from DMP analysis  
- ✅ Heatmap of the top 100 DMPs  
- ✅ Summary tables with raw, BH-adjusted, and Bonferroni p-values  

---

## 🎓 Academic Context

This project was developed within the course **DNA/RNA Dynamics (Module 2)**, taught by **Prof. Ravaioli**, as part of the **MSc in Bioinformatics** at the **University of Bologna**.

The course provides hands-on experience in the analysis of DNA methylation and gene expression using real experimental datasets and **R/Bioconductor** tools.  
This README serves as a methodological reference for the project and supports the development of the final report and analysis scripts.

