# 🧬 Single-cell DN2 B Cell Atlas  
**Revealing Novel Subtypes and Demographic Trends**

## 🧠 Motivation
Double-negative (DN2) B cells are a novel population of B cells that lack IgD and CD27. They are:
- Abundant in autoimmune diseases, elderly individuals, and severe COVID-19 patients.
- Poorly studied on a large scale across different disease types.

> 🔬 **Goal:** Establish a cross-population DN2 B cell atlas to understand their immune roles in health and disease.

---

## 🎯 Aims
- Create a **DN2 B cell atlas** using healthy PBMC datasets.
- Incorporate disease-specific PBMC data from immune-related diseases.
- Assess **demographic variations** (age, sex, ancestry).
- Compare DN2 subtypes across healthy and diseased conditions.

---

## 🧪 Methods

- **Data Source:** Publicly available single-cell RNA-seq datasets from the Gene Expression Omnibus (GEO).
- **Cell Types Analyzed:** ~300,000 B cells across:
  - **Healthy individuals**
  - **Crohn's Disease**
  - **HIV**
  - **Malaria**
  - **Multiple Sclerosis**
  - **Severe COVID-19**
  - **Systemic Lupus Erythematosus (SLE)**

- **Analysis Pipeline:**
  - Meta-atlas creation of harmonized single-cell RNA-seq data.
  - UMAP for clustering.
  - Marker gene visualization (e.g., `TBX21`, `GZMB`, `TIGIT`, `IGHG`).
  - DN2 subtype classification:
    - **Cytotoxic DN2**
    - **MHC-II+ DN2**
    - **Other B cells**

---

## 📊 Results

### 1. DN2 Subtypes Identified
- **Three major DN2 populations** emerged from unsupervised clustering and marker gene expression.
- UMAP plots revealed distinct distributions between healthy and disease samples.

### 2. Disease and Demographic Trends
- **Proportions of DN2 subtypes vary** significantly across:
  - **Disease states**
  - **Age groups (Young Adult, Middle-aged, Older Adult)**
  - **Sex**
  - **Ancestry (EUR, AFR, ASI)**

- **Notable findings:**
  - Cytotoxic DN2 cells increase in severe COVID-19 and multiple sclerosis.
  - HIV and malaria samples show high levels of MHC-II+ DN2 cells.
  - DN2:other B cell ratio increases with age.

### 3. Functional Differences
- **Gene signatures enriched** in DN2 subtypes:
  - Antigen presentation
  - Type I interferon (IFN) response
  - Exhaustion
  - Hypoxia

---

## ✅ Conclusions
- **Single-cell atlas** reveals **distinct DN2 B cell subpopulations**.
- DN2s are not a single homogeneous group.
- **Disease, age, sex, and ancestry** significantly affect DN2 subtype composition.
- DN2s express diverse **immune signaling pathways** across contexts.

---

## 📁 Repository Contents
```bash
├── COVID/
│ └── DN2_Covid_Youngs_Wilk.R
├── Crohns/
│ └── DN2_Crohns_Martin.R
├── HIV/
│ └── HIV_Holla.R
├── Healthy/
│ ├── DN2 Project Combined.R
│ ├── DN2 Project Female.R
│ ├── DN2 Project Male.R
│ ├── DN2_AIDA_Combined.R
│ ├── DN2_AIDA_Female.R
│ ├── DN2_AIDA_Male.R
│ └── DN2_Jimmie_Healthy.R
├── Malaria/
│ └── Malaria_Holla.R
├── Multiple Sclerosis
│ └── DN2_MS_Friese.R
├── SLE/
│ ├── DN2_Lupus_Flynn.R
│ └── DN2_Lupus_Webb.R
├── DN2_Healthy&Disease.R
├── MergedDN2.R
├── DN2_Presentation.pdf
└── README.md
```

`DN2 Project Combined.R`: R Script for Preprocessing of combined Powell healthy dataset from DOI: 10.1126/science.abf3041 (predominantly white elderly population)

`DN2 Project Female.R`: R Script for Preprocessing of healthy female dataset from DOI: 10.1126/science.abf3041 (predominantly white elderly population)

`DN2 Project Male.R`: R Script for Preprocessing of healthy male dataset from DOI: 10.1126/science.abf3041 (predominantly white elderly population)

`DN2_AIDA_Combined.R`: R Script for Preprocessing of combined AIDA healthy dataset from Asian Immune Diversity Atlas (predominantly Asian population)

`DN2_AIDA_Female.R`: R Script for Preprocessing of healthy female dataset from Asian Immune Diversity Atlas (predominantly Asian population)

`DN2_AIDA_Male.R`: R Script for Preprocessing of healthy male dataset from Asian Immune Diversity Atlas (predominantly Asian population)

`DN2_Covid_Youngs_Wilk.R`: R Script for Preprocessing of COVID dataset from DOI: 10.1016/j.isci.2023.108507 (Wilk) and DOI: 10.1371/journal.ppat.1009804 (Youngs)

`DN2_Crohns_Martin.R`: R Script for Preprocessing of Crohns dataset from DOI: 10.1038/s41467-023-37849-3

`DN2_Jimmie_Healthy`: R Script for Preprocessing of healthy dataset from DOI: 10.1126/science.abf1970

`DN2_Lupus_Flynn.R`: R Script for Preprocessing of lupus dataset from DOI: 10.1038/s41590-020-0743-0

`DN2_MS_Friese.R`: R Script for Preprocessing of multiple sclerosis dataset from DOI: 10.1016/j.medj.2021.01.006

`HIV_Holla.R`: R Script for Preprocessing of HIV dataset from DOI: 10.1126/sciadv.abg8384

`Malaria_Holla.R`: R Script for Preprocessing of Malaria dataset from DOI: 10.1126/sciadv.abg8384

`DN2_Healthy&Disease.R` : R Script for visualization after healthy and disease dataset integration using RPCA

`MergedDN2.R`: R Script for (healthy) dataset integration and visualization

`DN2_Presentation`: PDF Presentation of research findings including main visualizations

---


## ✨ Shiny Visualization

Here is an interactive shiny app link to explore the data even more:
- [**Shiny Link**](https://github.com/hmiura21/DN2_Atlas)  


---

## 🎨 Example Visuals

![mergedDN2_RPCA_Integration](https://github.com/user-attachments/assets/bd0f2837-209e-49f3-961d-6b7027b7230b)
- healthy DN2 datasets integrated using RPCA color coded by study ID


![mergedDN2_studyID](https://github.com/user-attachments/assets/5a21a6a5-5d4b-48c1-8e92-3404f9571700)
- healthy DN2 datasets color coded by cluster ID


![PropTest_Sex](https://github.com/user-attachments/assets/31780067-4fa3-46cb-99f2-549e117c673e)
- Proportion test between male vs female sex in healthy DN2 dataset


![StackedBarPlot_DN2_Powell_Female_IndividualEffect](https://github.com/user-attachments/assets/12e2ac02-3b7b-4e30-b79c-3008937a4247)
- stacked bar plot for female Powell data visualizing individual effects

![UMAP_Lupus_Webb_Annotated](https://github.com/user-attachments/assets/2c93d3f2-554a-4e59-a276-84d4da55704c)
- UMAP cell annotations for Lupus Webb dataset

---

## 🔐 Data Privacy Notice

This project uses potential **confidential raw data** from unpublished or protected datasets and **does not** include them in this repository.

➡️ To test the scripts, please substitute with your own pre-processed single-cell B cell expression matrix.

