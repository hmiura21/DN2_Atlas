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


---

## 🔐 Data Privacy Notice

This project uses potential **confidential raw data** from unpublished or protected datasets and **does not** include them in this repository.

➡️ To test the scripts, please substitute with your own pre-processed single-cell B cell expression matrix.

