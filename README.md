# Computational analysis of gene-expression asthma biomarkers: predicting severity, cellular specificity and drug identification

*reproduce all analyses using this codebase*

---

## Table of Contents

- [Overview](#overview)  
- [Getting Started](#getting-started)  
- [Running the Scripts](#running-the-scripts)  
- [How to Cite](#how-to-cite)  

---

## Overview

This repository contains scripts for the:

1. identification of blood-based mRNA biomarkers of the late-phase 
asthmatic response (LAR-mRNA) using the NanoString PanCancer platform
2. enrichment analysis of PanCancer and existing mRNA biomarkers
3. evaluation of LAR-mRNA biomarkers in predicting asthma exacerbations and
severity
4. determine cellular specificity of the LAR-mRNA biomarkers
5. identifying drugs that can reverse allergen-induced exprsesion of LAR-mRNA
biomarkers

---

## Getting Started

### 1) clone repo

```
git clone https://github.com/CompBio-Lab/asthma_LAR-mRNA_biomarkers.git
```

### 2) get data
> download data from this [zenodo link]()
> put data in asthma_LAR-mRNA_biomarkers folder

### 3) open .Rproj file
> asthma_latephase.Rproj

### 4) Install renv if you don't have it:

```
install.packages("renv")
```

### 5) Restore the project environment from renv.lock:

```
renv::restore()
```
> This will install the exact package versions used when this project was last saved.


## Running scripts

- Open .Rmd files in RStudio in the src folder

- Knit the files to HTML or run interactively


## How to Cite

- [preprint]()