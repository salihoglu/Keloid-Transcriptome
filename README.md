<div align="center">

# Keloid Transcriptome Meta-analysis

### A quality-stratified multi-cohort transcriptomic workflow for keloid disease-signature discovery

![Manuscript](https://img.shields.io/badge/manuscript-submitted-blue)
![Analysis](https://img.shields.io/badge/analysis-reproducible-brightgreen)
![Data](https://img.shields.io/badge/data-GEO-orange)
![Status](https://img.shields.io/badge/status-post--submission-lightgrey)

</div>

---

## About this repository

This repository contains the reproducible analysis workflow accompanying the submitted manuscript:

> **Multi-cohort transcriptomic meta-analysis identifies a fibrosis-centered keloid disease signature and mechanism-informed exploratory compound hypotheses**

The manuscript has been submitted for journal consideration.

This project re-analyses public keloid transcriptomic datasets from the Gene Expression Omnibus using a conservative, quality-stratified framework. The workflow separates accession retrieval, differential-expression computability, gene-symbol interpretability, and strict QC-passing biological evidence to reduce overinterpretation of heterogeneous public transcriptomic resources.

---

## Study design

| Evidence layer | Description | Count |
|---|---:|---:|
| GEO accession screening | Retrieved or attempted keloid-related GEO records | 13 accessions |
| DEG-complete landscape | Accessions producing complete differential-expression reports | 7 reports |
| Strict HGNC-like consensus layer | Gene-symbol-resolved contributors used for consensus integration | 4 contributors |
| Primary interpreted evidence layer | Strict QC-passing cohorts used for biological interpretation | 3 cohorts |

The strict interpreted layer included **13 modelled keloid samples** and **12 modelled control/comparator samples**, with effective balanced support of 12 samples per group.

---

## Main finding

The primary biological signal is a recurrent **keloid-up extracellular-matrix and fibrosis-centered transcriptomic program**.

The prioritized anchor genes are:

```text
SPARC, COL1A1, COL1A2, COL3A1, TGFBI, LOXL2, COL5A1,
FN1, POSTN, COL5A2, BGN, SERPINH1, COL6A1, TGFB2
