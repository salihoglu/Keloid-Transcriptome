# Keloid Transcriptome Meta-analysis

This repository contains the reproducible analysis workflow accompanying the submitted manuscript:

**Multi-cohort transcriptomic meta-analysis identifies a fibrosis-centered keloid disease signature and mechanism-informed exploratory compound hypotheses**

The manuscript has been submitted for journal consideration.

## Overview

This project re-analyses public keloid transcriptomic datasets from GEO using a conservative, quality-stratified framework. The analysis separates accession retrieval, differential-expression computability, gene-symbol interpretability, and strict QC-passing biological evidence to avoid over-counting public datasets as patient-level cohorts.

The workflow screened **13 GEO accessions**, generated **7 complete differential-expression reports**, integrated **4 strict HGNC-like consensus contributors**, and interpreted **3 strict QC-passing cohorts** as the primary evidence layer. The strict interpreted layer included **13 modelled keloid** and **12 modelled control/comparator** samples, with effective balanced support of 12 samples per group. :contentReference[oaicite:0]{index=0}

## Main finding

The analysis identifies a recurrent keloid-up extracellular-matrix and fibrosis-centered transcriptomic program involving:

**SPARC, COL1A1, COL1A2, COL3A1, TGFBI, LOXL2, COL5A1, FN1, POSTN, COL5A2, BGN, SERPINH1, COL6A1, and TGFB2**. :contentReference[oaicite:1]{index=1}

Functional enrichment converged on extracellular-matrix organisation, collagen biosynthesis, collagen fibril assembly, ECM proteoglycans, fibronectin matrix formation, and TGF-beta-associated pathway context. :contentReference[oaicite:2]{index=2}

## Repository contents

- GEO dataset eligibility and QC logic
- Metadata curation and phenotype assignment scripts
- Cohort-level differential-expression workflows
- Quality-weighted consensus integration
- GO, KEGG, and Reactome enrichment analyses
- STRINGdb-based network-prioritization outputs
- Sensitivity analyses
- Exploratory mechanism-informed compound prioritization

## Interpretation note

The compound candidates **pirfenidone, tranilast, doxycycline, triamcinolone, and ruxolitinib** are reported only as exploratory, mechanism-informed hypotheses. They should not be interpreted as validated therapeutics, clinical recommendations, or confirmed LINCS/CMap reversal hits. :contentReference[oaicite:3]{index=3}

## Status

**Manuscript status:** submitted for journal consideration.

This repository is intended to support transparency, reproducibility, and post-submission access to the computational workflow and analysis outputs.
