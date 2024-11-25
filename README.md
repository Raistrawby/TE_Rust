# TE_Rust
TE comparative genomic between Pucciniales and Other-Pucciniomycotina
This repository is link to the publication: XXXX  
Please cite us.

## **TRansposable Elements eXploration (TREX)** 
You will find in this repository all the script to large-scale genomes TE analysis, with curations steps provided.
This repository does not includes the following main tools: **REPET**, **MCHELPER**, and **CD-HIT**.

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [REPET](#repet)
  - [MCHELPER](#mchelper)
  - [CD-HIT](#cd-hit)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

---

## Introduction

This toolkit is designed for comprehensive genome analysis with a particular focus on repetitive element annotation. Whether you're working with small or large genomes, haploid or polyploid genomes, this toolkit offers tailored approaches to suit your analysis needs.

## Features

- **REPET**: Annotation and detection of repetitive elements in genomes. Please refer to : https://urgi.versailles.inra.fr/Tools/REPET

- **MCHELPER**: Help for TE curation an reclassify Conflict TEs. Please refer to : https://github.com/GonzalezLab/MCHelper/tree/main
- **CD-HIT**: Use for TE clustering and to classify the Unclassified. Please refer to: https://sites.google.com/view/cd-hit/home?authuser=0 
---

## Installation

To get started, clone the repository.

```bash
git clone https://github.com/Raistrawby/TE_Rust.git

cd TE_Rust
```

## PostRepet File:

First you will need to create the "big" dataframe of all species with supplementale informations for further analysis.

```bash
python ConfigDFRepet.py
```

## 1) Curation Step

Reclassify based on the highest % or in the percent from the threshold (50%).

``` bash
python ClassifyConflict.py 
```

## 2) Curation Step with MCHelper and CDHIT

``` bash
python Curation.py -R fileRepet.csv -M fileMCHelper -C clusterfileCDHIT
```

