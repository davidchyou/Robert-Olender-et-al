Association Rules and Hierarchical Clustering in R
--------------------------------------------------

**Overview**

This repository contains two R scripts implementing Association Rules Mining and Hierarchical Clustering.
- Association_Rules_Model.R: find frequent item sets in a dataset that meet a predetermined interestingness score.
- Hierarchical_Cluster_model.R: organize items in a dataset into a hierarchy of clusters based on pairwise distances.

**Prerequisites**

The following R packages are required:

Association Rules

        install.packages(c("readr", "dplyr", "tidyverse", "data.table", "ggplot2", "gridExtra", "grid", "cowplot", "arules", "reshape2", "epitools", "coxme", "forcats", "ggtext"))

Hierarchical Clustering

        install.packages(c("readr", "dplyr", "tidyverse", "klaR", "scatterplot3d", "reshape2", "ggplot2", "dendextend", "circlize"))

Note:
- For Association Rule, "arules" version 1.5-2 should be used. Further information can be found in the Association_Rules_Model script.
- R-version on which the initial analysis was performed is '4.4.1'.

**Data Requirements**

Initial analysis performed on the UK Biobank dataset. Data should be 'tidy', with each row representing a patient and each column representing a categorical feature.

**Customization**

The user may change several parameters, including but not limited to:
- Interestingness score for association rules.
- Number of distinct clusters ('k') for hierarchical clustering.

**Authorship and Contributions**

- Association Rules: Developemnt - Te-yuan Chyou, Implementation - Robert Olender.
- Hierarchical Clustering: Development - Robert Olender, Implementation - Robert Olender.
- Analysis of results: Robert Olender, Sandipan Roy, Te-yuan Chyou, Prasad Nishtala

**Project Information and Funding Body**

This analysis was completed as part of Robert Olender's PhD project, funded by the University Research Studentship Award at the University of Bath (project code: EA-PA1231). The Academic Ethics and Integrity Committee at the University of Bath has approved this project (form number: 6738). The project is supervised by Dr Prasad Nishtala (lead supervisor), and Dr Sandipan Roy (co-supervisor).
