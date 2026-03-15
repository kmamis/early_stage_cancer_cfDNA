# early_stage_cancer_cfDNA
Data files and codes for paper "Early-stage cancer results in a multiplicative increase in cell-free DNA originating from healthy tissue" (2026) by Konstantinos Mamis and Ivana Bozic.

== Data files (.mat) ==

All datasets analysed in this project were previously published; no new data were generated for this work.
For each .mat file included here, we report the original publication from which the data were obtained.
The .mat files provided in this repository are derived directly from those published datasets 
and are supplied only to organise the data into analysis‑ready vectors, 
following the classifications used in the manuscript 
(e.g., grouping by cancer type, stage, origin of plasma mutations, or MAF thresholds).

cfDNA_Cohen_2018.mat

Data on cfDNA concentrations (ng/mL) for healthy individuals and stage I-III cancer patients from 
Cohen JD, Li L, Wang Y, et al. Detection and localization of surgically resectable cancers with a multi-analyte blood test. Science 2018;359:926–30. https://doi.org/10.1126/SCIENCE.AAR3247.
This file contains the following vectors:

*healthy cohort*

- cfDNA_healthy: all healthy individuals with detectable cmDNA.
- cfDNA_healthy_no_cmDNA: healthy individuals with no detectable cmDNA.

*age‑stratified healthy groups*

- cfDNA_healthy_leq_40y: ≤ 40 years
- cfDNA_healthy_40y_to_50y: 40–50 years
- cfDNA_healthy_over_50y: ≥ 50 years

*all‑stage cancer cohorts (stage I–III mixed)*

- cfDNA_lung
- cfDNA_breast_I: stage I breast cancer
- cfDNA_breast_II_III: stage II + III breast cancer
- cfDNA_colorectal
- cfDNA_pancreatic
- cfDNA_ovarian
- cfDNA_esophageal
- cfDNA_stomach
- cfDNA_liver

*stage I only cohorts for the cancers where enough stage‑I cases exist*

- cfDNA_lung_I
- cfDNA_colorectal_I
- cfDNA_ovarian_I
- cfDNA_stomach_I
- cfDNA_esophageal_I
- cfDNA_liver_I

cmDNA_Cohen_2018.mat 
Data on cmDNA concentrations (fragments/mL) and MAFs (%) for healthy individuals and cancer patients from 
Cohen JD, Li L, Wang Y, et al. Detection and localization of surgically resectable cancers with a multi-analyte blood test. Science 2018;359:926–30. https://doi.org/10.1126/SCIENCE.AAR3247.
This file contains the following vectors:

*healthy cohort* 

- cmDNA_healthy: cmDNA concentrations for all healthy individuals with detectable cmDNA.
- cmDNA_healthy_MAF_001, cmDNA_healthy_MAF_0015, cmDNA_healthy_MAF_002: same data filtered to include only patients with plasma MAF > 0.01%, 0.015%, or 0.02%.
- MAF healthy: plasma MAFs for all healthy individuals with detectable cmDNA.

*For each cancer type (lung, breast, colorectal, pancreatic, ovarian, esophageal, stomach, liver):*

- cmDNA_<type>_0: cmDNA concentrations for patients whose top plasma mutation was not detected in tumor.
- cmDNA_<type>_1: cmDNA concentrations for patients whose top plasma mutation was also detected in the tumor.
- cmDNA_<type>_0_MAF_001, cmDNA_<type>_0_MAF_0015, cmDNA_<type>_0_MAF_002
- cmDNA_<type>_1_MAF_001, cmDNA_<type>_1_MAF_0015, cmDNA_<type>_1_MAF_002: same data filtered to include only patients with plasma MAF > 0.01%, 0.015%, or 0.02%.
- MAF_<type>_0: plasma MAFs for patients whose top plasma mutation was not detected in tumor.
- MAF_<type>_1: plasma MAFs for patients whose top plasma mutation was also detected in the tumor.
- plasma_MAF_for_tumor_MAF_geq_30: plasma MAF values for all cancer patients whose tumor contained at least one mutation with tumor MAF ≥ 30%.

cfDNA_Cohen_2017.mat

Data on cfDNA concentrations (ng/mL) for healthy individuals and stage I-III pancreatic cancer patients from 
Cohen JD, Javed AA, Thoburn C, et al. Combined circulating tumor DNA and protein biomarker-based liquid biopsy for the earlier detection of pancreatic cancers. Proceedings of the National Academy of Sciences 2017;114:10202–7. https://doi.org/10.1073/PNAS.1704961114.
This file contains the following vectors:

- cfDNA_healthy: cfDNA concentrations for healthy individuals.
- cfDNA_pancreas: cfDNA concentrations for pancreatic cancer patients.

cfDNA_Mattox_2023.mat

Data on the concentrations (ng/mL) for total and leukocyte-shed cfDNA for healthy individuals and stage I-III pancreatic and ovarian cancer patients from
Mattox AK, Douville C, Wang Y, et al. The Origin of Highly Elevated Cell-Free DNA in Healthy Individuals and Patients with Pancreatic, Colorectal, Lung, or Ovarian Cancer. Cancer Discov 2023;13:2166–79. https://doi.org/10.1158/2159-8290.CD-21-1252.
This file contains the following matrices, each one having 5 columns:

- healthy, pancreatic, ovarian
- healthy_no_outliers, ovarian_no_outliers: same data excluding outliers (see the ESM of the paper).

The first column of each of the above matrices is total cfDNA concentrations;
second column is the leukocyte cfDNA as calculated by the Sun et al. QP deconvolution method;
third column is the leukocyte cfDNA as calculated by the Sun et al. NNLS deconvolution method;
fourth column is the leukocyte cfDNA as calculated by the Moss et al. QP deconvolution method;
fifth column is the leukocyte cfDNA as calculated by the Moss et al. NNLS deconvolution method.

== Code files (.m) ==

MATLAB codes that generate the results and figures of the paper.

1) healthy_vs_cancer_cfDNA_distribution.m: returns Fig. 1 of the paper, the distributions of cfDNA concentration for healthy individuals and early-stage cancer patients.

2) multiplicative_factor_search.m: This script scans a grid of multiplicative factors α and performs a two‑sample Kolmogorov–Smirnov test comparing α × healthy vs a chosen cancer cohort.
For each α it records the KS statistic, p‑value, and reject/non‑reject flag. It then identifies
- the α that minimizes the KS distance (best match between distributions), and
- the full 5% KS non‑rejection interval (all α values where the test does not reject).
This is the procedure used in the manuscript to determine the multiplicative increase α between healthy and cancer cfDNA/cmDNA distributions.

3) cfDNA_CDF_plots.m: returns Fig. 2 of the paper, the multiplicative shift between healthy and cancer CDFs for cfDNA concentrations.

4) non_tumor_cmDNA_CDF_plots.m: returns Fig. 3 of the paper the multiplicative shift between healthy and cancer CDFs for cmDNA concentrations, considering only cancer patients with top plasma mutation not detected in the tumor.

5) ROC_curves.m: returns the ROC curves and performance metrics for the cfDNA-driven cancer detection across cancer types, shown in the ESM of the paper

6) saturation_model_fitting.m:  This script applies the saturation‑clearance cfDNA model using parameters calibrated from healthy individuals from Mattox et al. (2023).
For the ovarian and pancreatic cohorts in Mattox (2023), it applies the cancer‑specific multiplicative increases in leukocyte and non‑leukocyte shedding,
generates model‑predicted cfDNA distributions, and compares them with data using:
- Kolmogorov–Smirnov tests (requiring p > 0.05 for non‑rejection), and
- agreement of the leukocyte vs non‑leukocyte correlation between model and data.
It also returns the panels for Fig. 4 of the paper.



