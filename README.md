# Pòlya-Gamma Extended Transition Diagnostic Classification Model

[![DOI](https://zenodo.org/badge/DOI/10.1017/psy.2025.10031.svg)](https://doi.org/10.1017/psy.2025.10031)

This repository provides R code implementing an extension of the Transition Diagnostic Classification Model (TDCM; Matthew & Bradshaw, 2018) using Pòlya-Gamma data augmentation. The method enables flexible longitudinal modeling of student attribute mastery with covariate-dependent transitions, as described in:

**Resch, J., Baugh, S., Duan, H., Tang, J., Madison, M., Cotterell, M., & Jeon, M. (2025).**  
*Bayesian Transition Diagnostic Classification Models with Pòlya-Gamma Augmentation*. 
[https://doi.org/10.1017/psy.2025.10031](https://doi.org/10.1017/psy.2025.10031)

## Features

- Supports longitudinal diagnostic classification modeling
- Implements log-linear cognitive diagnosis model (LCDM) for the item response component
- Incorporates covariates into attribute transition modeling between time points
- Employs Pòlya-Gamma augmentation for efficient Gibbs sampling
- Demonstrated via simulation and empirical examples

## Working with Respository

To replicate the results in the manuscript or run the extended TDCM model on new data, see the main R script in the repository and refer to simulated/empirical data examples.
