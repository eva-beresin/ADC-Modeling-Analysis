# ADC-Modeling-Analysis
Dimensional Analysis and Scale Modeling of Antibody-Drug Conjugate Delivery Systems

``` markdown

## Author
Eva Beresin  
Spring 2026  

---

## Overview
This MATLAB project performs dimensional analysis, non-dimensional pi group calculations, model fitting (multi-linear, exponential, and power fits), and evaluation of critical physical parameters for experimental observation data.

---

## Project Structure & Workflow

### **Part 1: Dimensional Analysis & Pi Group Calculation**
1. **Data Generation:** Loads observation data using `ProjectDataGenerator(175210)`.
2. **Dimensional Matrix Formulation:** 
   - Defines physical fundamental dimensions (Mass $M$, Length $L$, Time $T$) for basis and non-basis parameters.
   - Applies Reduced Row Echelon Form (`rref`) to derive the Pi Exponent Table.
3. **Pi Group Construction:**
   - Computes 6 distinct $\pi$ groups ($\pi_1$ through $\pi_6$) using derived exponents and measurement data.
   - Outputs formatted scientific notation tables of all computed values.

---

### **Part 2: Model Fitting & Scaling**
1. **Data Normalization:** Scales $\pi$ group data relative to maximum observed values.
2. **Variable Designation:** Sets $\pi_1$ (Length 2) as the dependent target variable ($y$) and $\pi_2 \dots \pi_6$ as independent predictor variables ($x$).
3. **Regression Fitting:**
   - **Multi-Linear Fit:** Executed using both `fitlm` and nonlinear optimization `fitnlm`.
   - **Exponential Fit:** Fits $y = b_1 \cdot \exp(b_2 x_1 + b_3 x_2 + \dots)$ using log-transformed initial parameter estimates.
   - **Power Law Fit:** Fits $y = b_1 \cdot x_1^{b_2} \cdot x_2^{b_3} \dots$ using multi-variable log-log regression initial guesses.

---

### **Part 3: Model Evaluation & Critical Pi Group Identification**
1. **Coefficient Unscaling:** Translates scaled model parameters back to original physical scales.
2. **Model Selection:** Evaluates Akaike Information Criterion (AIC) values across Multi-Linear, Exponential, and Power models to identify the optimal predictor.
3. **Critical Parameter Extraction:**
   - Applies Bonferroni correction for statistical significance ($\alpha = 0.01 / 5 = 0.002$).
   - Filters statistically significant non-basis parameters based on $p$-values from the best-performing fit model.

---

## Requirements & Environment
- **MATLAB Version:** MATLAB R2024b or newer
- **Required Toolboxes:** Statistics and Machine Learning Toolbox (`fitlm`, `fitnlm`)
- **Required Files:** `ProjectDataGenerator.m` (must be present in the working path)

---

## Instructions to Run
1. Place the code into a `.m` file (e.g., `BMED3310_Computational_Project.m`).
2. Ensure `ProjectDataGenerator.m` is located in your MATLAB path/working directory.
3. Open MATLAB and execute the file script directly or press **Run** (F5).

```
