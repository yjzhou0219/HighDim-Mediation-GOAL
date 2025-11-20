# High-Dimensional Causal Mediation Analysis with GOAL

This repository contains the R code and simulation scripts for the paper **"Variable selection-combined causal mediation analysis for continuous treatments with application to high-dimensional biomedical data"**.

We implement a novel method that integrates Generalized Outcome-Adaptive Lasso (GOAL) with an IPW-based semi-parametric mediation estimator. This approach is specifically designed to estimate causal mediation effects for continuous exposures and mediators in the presence of high-dimensional covariates.


## Prerequisites

The code is written in R. To reproduce the analysis, you need to install the required R packages.

You can copy and run the following R script to install all dependencies (including CRAN packages and specific GitHub repositories):

```r
# 1. Install 'devtools' to load packages from GitHub
if (!require("devtools")) install.packages("devtools")

# 2. Install standard CRAN packages
packages <- c("stats", "MASS", "glmnet", "car", "lqa", "dplyr") # Add any other CRAN packages you use
install.packages(packages[!packages %in% installed.packages()])

# 3. Install specific packages from GitHub
library(devtools)
# Please ensure you have internet access for this step
install_github("strengejacke/sjstats")
install_github("JeffreyRacine/R-Package-np")
```


## Repository Structure

The repository is organized as follows:
* **`R/`**: This folder contains the main simulation script and the source code for the proposed method.
    * `GOAL_deviance.R`: Contains the function for the deviance-determined generalized outcome-adaptive Lasso (GOAL) variable selection algorithm.
    * `mediation_effect_continuous.R`: Contains the three functions for the IPW-based semi-parametric estimator for continuous exposures and mediators.
    * `simulation.R`: The main script to reproduce the simulation studies. It calls the necessary functions from the aforementioned files.
* **`LICENSE`**: The license agreement for using the code.
* **`README.md`**: Project documentation and instructions.


## How to Run


To replicate the simulation results, please follow these steps:

1.  **Clone or download** this repository to your local machine.
2.  Open the `R/simulation.R` file in R or RStudio.
3.  **Set the working directory** to the `R/` folder (where the script is located). This ensures that the source files are loaded correctly.
4.  **Run the script**. The `simulation.R` file is self-contained (assuming packages are installed). It will automatically source the helper functions (`GOAL_deviance.R` and `mediation_effect_continuous.R`) located in the same directory.

    ```r
    # Run the main simulation
    source("simulation.R")
    ```


