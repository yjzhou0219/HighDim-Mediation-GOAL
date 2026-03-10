# High-Dimensional Causal Mediation Analysis with GOAL

This repository contains the R code and simulation scripts for the paper **"Variable selection-combined causal mediation analysis for continuous treatments with application to high-dimensional biomedical data" (Yajing Zhou, et al.)**.

We implement a novel method that integrates Generalized Outcome-Adaptive Lasso (GOAL) with an IPW-based semi-parametric mediation estimator. This approach is specifically designed to estimate causal mediation effects for continuous exposures and mediators in the presence of high-dimensional covariates.


## Prerequisites

The source code was developed and successfully executed under R version 4.4.2 . To ensure full reproducibility of the results, we provide two methods for environment setup.


### Method 1: Manual Setup (Recommended for R 4.4.2 or earlier)

If you prefer to set up the environment manually, please ensure you are using R version 4.4.2 (2024-10-31) or an earlier version to maintain compatibility with all dependencies.


```r
# 1. Install 'devtools' to load packages from GitHub
if (!require("devtools")) install.packages("devtools")

# 2. Install standard CRAN packages
packages <- c("stats", "MASS", "glmnet", "car", "dplyr") # Add any other CRAN packages you use
install.packages(packages[!packages %in% installed.packages()])

# 3. Install specific packages from GitHub
library(devtools)
# Please ensure you have internet access for this step
install_github("strengejacke/sjstats")
install_github("JeffreyRacine/R-Package-np")

# 4. Install 'lqa' from the provided local source file
install.packages("renv/local/lqa_1.0-3.tar.gz", repos = NULL, type = "source")

```

### Method 2: Using **`renv`** (Automated Environment Recovery)
This repository includes a `renv` configuration to help you recreate the exact development environment.

```r
# 1. Prepare the Project
# set the working directory to the project folder 'R' or Clone the repository and open the project in RStudio
setwd("path/to/R/")

# 2. Initialize 'renv'. Ensure the renv package is installed on your system.
if (!require("renv")) install.packages("renv")


# 3. Restore the Environment
renv::restore()
```
- R should automatically source `renv/activate.R `upon opening the project. If the environment is not automatically activated, manually run the following step:


	```r
	## 2.2 Ativate the environment.
	source("renv/activate.R")
	```

- After running following `renv::restore()` command, all dependencies (including the local `lqa` source package from the `R/renv/local` folder) will be automatically installed.
- When prompted, confirm the installation. `renv` will synchronize your local library with the `renv.lock` file provided in this repository



#### **Note on `lqa` package**: 
The `lqa` package is required for GOAL but may **not be** directly available via `install.packages("lqa")` in recent versions of R. We have included the source file `lqa_1.0-3.tar.gz` in `renv/local/` of this repository. 

- If you are using Method 1 (Manual Setup), the `lqa` can be directly installed by setp 4:
	```r
	install.packages("renv/local/lqa_1.0-3.tar.gz", repos = NULL, type = "source")
	```

- If you are using Method 2 (`renv`): No additional action is needed. Run `renv::restore() ` and  `lqa` will be automatically installed from the provided local source file.




- If you are using Method 1 (Manual Setup), the `lqa` can be directly installed by setp 4:
	```r
	install.packages("renv/local/lqa_1.0-3.tar.gz", repos = NULL, type = "source")
	```

- If you are using Method 2 (`renv`): No additional action is needed. Run `renv::restore() ` and  `lqa` will be automatically installed from the provided local source file.






## Repository Structure 

All source code for the proposed method and simulation studies is located in the **`R/`** directory:


* **`GOAL_deviance.R`**: Contains the implementation of the Generalized Outcome-Adaptive Lasso (GOAL) algorithm for high-dimensional variable selection.
* **`mediation_effect_continuous.R`**: Contains the core function for the IPW-based semi-parametric estimator, specifically designed for continuous exposures and mediators.
* **`simulation.R`**: The main script to reproduce the simulation studies. It demonstrates how to generate data and apply the proposed estimator by calling the necessary functions from the aforementioned files.
* **`real-data analysis.R`**: The simple simplified demonstration script that outlines the data processing and analysis for the real-data application in this manuscript.
* **`plotting_Figs1-5.R`**: The plotting script for the core figures (Figures 1-5) presented in this manuscript.

You can copy and run the following R script to install all dependencies (including CRAN packages and specific GitHub repositories):



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




