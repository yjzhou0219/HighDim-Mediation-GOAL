# High-Dimensional Causal Mediation Analysis with GOAL

This repository contains the R code and simulation scripts for the paper **"Variable selection-combined causal mediation analysis for continuous treatments with application to high-dimensional biomedical data" (Yajing Zhou, et al.)**.

We implement a novel method that integrates Generalized Outcome-Adaptive Lasso (GOAL) with an IPW-based semi-parametric mediation estimator. This approach is specifically designed to estimate causal mediation effects for continuous exposures and mediators in the presence of high-dimensional covariates.


## Prerequisites

The source code was developed and successfully executed under R version 4.4.2 (2024-10-31). To ensure full reproducibility of the results, we provide two methods for environment setup.


### Method 1: Manual Setup (Recommended for R 4.4.2 or earlier)

If you prefer to set up the environment manually, please ensure you are using R version 4.4.2 (2024-10-31) or an earlier version to maintain compatibility with all dependencies.


```r
# Step 1: Install 'devtools' to load packages from GitHub
if (!require("devtools")) install.packages("devtools")

# Step 2: Install standard CRAN packages
packages <- c("stats", "MASS", "glmnet", "car", "dplyr") # Add any other CRAN packages you use
install.packages(packages[!packages %in% installed.packages()])

# Step 3: Install specific packages from GitHub
library(devtools)
# Please ensure you have internet access for this step
install_github("strengejacke/sjstats")
install_github("JeffreyRacine/R-Package-np")

# Step 4: Install 'lqa' from the provided local source file
setwd("path/to/HighDim-Mediation-GOAL/R")
install.packages("renv/local/lqa_1.0-3.tar.gz", repos = NULL, type = "source")

```


### Method 2: Using **`renv`** (Automated Environment Recovery)
This repository includes a `renv` configuration to help you recreate the exact development environment.

```r
# Step 1: Prepare the Project
# Clone the repository, then set the working directory to the 'R' folder, or Clone the repository and open the project in RStudio
setwd("path/to/HighDim-Mediation-GOAL/R")

# Step 2: Initialize 'renv'. Ensure the renv package is installed on your system.
if (!require("renv")) install.packages("renv", repos = "https://cloud.r-project.org")

# Step 3: Restore the project library
renv::restore(project = getwd(), prompt = FALSE)
```
- R should automatically source `renv/activate.R `upon opening the project. If the environment is not automatically activated, manually run the following step:


	```r
	## Activate the environment
	source("renv/activate.R")
	```

- After running following `renv::restore()` command, all dependencies (including the local `lqa` source package from the `R/renv/local` folder) will be automatically installed.
- When prompted, confirm the installation. `renv` will synchronize your local library with the `renv.lock` file provided in this repository



#### **Note on `lqa` package**: 
The `lqa` package is required for GOAL but may **not be** directly available via `install.packages("lqa")` in recent versions of R. We have included the source file `lqa_1.0-3.tar.gz` in `R/renv/local/` of this repository. 

- If you are using **Method 1** (Manual Setup), the `lqa` can be directly installed by setp 4:
	```r
   # set the working directory
   setwd("path/to/HighDim-Mediation-GOAL/R")

   # install from the local source file
   install.packages("renv/local/lqa_1.0-3.tar.gz", repos = NULL, type = "source")

   # Check lqa installation
   library(lqa)
   packageVersion("lqa")
	```
	- In addition to downloading the package from the provided local source file, the `lqa` package can also be downloaded directly from the official website from CRAN Archive: https://cran.r-project.org/src/contrib/Archive/lqa/


- If you are using **Method 2** (`renv`): No additional action is needed. After following the steps in Method 2 and running `renv::restore() `, `lqa` package will be automatically installed from the provided local source file.




### Method 3: Using Docker (Recommended for Full Reproducibility)

To eliminate any operating system conflicts or dependency issues (such as difficulties installing the `lqa` package), we provide a Docker environment configured with `renv`. This ensures the simulation runs in an isolated container identically matching our exact computational environment. We have provided a `Dokerfile` in this repository alongwith the `R` file.

**Prerequisites**: Ensure [Docker](https://www.docker.com/products/docker-desktop) is installed and running.

**Step 1: Build the Docker Image** Open your Terminal or Command Prompt (CMD), navigate to the root directory of this repository (e.g., `cd path/to/HighDim-Mediation-GOAL/`), and run the following command. This will automatically set up R 4.4.2, install system libraries, and restore exact package versions via `renv.lock` (this process may take a few minutes):

```r
docker build -t goal-mediation-env .
```
**Step 2: Run the Docker Container** Launch an interactive container from the built image:
```r
docker run -it --rm goal-mediation-env
```


**Step 3: Execute the Simulation** You will now be inside the R console within the R/ directory. The `lqa` package and all other dependencies have been pre-installed. You can then run the provided script as runing the regular R code:
```r
source("simulation.R")
```





## Repository Structure 

All source code for the proposed method and simulation studies is located in the **`R/`** directory:


* **`GOAL_deviance.R`**: Contains the implementation of the Generalized Outcome-Adaptive Lasso (GOAL) algorithm for high-dimensional variable selection.
* **`mediation_effect_continuous.R`**: Contains the core function for the IPW-based semi-parametric estimator, specifically designed for continuous exposures and mediators.
* **`simulation.R`**: The main script to reproduce the simulation studies. It demonstrates how to generate data and apply the proposed estimator by calling the necessary functions from the aforementioned files.
* **`real-data analysis.R`**: The simple simplified demonstration script that outlines the data processing and analysis for the real-data application in this manuscript.
* **`plotting_Figs1-5.R`**: The plotting script for the core figures (Figures 1-5) presented in this manuscript.

You can copy and run the following R script to install all dependencies (including CRAN packages and specific GitHub repositories):





## How to Run (Non-Docker version)


To replicate the simulation results, please follow these steps:

1.  **Clone or download** this repository to your local machine.
2.  **Open** the `R/simulation.R` file in R or RStudio.
3.  **Set the working directory** to the `R/` folder (where the script is located). This ensures that the source files are loaded correctly.
4.  **Run the script**. The `simulation.R` file is self-contained (assuming packages are installed). It will automatically source the helper functions (`GOAL_deviance.R` and `mediation_effect_continuous.R`) located in the same directory.

    ```r
    # Run the main simulation
    source("simulation.R")
    ```


