# DiSTect

## High-dimensional Bayesian Model for Disease-Specific Gene Detection in Spatial Transcriptomics
<p align="center">
  <img width="300" height="300" alt="DiSTect_logo"
       src="https://github.com/user-attachments/assets/58053101-387b-40de-8967-2408396d6d69" />
</p>


Identifying disease-indicative genes is critical for deciphering disease mechanisms and continues to attract significant interest. Spatial transcriptomics offers unprecedented insights for
the detection of disease-specific genes by enabling within-tissue contrasts. However, this new
technology poses challenges for conventional statistical models developed for RNA-seq, as
these models often neglect the spatial organization of tissue spots. In this article, we propose
a Bayesian shrinkage model to characterize the relationship between high-dimensional gene
expressions and the disease status of tissue spots, incorporating spatial correlation among
these spots through autoregressive terms. Our model adopts a hierarchical structure to
accommodate for the missing data within tissues and is further extended to facilitate the
analysis of multiple correlated samples. To ensure the model’s applicability to datasets of
varying sizes, we carry out two computational frameworks for Bayesian parameter estimation, tailored to both small and large sample scenarios. 

## Installation

The **DiSTect** package depends on two core packages: **rstan** and **rjags**.  
Please ensure both dependencies are installed before installing DiSTect.

### Step 1. Install Dependencies

#### For all systems
To install **rstan**, run the following command in R:

```r
install.packages("rstan", dependencies = TRUE)
```

#### For macOS users
First, install JAGS via Homebrew in your terminal:
```bash
brew install jags
```

Then install rjags from source in R:
```r
install.packages("rjags", type = "source")
```

#### For Windows users
Please visit the JAGS official website
 to download and install JAGS manually.
After that, install rjags in R:
```r
install.packages("rjags")
```

### Step 2. Install DiSTect
Once all dependencies are installed, use the following commands to install DiSTect from GitHub:
```r
# install devtools if not already installed
install.packages("devtools")

# install the DiSTect package
devtools::install_github("StaGill/DiSTect")

# load the package
library(DiSTect)
```


## Tutorial

A R Markdown of the tutorial is accessible from: [https://qihuangzhang.github.io/software/DiSTect_tutorial](https://qihuangzhang.github.io/software/DiSTect_tutorial).

## Data

The dataset about HER2-positive breast is accessible from https://github.com/almaan/her2st.

## How to cite DiSTect

Use the following BibTeX entry:

```bibtex
@article{zhao2025distect,
  title={DiSTect: a Bayesian model for disease-associated gene discovery and prediction in spatial transcriptomics},
  author={Zhao, Qicheng and Deng, Anji and Zhang, Qihuang},
  journal={Bioinformatics},
  pages={btaf530},
  year={2025},
  publisher={Oxford University Press}
}

