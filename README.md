# groupedHG

**Grouped Hypergeometric and Binomial Sampling with Sensitivity and Specificity**

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

`groupedHG` is an R package implementing probability mass functions, expectations, variances, and information measures for group sampling processes where detection tests may be imperfect. The package was developed to support research and applications in epidemiology, biosecurity, and quality control, as described in:

> **Barnes, B., Parsa, M., Das, S., & Clark, R. (2024). Hypergeometric and binomial group sampling with sensitivity and specificity.**  
> *Journal of Applied Statistics*, [DOI/link if available]

The methodology is described in detail in the accompanying paper and supplementary material.

> 📘 **Repository:**  
> [https://github.com/sumon148/groupedHG](https://github.com/sumon148/groupedHG)

## Overview

Group (pooled) sampling is widely used to detect contamination, infection, or defects while reducing cost and effort compared to individual testing. However, standard models often assume perfect test accuracy. In practice, sensitivity and specificity can vary with contamination levels or group sizes (e.g., dilution effects).

This package provides:

- Analytical distributions for grouped hypergeometric and binomial sampling, incorporating imperfect sensitivity and specificity.
- Support for **group-level** or **item-level** detection functions, including PCR and serological test models.
- Functions to compute expectations, variances, and Hellinger (Fisher) information to guide sampling design.

## Installation

The package is **not** available on CRAN. Install directly from GitHub using `devtools`:

```r
# install.packages("devtools") if needed
devtools::install_github("sumon148/groupedHG")
```

## Getting Started

For a step-by-step tutorial—including examples, plots, and code to reproduce the main results—please see the **vignette**:

```r
# Load the package
library(groupedHG)

# View the vignette
vignette("groupedHG")
```

The vignette covers:

* Computing probability mass functions under different sensitivity/specificity scenarios.
* Item-level versus group-level detection.
* Variance and covariance estimation.
* Hellinger information to evaluate sampling designs.
* Example plots and interpretations.

## Example (Quick Preview)

Below is a minimal illustration. For detailed use cases, refer to the vignette.

```{r example, eval=TRUE}
# Load library
library(groupedHG)

# Example parameters
N <- 100
TX <- 5
b <- 10
n_bar <- 5
Delta <- 0.95
Lambda <- 0.99

# Compute probability mass function
pmf <- groupedHG_pmf(N, TX, b, n_bar, Delta, Lambda)

pmf
```
## Methods

The package implements the analytical results described in:

* Barnes et al. (2025). *Hypergeometric and binomial group sampling with sensitivity and specificity*. 

## Citation

If you use this package in publications, please cite:

Barnes, B., Parsa, M., Das, S., & Clark, R. (2025). *Hypergeometric and binomial group sampling with sensitivity and specificity*. 

## License

This package is released under the MIT License.

