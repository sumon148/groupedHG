# groupedHG

**Grouped Hypergeometric and Binomial Sampling with Sensitivity and Specificity**

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

`groupedHG` is an R package implementing probability mass functions, expectations, variances, and information measures for group sampling processes where detection tests may be imperfect. The package was developed to support research and applications in epidemiology, biosecurity, and quality control, as described in:

> **Barnes, B., Parsa, M., Das, S., & Clark, R. (2025). Hypergeometric and binomial group sampling with sensitivity and specificity.**
> *Communications in Statistics*, 

The methodology is described in detail in the accompanying paper and supplementary material.

> 📘 **Repository:**
> [https://github.com/sumon148/groupedHG](https://github.com/sumon148/groupedHG)

## Overview

Group (pooled) sampling is widely used to detect contamination, infection, or defects while reducing cost and effort compared to individual testing. However, standard models often assume perfect test accuracy. In practice, sensitivity and specificity can vary with contamination levels or group sizes (e.g., dilution effects).

This package provides:

* Analytical distributions for grouped hypergeometric and binomial sampling, incorporating imperfect sensitivity and specificity.
* Support for **group-level** or **item-level** detection functions, including PCR and serological test models.
* Functions to compute expectations, variances, and Hellinger (Fisher) information to guide sampling design.

## Installation

The package is **not** available on CRAN. You can install it directly from GitHub.

### Install without vignettes

This is faster if you do not need the tutorial vignette:

```r
# install.packages("devtools") if needed
devtools::install_github("sumon148/groupedHG")
```

### Install with vignette (recommended)

To install **and build the vignette**, use:

```r
devtools::install_github("sumon148/groupedHG", build_vignettes = TRUE)
```

> **Note:** Building vignettes requires `rmarkdown` and `knitr`:
>
> ```r
> install.packages(c("rmarkdown", "knitr"))
> ```

---

## Getting Started

After installation, load the package:

```r
library(groupedHG)
```

To see the vignette tutorial:

1. **List available vignettes:**

```r
vignette(package = "groupedHG")
```

2. **Open the vignette:**

```r
vignette("groupedHG")
```

The vignette includes:

* Computing probability mass functions under different sensitivity/specificity scenarios
* Item-level versus group-level detection
* Variance and covariance estimation
* Hellinger information to evaluate sampling designs
* Example plots and interpretations

---

## Example (Quick Preview)

Below is a minimal illustration. For detailed use cases, refer to the vignette.

```r
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

* Barnes, B., Parsa, M., Das, S., & Clark, R. (2025). Hypergeometric and binomial group sampling with sensitivity and specificity. *Communications in Statistics – Theory and Methods* (under review).



## Citation

If you use this package in publications, please cite:

* Barnes, B., Parsa, M., Das, S., & Clark, R. (2025). Hypergeometric and binomial group sampling with sensitivity and specificity. *Communications in Statistics – Theory and Methods* (under review).


## License

This package is released under the MIT License.

---

✅ **Tip:**
If you want a **shorter installation example**, you can remove the "Install without vignettes" section, but many users appreciate the clarity.

---

**Let me know if you’d like help adding badges, more examples, or tweaking formatting!**
