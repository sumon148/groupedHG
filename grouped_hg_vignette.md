---

## title: "groupedHG: Hypergeometric & Binomial Group Sampling" author: "Your Name" date: "`r Sys.Date()`" output: rmarkdown::html\_vignette: toc: true toc\_depth: 3 vignette: > %\VignetteIndexEntry{groupedHG overview and examples} %\VignetteEngine{knitr::rmarkdown} %\VignetteEncoding{UTF-8}

## 1. Introduction

- **Problem statement:** pooled (group) sampling when tests have imperfect sensitivity/specificity
- **Applications:** biosecurity, public health, ecology
- **Key features:** exact hypergeometric formulas, binomial approximations, Hellinger information for design

## 2. Installation

```r
# From CRAN (when available)
install.packages("groupedHG")

# Or development version from GitHub
remotes::install_github("sumon148/groupedHG")
```

## 3. Mathematical Background

- **Notation:**\
  \(N\): population size\
  \(T_X\): number of contaminated items\
  \(b\), \(\barN\): number and size of pools
- **Core formulas:**
  - Hypergeometric PMF (perfect test)
  - Group-level \(\Delta/\Lambda\) and item-level \(\delta/\lambda\) adjustments
  - Binomial approximations
  - Hellinger information

*(Include key equations from the paper here.)*

## 4. Basic Usage

### 4.1 Perfect-Test PMF

```r
library(groupedHG)
pmf <- pmf_hg(tx = 20, N = 100, b = 4, nbar = 4)
plot(0:4, pmf, type = "h", lwd = 2,
     xlab = "Number of positive pools",
     ylab = "P(Y = y)")
```

### 4.2 Imperfect Tests

```r
# Group-level sensitivity/specificity
pmf_gl <- pmf_hg_group(tx = 20, N = 100, b = 4, nbar = 4,
                       Delta = 0.7, Lambda = 0.8)

# Item-level sensitivity/specificity
pmf_il <- pmf_hg_item(tx = 20, N = 100, b = 4, nbar = 4,
                      delta = 0.7, lambda = 0.8)
```

## 5. Binomial Approximations

```r
pmf_gl_bn <- pmf_bn_group(tx = 20, N = 100, b = 4, nbar = 4,
                          Delta = 0.7, Lambda = 0.8)
pmf_il_bn <- pmf_bn_item(tx = 20, N = 100, b = 4, nbar = 4,
                         delta = 0.7, lambda = 0.8)
```

*(Compare exact vs. approximate curves.)*

## 6. Hellinger Information

```r
tx_vals <- seq(0, 50, by = 1)
info <- sapply(tx_vals, function(x) {
  hellinger_info(tx = x, N = 100, b = 4, nbar = 4,
                 method = "group", Delta = 0.7, Lambda = 0.8)
})
plot(tx_vals, info, type = "l",
     xlab = "T_X", ylab = "Hellinger information")
```

## 7. Case Study: Salmonella in Pig Pens

1. **Setup:** 60 pigs, pooled fecal samples
2. **Sensitivity model:** Eq. 21 from the paper
3. **Design exploration:** grid search over \(b\) and \(\bar n\)

```r
design_grid <- expand.grid(b = 1:20, nbar = c(1, 2, 4, 5, 10, 20))
design_grid$info <- mapply(function(b, nbar) {
  hellinger_info(
    tx = 0.15 * 60, N = 60, b = b, nbar = nbar,
    method = "item", delta_fun = sensitivity_salmonella,
    lambda = 1
  )
}, design_grid$b, design_grid$nbar)
# Visualize with a heatmap or contour plot
```


## 8. Session Info & References

```r
sessionInfo()
```

> **References**\
> – Barnes *et al.* (2025), *Hypergeometric & Binomial Group Sampling*\
> – Arnold *et al.* (2005), *Salmonella in pigs*\
> – Theobald & Davie (2014), *Pooled hypergeometric designs*

