---

## title: "groupedHG: Figure Gallery" author: "Your Name" date: "`r Sys.Date()`" output: rmarkdown::html\_vignette: toc: true toc\_depth: 2 vignette: > %\VignetteIndexEntry{groupedHG figure gallery} %\VignetteEngine{knitr::rmarkdown} %\VignetteEncoding{UTF-8}

```{r
library(knitr)
opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 6, fig.height = 4,
  dev = "png"
)
library(groupedHG)
```

## 1. Figure 1: Hypergeometric PMF (Perfect Test)

**Paper reference:** "Figure 1. Distribution of the number of positive pools under exact hypergeometric sampling."

```r
N    <- 100
tx   <- 20
b    <- 4
nbar <- N / b

pmf <- pmf_hg(tx = tx, N = N, b = b, nbar = nbar)

plot(0:b, pmf, type = "h", lwd = 2,
     xlab = "Number of positive pools",
     ylab = "P(Y = y)",
     main = "Exact Hypergeometric PMF")
```

> **Caption:** As contamination increases, the mass shifts right; note the discrete support 0…b.

## 2. Figure 2: Exact vs. Binomial Approximation

**Paper reference:** "Figure 2. Comparison of exact HG vs. binomial-approximate PMFs."

```r
pmf_exact <- pmf_hg(tx = tx, N = N, b = b, nbar = nbar)
pmf_bin   <- pmf_bn_group(tx = tx, N = N, b = b, nbar = nbar,
                          Delta = 1, Lambda = 1)

plot(0:b, pmf_exact, type = "h", lwd = 2,
     xlab = "Positive pools", ylab = "Probability",
     main = "Exact HG vs. Binomial Approximation")
points(0:b, pmf_bin, pch = 19, col = "red")
legend("topright", legend = c("Exact","Binomial"),
       lty = c(1, NA), pch = c(NA, 19), col = c("black","red"))
```

> **Caption:** Even with perfect tests, binomial can under/overestimate tails.

## 3. Figure 3: Hellinger Information Curve

**Paper reference:** "Figure 3. Hellinger information as a function of true contamination level."

```r
tx_vals <- seq(0, 50, by = 1)
info    <- sapply(tx_vals, function(x)
  hellinger_info(tx = x, N = N, b = b, nbar = nbar,
                 method = "group", Delta = 0.7, Lambda = 0.8))

plot(tx_vals, info, type = "l", lwd = 2,
     xlab = expression(T[X]),
     ylab = "Hellinger Information",
     main = "Design Information vs. Contamination")
```

> **Caption:** Peaks where sampling is most informative; shifts with Δ/Λ.

## 4. Figure 4: Sensitivity Curve (Salmonella Model)

**Paper reference:** "Figure 4. Detection probability vs. number contaminated in pool (item-level model)."

```r
sensitivity_salmonella <- function(tx, ws = 20, C = 1, rho = 0.1, nbar) {
  1 - exp(-ws * C * (1 - exp(-rho / (ws * nbar))) * tx / nbar)
}

tx_seq <- seq(0, 60, by = 1)
sens   <- sensitivity_salmonella(tx = tx_seq, nbar = 5)

plot(tx_seq, sens, type = "l", lwd = 2,
     xlab = expression(T[X] / bar(n)),
     ylab = "P(detect)",
     main = "Item-Level Sensitivity vs. Pool Contamination")
```

> **Caption:** Models dilution: more positives → higher detect probability, asymptoting at 1.

---

### Usage

1. Copy to **vignettes/groupedHG-figures.Rmd**.
2. Run `devtools::build_vignettes()` to compile.

