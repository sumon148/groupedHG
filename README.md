
<!-- README.md is generated from README.Rmd. Please edit that file -->

# groupedHG

<!-- badges: start -->
<!-- badges: end -->

The goal of **groupedHG** is to estimate probability mass function (PMF)
for grouped hypergeometric distribution considering adjustment for
sensitivity and specificity of the test. In addition to this, Fisher
Information of PMF with respect to the total number of contaminated
items assumed in the population (denoted as $T_x$) can be calculated
using two alternative methods - analytic derivative (AD) and PMF based.
The PMF based FI calculation is based on the Sánchez-Moreno, Yánez and
Dehesa (2009) paper on discrete densities and Fisher information.

Sánchez-Moreno, P., Yánez, R. J., & Dehesa, J. S. (2009, October).
Discrete densities and Fisher information. In Proceedings of the 14th
International Conference on Difference Equations and Applications.
Difference Equations and Applications. Istanbul, Turkey: Bahçesehir
University Press (pp. 291-298).

## Installation

You can install the development version of mypackage like so:

``` r
install.packages("devtools")
devtools::install_github("sumon148/groupedHG")
```

## Example

Let assume a population containing $N$ items of which $Tx$ are
contaminated. We can calculate PMF of having $ty$ contaminated groups of
size $barN$ after inspecting a sample of $b$ groups as below assuming
group-level imperfect sensitivity ($\Delta=0.7$) and specificity
($\Lambda=0.8$).

``` r
library(groupedHG)
b=4
ty.values <- c(0:b)
PRty.HG.group <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N=100, barN=4, Tx=20, b=b, delta=0.7, lambda=0.8))
plot(0:b, PRty.HG.group, type = "l", col = "black", xlab=expression(paste("Number of group detections, ", t[y])),
     ylab = "Probability mass", main = "Grouped Hypergeometric model",lty=2)
points(0:b,PRty.HG.group,type = "b", pch = 19, col = "black")
```

<img src="man/figures/README-PMF HG Group-1.png" width="100%" />

In the similar way, the PMF can be calculated assuming item-level
imperfect test as below:

``` r
PRty.HG.item <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N=100, barN=4, Tx=20, b=b, delta=0.7, lambda=0.8))
plot(0:b, PRty.HG.item, type = "l", col = "black", xlab=expression(paste("Number of group detections, ", t[y])),
     ylab = "Probability mass", main = "Grouped Hypergeometric model",lty=2)
points(0:b,PRty.HG.item,type = "b", pch = 19, col = "black")
```

<img src="man/figures/README-PMF HG Item-1.png" width="100%" />

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
