
# groupedHG

The goal of **groupedHG** is to estimate probability mass function (PMF)
for grouped hypergeometric distribution considering adjustment for
sensitivity and specificity of the test. The sensitivity and specificity
can be considered at the item and group level. In addition to this,
Fisher Information (FI) of PMF with respect to the total number of
contaminated items assumed in the population (denoted as $T_x$) can be
calculated using two alternative methods - analytic derivative (AD) and
PMF-based. The PMF-based FI calculation is based on the Sánchez-Moreno,
Yánez and Dehesa (2009).

## Reference

Sánchez-Moreno, P., Yánez, R. J., & Dehesa, J. S. (2009, October).
Discrete densities and Fisher information. In Proceedings of the 14th
International Conference on Difference Equations and Applications.
Difference Equations and Applications. Istanbul, Turkey: Bahçesehir
University Press (pp. 291-298).

## Installation

You can install the development version of `groupedHG` like so:

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
PRty.HG.group <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N=100, barN=4, Tx=20, b=b, delta=0.7, lambda=0.8,verbose = FALSE))
plot(0:b, PRty.HG.group, type = "l", col = "black", xlab=expression(paste("Number of group detections, ", t[y])),
     ylab = "Probability mass", main = "Grouped Hypergeometric model",lty=2)
points(0:b,PRty.HG.group,type = "b", pch = 19, col = "black")
```

<img src="man/figures/README-PMF HG Group-1.png" width="100%" />

In the similar way, the PMF can be calculated assuming item-level
imperfect test as below:

``` r
PRty.HG.item <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N=100, barN=4, Tx=20, b=b, delta=0.7, lambda=0.8,verbose =FALSE))
PRty.HG.item
#> [1] 0.002993711 0.040346066 0.198811210 0.424987223 0.332861789
```

Expectation and variance of the grouped HG distribution for the number
of contaminated groups can be obtained as below:

``` r
ExpVarHG.imperfect.item <- ExpVarHG.imperfect(N=100, Tx=20, barN=4, b=4, delta=0.70, lambda=0.8,type="item")
ExpVarHG.imperfect.item$expectation
#> [1] 3.044377
```

``` r
ExpVarHG.imperfect.item$variance
#> [1] 0.7180313
```

``` r

ExpVarHG.perfect.group <- ExpVarHG.imperfect(N=100, Tx=20, barN=4, b=4, delta=1, lambda=1,type="group")
ExpVarHG.perfect.group$expectation
#> [1] 2.386647
```

``` r
ExpVarHG.perfect.group$variance
#> [1] 0.8797255
```

The Fisher Information of PMF with respect to given number of
contaminated items in the population $T_x$ assuming group and item level
imperfect sensitivity and specificity can be calculated as below:

``` r
Tx_values <- c(0:80)
# Different cases: Item level sensitivity
FI_case1_AD_Based_group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "group",verbose=FALSE))
FI_case1_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "item",verbose=FALSE))
plot(Tx_values, FI_case1_AD_Based_group, type = "l", col = "black", xlab=expression(paste("Number of Contaminated Items, ", T[x])),
     ylab = "Fisher Information: Hypergeometric", main = "FI under HG: Group and item level imperfect test",lty=2)
lines(Tx_values,FI_case1_AD_Based_Item,type = "l", lty = 2, col = "blue")
legend("topleft",legend=c("Group","Item"),col=c("black","blue"),lty=c(1,2))
```

<img src="man/figures/README-FI Item-1.png" width="100%" />
