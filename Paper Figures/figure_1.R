# Paper Figure 1 -----------------------
library(groupedHG)

N <- 100  # Total population size
Tx <- 20  # Number of contaminated items
barN <- 4  # Group size
b <- 4     # Number of groups
delta <- 0.7  # Test sensitivity
lambda <- 0.8  # Test Specificity


ty.values <- c(0:b)


# Based on Grouped-HG Model
# Different cases: Group level sensitivity
PRty.HG.group <- sapply(ty.values, function(ty) pmf_hg_group_imperfect(ty, N, barN, Tx, b, delta=1.0, lambda=1.0))
PRtyDelta.HG <- sapply(ty.values, function(ty) pmf_hg_group_imperfect(ty, N, barN, Tx, b, delta=delta, lambda=1.0))
PRtyDeltaLambda.HG <- sapply(ty.values, function(ty) pmf_hg_group_imperfect(ty, N, barN, Tx, b, delta=delta, lambda=lambda))

# Different cases: Item level sensitivity
PRty.HG.item <- sapply(ty.values, function(ty) pmf_hg_item_imperfect(ty, N, barN, Tx, b, delta=1.0, lambda=1.0))
PRty.HG.item.perfect <- sapply(ty.values, function(ty) pmf_hg_perfect(ty, N, barN, Tx, b))
PRtydelta.HG <- sapply(ty.values, function(ty) pmf_hg_item_imperfect(ty, N, barN, Tx, b, delta=delta, lambda=1.0))
PRtydeltalambda.HG <- sapply(ty.values, function(ty) pmf_hg_item_imperfect(ty, N, barN, Tx, b, delta=delta, lambda=lambda))

# Based on Grouped-BN Model
# Different cases: Group level sensitivity
PRty.BN.group <- sapply(ty.values, function(ty) pmf_bn_group_imperfect(ty, Tx, b, barN, N, delta=1.0, lambda=1.0))
PRtyDelta.BN <- sapply(ty.values, function(ty) pmf_bn_group_imperfect(ty, Tx, b, barN, N, delta=delta, lambda=1.0))
PRtyDeltaLambda.BN <- sapply(ty.values, function(ty) pmf_bn_group_imperfect(ty, Tx, b, barN, N, delta=delta, lambda=lambda))

# Different cases: Item level sensitivity
PRty.BN.item <- sapply(ty.values, function(ty) pmf_bn_item_imperfect(ty, Tx, b, barN, N, delta=1.0, lambda=1.0))
PRtydelta.BN <- sapply(ty.values, function(ty) pmf_bn_item_imperfect(ty, Tx, b, barN, N, delta=delta, lambda=1.0))
PRtydeltalambda.BN <- sapply(ty.values, function(ty) pmf_bn_item_imperfect(ty, Tx, b, barN, N, delta=delta, lambda=lambda))





# FI: Grouped-HG model for Tx : 0 - N
Tx_values <- c(0:100)
# Different cases: Item level sensitivity
FI_case1_PMF_Based_Item_HG <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF-HI", type = "item"))
FI_case2_PMF_Based_Item_HG <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF-HI", type = "item"))
FI_case4_PMF_Based_Item_HG <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF-HI", type = "item"))
# Different cases: Group level sensitivity
FI_case1_PMF_Based_Group_HG <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF-HI", type = "group"))
FI_case2_PMF_Based_Group_HG <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF-HI", type = "group"))
FI_case4_PMF_Based_Group_HG <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF-HI", type = "group"))


# FI: Grouped-BN model for Tx : 0 - N
# Different cases: Item level sensitivity
FI_case1_PMF_Based_Item_BN <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF-HI", type = "item"))
FI_case2_PMF_Based_Item_BN <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF-HI", type = "item"))
FI_case4_PMF_Based_Item_BN <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF-HI", type = "item"))
# Different cases: Group level sensitivity
FI_case1_PMF_Based_Group_BN <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF-HI", type = "group"))
FI_case2_PMF_Based_Group_BN <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF-HI", type = "group"))
FI_case4_PMF_Based_Group_BN <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF-HI", type = "group"))





png("Paper Figures/Figue 1.png",width = 10,height = 8,units = "in", res = 300)
par(mfrow=c(2,2))

range_pr_ty <- range(PRty.HG.group,PRtyDelta.HG,PRtyDeltaLambda.HG,PRtydelta.HG,PRtydeltalambda.HG)
plot(0:b, PRty.HG.group, type = "l", pch = 21, col = "black", xlab=expression(paste("Number of group detections, ", t[y])),
     ylab = "Probability Mass Function (PMF)", main = "PMF: Grouped-hypergeometric model",ylim=range_pr_ty, xaxt = "n")
axis(side = 1, at = seq(0, b, by = 0.5))
points(0:b,PRty.HG.group,type = "b", pch = 19, col = "black",lty=2)
lines(0:b,PRtyDelta.HG,type = "l", col = "red",lty=1)
points(0:b,PRtyDelta.HG,pch = 19, col = "red")
# abline(v=ExpVarHG.imperfect.sen.group.test$expectation, col = "red",lty=1)
lines(0:b,PRtyDeltaLambda.HG,type = "l", col = "red",lty=2)
points(0:b,PRtyDeltaLambda.HG,pch = 19, col = "red")
# abline(v=ExpVarHG.imperfect.group.test$expectation, col = "red",lty=2)
lines(0:b,PRtydelta.HG,type = "l", col = "blue",lty=1)
points(0:b,PRtydelta.HG,pch = 19, col = "blue")
# abline(v=ExpVarHG.imperfect.sen.item.test$expectation, col = "blue",lty=1)
lines(0:b,PRtydeltalambda.HG,type = "l", col = "blue",lty=2)
points(0:b,PRtydeltalambda.HG,pch = 19, col = "blue")
# abline(v=ExpVarHG.imperfect.item.test$expectation, col = "blue",lty=2)

legend("topleft",
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[delta]),
         bquote(Imperfect ~ D[delta] ~ "," ~ G[lambda]),
         bquote(Imperfect ~ D[Delta]),
         bquote(Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black", "blue", "blue", "red", "red"),
       lty = c(1, 2, 3, 2, 3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size



# png("Figue1 Binomial Only PMF.png",width = 8,height = 6,units = "in", res = 300)
range_pr_ty <- range(PRty.BN.group,PRtyDelta.BN,PRtyDeltaLambda.BN,PRtydelta.BN,PRtydeltalambda.BN)
plot(0:b, PRty.BN.group, type = "l", pch = 21, col = "black", xlab=expression(paste("Number of group detections ", t[y])),
     ylab = "Probability Mass Function (PMF)", main = "PMF: Grouped-binomial model",ylim=range_pr_ty,xaxt = "n")
axis(side = 1, at = seq(0, b, by = 0.5))
points(0:b,PRty.BN.group,type = "b", pch = 19, col = "black",lty=2)
lines(0:b,PRtyDelta.BN,type = "l", col = "red",lty=1)
points(0:b,PRtyDelta.BN,pch = 19, col = "red")
lines(0:b,PRtyDeltaLambda.BN,type = "l", col = "red",lty=2)
points(0:b,PRtyDeltaLambda.BN,pch = 19, col = "red")
lines(0:b,PRtydelta.BN,type = "l", col = "blue",lty=1)
points(0:b,PRtydelta.BN,pch = 19, col = "blue")
lines(0:b,PRtydeltalambda.BN,type = "l", col = "blue",lty=2)
points(0:b,PRtydeltalambda.BN,pch = 19, col = "blue")

legend("topleft",
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[delta] ),
         bquote(Imperfect ~ D[delta] ~ "," ~ G[lambda]),
         bquote(Imperfect ~ D[Delta] ),
         bquote(Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","blue", "blue", "red", "red"),
       lty = c(1,2,3,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size


FI_Range <- range(FI_case1_PMF_Based_Item_HG,FI_case2_PMF_Based_Item_HG,FI_case4_PMF_Based_Item_HG,
                  FI_case1_PMF_Based_Group_HG,FI_case2_PMF_Based_Group_HG,FI_case4_PMF_Based_Group_HG)

plot(Tx_values, FI_case1_PMF_Based_Item_HG, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Hellinger information (HI)",
     main = "HI: Grouped-hypergeometric model")
lines(Tx_values, FI_case2_PMF_Based_Item_HG, col = "blue", lwd = 2, lty = 1)
lines(Tx_values, FI_case4_PMF_Based_Item_HG, col = "blue", lwd = 2, lty = 2)

lines(Tx_values, FI_case2_PMF_Based_Group_HG, col = "red", lwd = 2, lty = 1)
lines(Tx_values, FI_case4_PMF_Based_Group_HG, col = "red", lwd = 2, lty = 2)


legend(70,0.085,
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[delta] ),
         bquote(Imperfect ~ D[delta] ~ "," ~ G[lambda]),
         bquote(Imperfect ~ D[Delta] ),
         bquote(Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","blue", "blue", "red", "red"),
       lty = c(1,1,2,1,2),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size


FI_Range <- range(FI_case1_PMF_Based_Item_BN,FI_case2_PMF_Based_Item_BN,FI_case4_PMF_Based_Item_BN,
                  FI_case1_PMF_Based_Group_BN,FI_case2_PMF_Based_Group_BN,FI_case4_PMF_Based_Group_BN)

plot(Tx_values, FI_case1_PMF_Based_Item_BN, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Hellinger information (HI)",
     main = "HI: Grouped-binomial model")
lines(Tx_values, FI_case2_PMF_Based_Item_BN, col = "blue", lwd = 2, lty = 1)
lines(Tx_values, FI_case4_PMF_Based_Item_BN, col = "blue", lwd = 2, lty = 2)

lines(Tx_values, FI_case2_PMF_Based_Group_BN, col = "red", lwd = 2, lty = 1)
lines(Tx_values, FI_case4_PMF_Based_Group_BN, col = "red", lwd = 2, lty = 2)


legend(70,0.075,
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[delta] ),
         bquote(Imperfect ~ D[delta] ~ "," ~ G[lambda]),
         bquote(Imperfect ~ D[Delta] ),
         bquote(Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","blue", "blue", "red", "red"),
       lty = c(1,1,2,1,2),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size

dev.off()
