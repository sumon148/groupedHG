# PMF Figures : HG Distribution -------
N <- 100  # Total population size
Tx <- 20  # Number of contaminated items
barN <- 4  # Group size
b <- 4     # Number of groups
delta <- 0.7  # Test sensitivity
lambda <- 0.8  # Test Specificity

ExpVarHG.perfect.item.test <- ExpVarHG.imperfect(N=100, Tx=20, barN=4, b=4, delta=1, lambda=1,type="item")
ExpVarHG.imperfect.sen.item.test <- ExpVarHG.imperfect(N=100, Tx=20, barN=4, b=4, delta=0.70, lambda=1,type="item")
ExpVarHG.imperfect.item.test <- ExpVarHG.imperfect(N=100, Tx=20, barN=4, b=4, delta=0.70, lambda=0.8,type="item")

ExpVarHG.perfect.group.test <- ExpVarHG.imperfect(N=100, Tx=20, barN=4, b=4, delta=1, lambda=1,type="group")
ExpVarHG.imperfect.sen.group.test <- ExpVarHG.imperfect(N=100, Tx=20, barN=4, b=4, delta=0.70, lambda=1,type="group")
ExpVarHG.imperfect.group.test <- ExpVarHG.imperfect(N=100, Tx=20, barN=4, b=4, delta=0.70, lambda=0.8,type="group")


ty.values <- c(0:b)
# Different cases: Group level sensitivity
PRty.HG.group <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N, barN, Tx, b, delta=1.0, lambda=1.0))
PRtyDelta.HG <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N, barN, Tx, b, delta=delta, lambda=1.0))
PRtyDeltaLambda.HG <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N, barN, Tx, b, delta=delta, lambda=lambda))

# Different cases: Item level sensitivity
PRty.HG.item <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N, barN, Tx, b, delta=1.0, lambda=1.0))
PRtydelta.HG <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N, barN, Tx, b, delta=delta, lambda=1.0))
PRtydeltalambda.HG <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N, barN, Tx, b, delta=delta, lambda=lambda))

png("Figue1 HG Only PMF.png",width = 8,height = 6,units = "in", res = 300)
range_pr_ty <- range(PRty.HG.group,PRtyDelta.HG,PRtyDeltaLambda.HG,PRtydelta.HG,PRtydeltalambda.HG)
plot(0:b, PRty.HG.group, type = "l", pch = 21, col = "black", xlab=expression(paste("Number of group detections, ", t[y])),
     ylab = "Probability mass", main = "Grouped Hypergeometric model",ylim=range_pr_ty, xaxt = "n")
axis(side = 1, at = seq(0, b, by = 0.5))
points(0:b,PRty.HG.group,type = "b", pch = 19, col = "black",lty=2)
abline(v=ExpVarHG.perfect.group.test$expectation,col="black",lty=1)
# lines(0:b,pmf_perfect_test_group_binomial,type = "p", pch = 19, col = "black",cex=0.75)
lines(0:b,PRtyDelta.HG,type = "l", col = "red",lty=1)
points(0:b,PRtyDelta.HG,pch = 19, col = "red")
abline(v=ExpVarHG.imperfect.sen.group.test$expectation, col = "red",lty=1)
lines(0:b,PRtyDeltaLambda.HG,type = "l", col = "red",lty=2)
points(0:b,PRtyDeltaLambda.HG,pch = 19, col = "red")
abline(v=ExpVarHG.imperfect.group.test$expectation, col = "red",lty=2)
lines(0:b,PRtydelta.HG,type = "l", col = "blue",lty=1)
points(0:b,PRtydelta.HG,pch = 19, col = "blue")
abline(v=ExpVarHG.imperfect.sen.item.test$expectation, col = "blue",lty=1)
lines(0:b,PRtydeltalambda.HG,type = "l", col = "blue",lty=2)
points(0:b,PRtydeltalambda.HG,pch = 19, col = "blue")
abline(v=ExpVarHG.imperfect.item.test$expectation, col = "blue",lty=2)

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


dev.off()



# PMF Figures : BN Distribution -------

N <- 100  # Total population size
Tx <- 20  # Number of contaminated items
barN <- 4  # Group size
b <- 4     # Number of groups
delta <- 0.7  # Test sensitivity
lambda <- 0.8  # Test Specificity

ExpVarBN.perfect.item.test <- ExpVarBN.imperfect(N=100, Tx=20, barN=4, b=4, delta=1, lambda=1,type="item")
ExpVarBN.imperfect.sen.item.test <- ExpVarBN.imperfect(N=100, Tx=20, barN=4, b=4, delta=0.70, lambda=1,type="item")
ExpVarBN.imperfect.item.test <- ExpVarBN.imperfect(N=100, Tx=20, barN=4, b=4, delta=0.70, lambda=0.8,type="item")

ExpVarBN.perfect.group.test <- ExpVarBN.imperfect(N=100, Tx=20, barN=4, b=4, delta=1, lambda=1,type="group")
ExpVarBN.imperfect.sen.group.test <- ExpVarBN.imperfect(N=100, Tx=20, barN=4, b=4, delta=0.70, lambda=1,type="group")
ExpVarBN.imperfect.group.test <- ExpVarBN.imperfect(N=100, Tx=20, barN=4, b=4, delta=0.70, lambda=0.8,type="group")

ty.values <- c(0:b)
# Different cases: Group level sensitivity
PRty.BN.group <- sapply(ty.values, function(ty) pmfBN.imperfect.group(ty, Tx, b, barN, N, delta=1.0, lambda=1.0))
PRtyDelta.BN <- sapply(ty.values, function(ty) pmfBN.imperfect.group(ty, Tx, b, barN, N, delta=delta, lambda=1.0))
PRtyDeltaLambda.BN <- sapply(ty.values, function(ty) pmfBN.imperfect.group(ty, Tx, b, barN, N, delta=delta, lambda=lambda))

# Different cases: Item level sensitivity
PRty.BN.item <- sapply(ty.values, function(ty) pmfBN.imperfect.item(ty, Tx, b, barN, N, delta=1.0, lambda=1.0))
PRtydelta.BN <- sapply(ty.values, function(ty) pmfBN.imperfect.item(ty, Tx, b, barN, N, delta=delta, lambda=1.0))
PRtydeltalambda.BN <- sapply(ty.values, function(ty) pmfBN.imperfect.item(ty, Tx, b, barN, N, delta=delta, lambda=lambda))

png("Figue1 Binomial Only PMF.png",width = 8,height = 6,units = "in", res = 300)
range_pr_ty <- range(PRty.BN.group,PRtyDelta.BN,PRtyDeltaLambda.BN,PRtydelta.BN,PRtydeltalambda.BN)
plot(0:b, PRty.BN.group, type = "l", pch = 21, col = "black", xlab=expression(paste("Number of group detections ", t[y])),
     ylab = "Probability mass", main = "Grouped Binomial model",ylim=range_pr_ty,xaxt = "n")
axis(side = 1, at = seq(0, b, by = 0.5))
points(0:b,PRty.BN.group,type = "b", pch = 19, col = "black",lty=2)
abline(v=ExpVarBN.perfect.group.test$expectation,col="black",lty=1)
# lines(0:b,pmf_perfect_test_group_binomial,type = "p", pch = 19, col = "black",cex=0.75)
lines(0:b,PRtyDelta.BN,type = "l", col = "red",lty=1)
points(0:b,PRtyDelta.BN,pch = 19, col = "red")
abline(v=ExpVarBN.imperfect.sen.group.test$expectation, col = "red",lty=1)
lines(0:b,PRtyDeltaLambda.BN,type = "l", col = "red",lty=2)
points(0:b,PRtyDeltaLambda.BN,pch = 19, col = "red")
abline(v=ExpVarBN.imperfect.group.test$expectation, col = "red",lty=2)
lines(0:b,PRtydelta.BN,type = "l", col = "blue",lty=1)
points(0:b,PRtydelta.BN,pch = 19, col = "blue")
abline(v=ExpVarBN.imperfect.sen.item.test$expectation, col = "blue",lty=1)
lines(0:b,PRtydeltalambda.BN,type = "l", col = "blue",lty=2)
points(0:b,PRtydeltalambda.BN,pch = 19, col = "blue")
abline(v=ExpVarBN.imperfect.item.test$expectation, col = "blue",lty=2)

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

dev.off()



# FI Figures : HG Distribution using Analytic Derivative (AD) Form ------------------
# As per the paper : HG Distribution  ---------#

Tx_values <- c(0:80)
# Different cases: Item level sensitivity
FI_case1_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "item"))
FI_case2_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "item"))
FI_case4_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "item"))
# Different cases: Group level sensitivity
FI_case1_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "group"))
FI_case2_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "group"))
FI_case4_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "group"))


png("Figure/Fig 2 HG AD Based.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,3))


# Plot results: Group level

FI_Range <- range(FI_case1_AD_Based_Item,FI_case2_AD_Based_Item,FI_case4_AD_Based_Item,
                  FI_case1_AD_Based_Group,FI_case2_AD_Based_Group,FI_case4_AD_Based_Group)

plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD based FI: All Cases under HG")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_AD_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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



# Plot results: Item level
# Plot results: Group level

FI_Range <- c(0,0.05)

plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD based FI: All Cases under HG")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_AD_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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



# Plot results: Group level

FI_Range <- c(0,0.02)

plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD based FI: All Cases under HG")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_AD_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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



dev.off()



# FI Figures : HG Distribution using Probability Mass Function (PMF) Based ------------------
# As per the paper : HG Distribution  ---------#


Tx_values <- c(0:80)
# Different cases: Item level sensitivity
FI_case1_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "item"))
FI_case2_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "item"))
FI_case4_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "item"))
# Different cases: Group level sensitivity
FI_case1_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "group"))
FI_case2_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "group"))
FI_case4_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "group"))


png("Figure/Fig 2 HG PMF Based.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,3))


# Plot results: Group level

FI_Range <- range(FI_case1_PMF_Based_Item,FI_case2_PMF_Based_Item,FI_case4_PMF_Based_Item,
                  FI_case1_PMF_Based_Group,FI_case2_PMF_Based_Group,FI_case4_PMF_Based_Group)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "PMF Based FI: All Cases under HG")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_PMF_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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



# Plot results: Item level
# Plot results: Group level

FI_Range <- c(0,0.05)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "PMF based FI: All Cases under HG")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_PMF_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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



# Plot results: Group level

FI_Range <- c(0,0.005)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "PMF based FI: All Cases under HG")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_PMF_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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



dev.off()






# FI Figures : BN Distribution using Analytic Derivative (AD) Form ------------------

Tx_values <- c(0:80)

# Different cases: Item level sensitivity
FI_case1_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "item"))
FI_case2_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "item"))
FI_case4_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "item"))
# Different cases: Group level sensitivity
FI_case1_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "group"))
FI_case2_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "group"))
FI_case4_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "group"))



png("Figure/Fig 2 BN AD Based.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,3))


# Plot results: Group level

FI_Range <- range(FI_case1_AD_Based_Item,FI_case2_AD_Based_Item,FI_case4_AD_Based_Item,
                  FI_case1_AD_Based_Group,FI_case2_AD_Based_Group,FI_case4_AD_Based_Group)

plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "Fisher Information: All Cases under binomial approximation")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_AD_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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



FI_Range <- c(0,0.02)

plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "Fisher Information: All Cases under binomial approximation")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_AD_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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


FI_Range <- c(0,0.005)

plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "Fisher Information: All Cases under binomial approximation")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_AD_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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



dev.off()



# FI Figures : BN Distribution using Analytic Derivative (PMF) Form ------------------


Tx_values <- c(0:80)

# Different cases: Item level sensitivity
FI_case1_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "item"))
FI_case2_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "item"))
FI_case4_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "item"))
# Different cases: Group level sensitivity
FI_case1_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "group"))
FI_case2_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "group"))
FI_case4_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "group"))



png("Figure/Fig 2 BN PMF Based.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,3))


# Plot results: Group level

FI_Range <- range(FI_case1_PMF_Based_Item,FI_case2_PMF_Based_Item,FI_case4_PMF_Based_Item,
                  FI_case1_PMF_Based_Group,FI_case2_PMF_Based_Group,FI_case4_PMF_Based_Group)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "Fisher Information: All Cases under binomial approximation")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_PMF_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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



FI_Range <- c(0,0.02)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "Fisher Information: All Cases under binomial approximation")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_PMF_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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


FI_Range <- c(0,0.005)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "Fisher Information: All Cases under binomial approximation")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_PMF_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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



dev.off()





# FI for HG Distribution  : AD vs PMF ---------#


Tx_values <- c(0:80)
# Different cases: Item level sensitivity
FI_case1_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "item"))
FI_case2_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "item"))
FI_case4_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "item"))
# Different cases: Group level sensitivity
FI_case1_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "group"))
FI_case2_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "group"))
FI_case4_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "group"))

# Different cases: Item level sensitivity
FI_case1_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "item"))
FI_case2_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "item"))
FI_case4_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "item"))
# Different cases: Group level sensitivity
FI_case1_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "group"))
FI_case2_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "group"))
FI_case4_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "group"))

png("Figure/Fig 2 HG AD PMF.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))


# Plot results: Group level

FI_Range <- range(FI_case1_AD_Based_Item,FI_case2_AD_Based_Item,FI_case4_AD_Based_Item,
                  FI_case1_AD_Based_Group,FI_case2_AD_Based_Group,FI_case4_AD_Based_Group)

FI_Range <- range(FI_case1_PMF_Based_Item,FI_case2_PMF_Based_Item,FI_case4_PMF_Based_Item,
                  FI_case1_PMF_Based_Group,FI_case2_PMF_Based_Group,FI_case4_PMF_Based_Group)

plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD based FI: All Cases under HG")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_AD_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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




FI_Range <- range(FI_case1_PMF_Based_Item,FI_case2_PMF_Based_Item,FI_case4_PMF_Based_Item,
                  FI_case1_PMF_Based_Group,FI_case2_PMF_Based_Group,FI_case4_PMF_Based_Group)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "PMF based FI: All Cases under HG")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_PMF_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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




dev.off()



# FI for BN Distribution  : AD vs PMF ---------#


Tx_values <- c(0:80)
# Different cases: Item level sensitivity
FI_case1_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "item"))
FI_case2_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "item"))
FI_case4_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "item"))
# Different cases: Group level sensitivity
FI_case1_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "group"))
FI_case2_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "group"))
FI_case4_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "group"))

# Different cases: Item level sensitivity
FI_case1_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "item"))
FI_case2_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "item"))
FI_case4_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "item"))
# Different cases: Group level sensitivity
FI_case1_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "group"))
FI_case2_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "group"))
FI_case4_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfBN.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "group"))

png("Figure/Fig 2 BN AD PMF.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))


# Plot results: Group level

FI_Range <- range(FI_case1_AD_Based_Item,FI_case2_AD_Based_Item,FI_case4_AD_Based_Item,
                  FI_case1_AD_Based_Group,FI_case2_AD_Based_Group,FI_case4_AD_Based_Group)

FI_Range <- range(FI_case1_PMF_Based_Item,FI_case2_PMF_Based_Item,FI_case4_PMF_Based_Item,
                  FI_case1_PMF_Based_Group,FI_case2_PMF_Based_Group,FI_case4_PMF_Based_Group)


plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD based FI: All Cases under BN")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_AD_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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




FI_Range <- range(FI_case1_PMF_Based_Item,FI_case2_PMF_Based_Item,FI_case4_PMF_Based_Item,
                  FI_case1_PMF_Based_Group,FI_case2_PMF_Based_Group,FI_case4_PMF_Based_Group)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "PMF based FI: All Cases under BN")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_PMF_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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




dev.off()



png("Figure/Fig 2 BN AD PMF BB.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))


# Plot results: Group level

FI_Range <- range(FI_case1_AD_Based_Item,FI_case2_AD_Based_Item,FI_case4_AD_Based_Item,
                  FI_case1_AD_Based_Group,FI_case2_AD_Based_Group,FI_case4_AD_Based_Group)

FI_Range <- range(FI_case1_PMF_Based_Item,FI_case2_PMF_Based_Item,FI_case4_PMF_Based_Item,
                  FI_case1_PMF_Based_Group,FI_case2_PMF_Based_Group,FI_case4_PMF_Based_Group)

FI_Range <- c(0,0.02)
plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD based FI: All Cases under BN")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_AD_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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




FI_Range <- range(FI_case1_PMF_Based_Item,FI_case2_PMF_Based_Item,FI_case4_PMF_Based_Item,
                  FI_case1_PMF_Based_Group,FI_case2_PMF_Based_Group,FI_case4_PMF_Based_Group)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "PMF based FI: All Cases under BN")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "blue", lwd = 1, lty = 3)

lines(Tx_values, FI_case2_PMF_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
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




dev.off()
