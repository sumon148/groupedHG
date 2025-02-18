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



# Check Probability Functions and their derivatives for lower Tx and Imperfect Specificity:

pmfHG.imperfect.group(ty=0, N=100, barN=4, Tx=0, b=4, delta=1, lambda=1)
pmfHG.imperfect.group(ty=0, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=1)
pmfHG.imperfect.group(ty=0, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=0.7)

pmfHG.imperfect.group(ty=0, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=0.7)
pmfHG.imperfect.group(ty=1, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=0.7)
pmfHG.imperfect.group(ty=2, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=0.7)
pmfHG.imperfect.group(ty=3, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=0.7)
pmfHG.imperfect.group(ty=4, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=0.7)

ty.values <- c(0:b)
PRtyDeltaLambda.HG.Tx.0 <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=0.7))
sum(PRtyDeltaLambda.HG.Tx.0)
dPRtyDeltaLambda.HG.Tx.0 <- sapply(ty.values, function(ty) dpmfHG.imperfect.group(ty, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=0.7))
cbind(PRtyDeltaLambda.HG.Tx.0,dPRtyDeltaLambda.HG.Tx.0)

PRtyDeltaLambda.HG.Tx.1 <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N=100, barN=4, Tx=1, b=4, delta=0.75, lambda=0.7))
sum(PRtyDeltaLambda.HG.Tx.1)
dPRtyDeltaLambda.HG.Tx.1 <- sapply(ty.values, function(ty) dpmfHG.imperfect.group(ty, N=100, barN=4, Tx=1, b=4, delta=0.75, lambda=0.7))
cbind(PRtyDeltaLambda.HG.Tx.1,dPRtyDeltaLambda.HG.Tx.1)

PRtyDeltaLambda.HG.Tx.2 <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N=100, barN=4, Tx=2, b=4, delta=0.75, lambda=0.7))
sum(PRtyDeltaLambda.HG.Tx.2)
dPRtyDeltaLambda.HG.Tx.2 <- sapply(ty.values, function(ty) dpmfHG.imperfect.group(ty, N=100, barN=4, Tx=2, b=4, delta=0.75, lambda=0.7))
cbind(PRtyDeltaLambda.HG.Tx.2,dPRtyDeltaLambda.HG.Tx.2)

PRtyDeltaLambda.HG.Tx.3 <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N=100, barN=4, Tx=3, b=4, delta=0.75, lambda=0.7))
sum(PRtyDeltaLambda.HG.Tx.3)
dPRtyDeltaLambda.HG.Tx.3 <- sapply(ty.values, function(ty) dpmfHG.imperfect.group(ty, N=100, barN=4, Tx=3, b=4, delta=0.75, lambda=0.7))
cbind(PRtyDeltaLambda.HG.Tx.3,dPRtyDeltaLambda.HG.Tx.3)

PRtyDeltaLambda.HG.Tx.4 <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N=100, barN=4, Tx=4, b=4, delta=0.75, lambda=0.7))
sum(PRtyDeltaLambda.HG.Tx.4)
dPRtyDeltaLambda.HG.Tx.4 <- sapply(ty.values, function(ty) dpmfHG.imperfect.group(ty, N=100, barN=4, Tx=4, b=4, delta=0.75, lambda=0.7))
cbind(PRtyDeltaLambda.HG.Tx.4,dPRtyDeltaLambda.HG.Tx.4)

PRtyDeltaLambda.HG.Tx.5 <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N=100, barN=4, Tx=5, b=4, delta=0.75, lambda=0.7))
sum(PRtyDeltaLambda.HG.Tx.5)
dPRtyDeltaLambda.HG.Tx.5 <- sapply(ty.values, function(ty) dpmfHG.imperfect.group(ty, N=100, barN=4, Tx=5, b=4, delta=0.75, lambda=0.7))
cbind(PRtyDeltaLambda.HG.Tx.5,dPRtyDeltaLambda.HG.Tx.5)

PRtyDeltaLambda.HG.Tx.20 <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N=100, barN=4, Tx=20, b=4, delta=0.75, lambda=0.7))
sum(PRtyDeltaLambda.HG.Tx.20)
dPRtyDeltaLambda.HG.Tx.20 <- sapply(ty.values, function(ty) dpmfHG.imperfect.group(ty, N=100, barN=4, Tx=20, b=4, delta=0.75, lambda=0.7))
cbind(PRtyDeltaLambda.HG.Tx.20,dPRtyDeltaLambda.HG.Tx.20)


ty.values <- c(0:b)
PRtydeltalambda.HG.Tx.0 <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=0.7))
sum(PRtydeltalambda.HG.Tx.0)
dPRtydeltalambda.HG.Tx.0 <- sapply(ty.values, function(ty) dpmfHG.imperfect.item(ty, N=100, barN=4, Tx=0, b=4, delta=0.75, lambda=0.7))
cbind(PRtydeltalambda.HG.Tx.0,dPRtydeltalambda.HG.Tx.0)

PRtydeltalambda.HG.Tx.1 <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N=100, barN=4, Tx=1, b=4, delta=0.75, lambda=0.7))
sum(PRtydeltalambda.HG.Tx.1)
dPRtydeltalambda.HG.Tx.1 <- sapply(ty.values, function(ty) dpmfHG.imperfect.item(ty, N=100, barN=4, Tx=1, b=4, delta=0.75, lambda=0.7))
cbind(PRtydeltalambda.HG.Tx.1,dPRtydeltalambda.HG.Tx.1)

PRtydeltalambda.HG.Tx.2 <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N=100, barN=4, Tx=2, b=4, delta=0.75, lambda=0.7))
sum(PRtydeltalambda.HG.Tx.2)
dPRtydeltalambda.HG.Tx.2 <- sapply(ty.values, function(ty) dpmfHG.imperfect.item(ty, N=100, barN=4, Tx=2, b=4, delta=0.75, lambda=0.7))
cbind(PRtydeltalambda.HG.Tx.2,dPRtydeltalambda.HG.Tx.2)

PRtydeltalambda.HG.Tx.3 <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N=100, barN=4, Tx=3, b=4, delta=0.75, lambda=0.7))
sum(PRtydeltalambda.HG.Tx.3)
dPRtydeltalambda.HG.Tx.3 <- sapply(ty.values, function(ty) dpmfHG.imperfect.item(ty, N=100, barN=4, Tx=3, b=4, delta=0.75, lambda=0.7))
cbind(PRtydeltalambda.HG.Tx.3,dPRtydeltalambda.HG.Tx.3)

PRtydeltalambda.HG.Tx.4 <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N=100, barN=4, Tx=4, b=4, delta=0.75, lambda=0.7))
sum(PRtydeltalambda.HG.Tx.4)
dPRtydeltalambda.HG.Tx.4 <- sapply(ty.values, function(ty) dpmfHG.imperfect.item(ty, N=100, barN=4, Tx=4, b=4, delta=0.75, lambda=0.7))
cbind(PRtydeltalambda.HG.Tx.4,dPRtydeltalambda.HG.Tx.4)


PRtydeltalambda.HG.Tx.5 <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N=100, barN=4, Tx=5, b=4, delta=0.75, lambda=0.7))
sum(PRtydeltalambda.HG.Tx.5)
dPRtydeltalambda.HG.Tx.5 <- sapply(ty.values, function(ty) dpmfHG.imperfect.item(ty, N=100, barN=4, Tx=5, b=4, delta=0.75, lambda=0.7))
cbind(PRtydeltalambda.HG.Tx.5,dPRtydeltalambda.HG.Tx.5)

PRtydeltalambda.HG.Tx.20 <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N=100, barN=4, Tx=20, b=4, delta=0.75, lambda=0.7))
sum(PRtydeltalambda.HG.Tx.20)
dPRtydeltalambda.HG.Tx.20 <- sapply(ty.values, function(ty) dpmfHG.imperfect.item(ty, N=100, barN=4, Tx=20, b=4, delta=0.75, lambda=0.7))
cbind(PRtydeltalambda.HG.Tx.20,dPRtydeltalambda.HG.Tx.20)






ty.values <- c(0:b)
# Different cases: Group level sensitivity
PRty.HG.group <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N, barN, Tx, b, delta=1.0, lambda=1.0))
PRtyDelta.HG <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N, barN, Tx, b, delta=delta, lambda=1.0))
PRtyDeltaLambda.HG <- sapply(ty.values, function(ty) pmfHG.imperfect.group(ty, N, barN, Tx, b, delta=delta, lambda=lambda))

# Different cases: Item level sensitivity
PRty.HG.item <- sapply(ty.values, function(ty) pmfHG.imperfect.item(ty, N, barN, Tx, b, delta=1.0, lambda=1.0))
PRty.HG.item.perfect <- sapply(ty.values, function(ty) pmfHG.perfect(ty, N, barN, Tx, b))
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
FI_case1_AD_Based_Perfect <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.perfect(Tx, N=100, b=4, barN=4, method = "AD"))
FI_case2_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "item",verbose = FALSE))
FI_case4_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "item",verbose = FALSE))
# Different cases: Group level sensitivity
FI_case1_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "group",verbose = FALSE))
FI_case2_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "group",verbose = FALSE))
FI_case4_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "group",verbose = FALSE))

# Using detail estimation by ty values

Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_case1_AD_Based_Item <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N = 100, b = 4, barN = 4, delta = 1.0, lambda = 1.0, method = "AD", type = "item")

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_case1_AD_Based_Item_detail <- bind_rows(FI_results_case1_AD_Based_Item)

FI_case1_AD_Based_Item_detail <- data.frame(cbind(Tx=FI_case1_AD_Based_Item_detail$Tx,ty=FI_case1_AD_Based_Item_detail$ty,p_ty=FI_case1_AD_Based_Item_detail$p_ty,dp_ty=FI_case1_AD_Based_Item_detail$dp_ty,FI_Tx=FI_case1_AD_Based_Item_detail$FI_Tx))
View(FI_case1_AD_Based_Item_detail)

tapply(FI_case1_AD_Based_Item_detail$p_ty,FI_case1_AD_Based_Item_detail$Tx,sum)


Tx_values <- c(0:100)
# Apply function over all Tx values
FI_results_case3_AD_Based_item <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N = 100, b = 4, barN = 4, delta = 0.7, lambda = 0.8, method = "AD", type = "item")

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_case3_AD_Based_item <- bind_rows(FI_results_case3_AD_Based_item)

FI_results_case3_AD_Based_item <- data.frame(cbind(Tx=FI_results_case3_AD_Based_item$Tx,ty=FI_results_case3_AD_Based_item$ty,p_ty=FI_results_case3_AD_Based_item$p_ty,dp_ty=FI_results_case3_AD_Based_item$dp_ty,FI_Tx=FI_results_case3_AD_Based_item$FI_Tx))
View(FI_results_case3_AD_Based_item)

tapply(FI_results_case3_AD_Based_item$p_ty,FI_results_case3_AD_Based_item$Tx,sum)



Tx_values <- c(0:100)
# Apply function over all Tx values
FI_results_case3_PMF_Based_item <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N = 100, b = 4, barN = 4, delta = 0.7, lambda = 0.8, method = "PMF", type = "item")

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_case3_PMF_Based_item <- bind_rows(FI_results_case3_PMF_Based_item)

FI_results_case3_PMF_Based_item <- data.frame(cbind(Tx=FI_results_case3_PMF_Based_item$Tx,ty=FI_results_case3_PMF_Based_item$ty,p_ty=FI_results_case3_PMF_Based_item$p_ty,dp_ty=FI_results_case3_PMF_Based_item$dp_ty,FI_Tx=FI_results_case3_PMF_Based_item$FI_Tx))
View(FI_results_case3_PMF_Based_item)

tapply(FI_results_case3_PMF_Based_item$p_ty,FI_results_case3_PMF_Based_item$Tx,sum)


plot(Tx_values,tapply(FI_results_case3_PMF_Based_item$FI_Tx,FI_results_case3_PMF_Based_item$Tx,mean))
plot(Tx_values,tapply(FI_results_case3_AD_Based_item$FI_Tx,FI_results_case3_AD_Based_item$Tx,mean))

# Using detail estimation by ty values

Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_case1_AD_Based_group <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N = 100, b = 4, barN = 4, delta = 1.0, lambda = 1.0, method = "AD", type = "group")

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_case1_AD_Based_group_detail <- bind_rows(FI_results_case1_AD_Based_group)

FI_case1_AD_Based_group_detail <- data.frame(cbind(Tx=FI_case1_AD_Based_group_detail$Tx,ty=FI_case1_AD_Based_group_detail$ty,p_ty=FI_case1_AD_Based_group_detail$p_ty,dp_ty=FI_case1_AD_Based_group_detail$dp_ty,FI_Tx=FI_case1_AD_Based_group_detail$FI_Tx))
View(FI_case1_AD_Based_group_detail)

tapply(FI_case1_AD_Based_group_detail$p_ty,FI_case1_AD_Based_group_detail$Tx,sum)





# Using detail estimation by ty values

Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_case2_AD_Based_group <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N = 100, b = 4, barN = 4, delta = 0.7, lambda = 1.0, method = "AD", type = "group")

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_case2_AD_Based_group <- bind_rows(FI_results_case2_AD_Based_group)

FI_results_case2_AD_Based_group <- data.frame(cbind(Tx=FI_results_case2_AD_Based_group$Tx,ty=FI_results_case2_AD_Based_group$ty,p_ty=FI_results_case2_AD_Based_group$p_ty,dp_ty=FI_results_case2_AD_Based_group$dp_ty,FI_Tx=FI_results_case2_AD_Based_group$FI_Tx))
View(FI_case1_AD_Based_group_detail)

tapply(FI_results_case2_AD_Based_group$p_ty,FI_results_case2_AD_Based_group$Tx,sum)



# Using detail estimation by ty values

Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_case3_AD_Based_group <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N = 100, b = 4, barN = 4, delta = 0.7, lambda = 0.8, method = "AD", type = "group")

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_case3_AD_Based_group <- bind_rows(FI_results_case3_AD_Based_group)

FI_results_case3_AD_Based_group <- data.frame(cbind(Tx=FI_results_case3_AD_Based_group$Tx,ty=FI_results_case3_AD_Based_group$ty,p_ty=FI_results_case3_AD_Based_group$p_ty,dp_ty=FI_results_case3_AD_Based_group$dp_ty,FI_Tx=FI_results_case3_AD_Based_group$FI_Tx))
View(FI_results_case3_AD_Based_group)

tapply(FI_results_case3_AD_Based_group$p_ty,FI_results_case3_AD_Based_group$Tx,sum)



plot(FI_case1_AD_Based_group_detail$Tx,FI_case1_AD_Based_group_detail$FI_Tx)
lines(FI_results_case2_AD_Based_group$Tx,FI_results_case2_AD_Based_group$FI_Tx)



# Check FI for Tx=0 based on AD and PMF based ----------------#

p_ty_0 <- pmfHG.imperfect.item(ty=0, N=100, barN=4, Tx=0, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
p_ty_1 <- pmfHG.imperfect.item(ty=1, N=100, barN=4, Tx=0, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
p_ty_2 <- pmfHG.imperfect.item(ty=2, N=100, barN=4, Tx=0, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
p_ty_3 <- pmfHG.imperfect.item(ty=3, N=100, barN=4, Tx=0, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
p_ty_4 <- pmfHG.imperfect.item(ty=4, N=100, barN=4, Tx=0, b=4, delta=0.7, lambda=0.8, verbose=TRUE)

p_ty_0_tx_1 <- pmfHG.imperfect.item(ty=0, N=100, barN=4, Tx=1, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
p_ty_1_tx_1 <- pmfHG.imperfect.item(ty=1, N=100, barN=4, Tx=1, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
p_ty_2_tx_1 <- pmfHG.imperfect.item(ty=2, N=100, barN=4, Tx=1, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
p_ty_3_tx_1 <- pmfHG.imperfect.item(ty=3, N=100, barN=4, Tx=1, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
p_ty_4_tx_1 <- pmfHG.imperfect.item(ty=4, N=100, barN=4, Tx=1, b=4, delta=0.7, lambda=0.8, verbose=TRUE)

4 * ((sqrt(p_ty_0_tx_1)-sqrt(p_ty_0))^2+
(sqrt(p_ty_1_tx_1)-sqrt(p_ty_1))^2+
(sqrt(p_ty_2_tx_1)-sqrt(p_ty_2))^2+
(sqrt(p_ty_3_tx_1)-sqrt(p_ty_3))^2)


dp_ty_0 <- dpmfHG.imperfect.item(ty=0, N=100, barN=4, Tx=0, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
dp_ty_1 <- dpmfHG.imperfect.item(ty=1, N=100, barN=4, Tx=0, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
dp_ty_2 <- dpmfHG.imperfect.item(ty=2, N=100, barN=4, Tx=0, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
dp_ty_3 <- dpmfHG.imperfect.item(ty=3, N=100, barN=4, Tx=0, b=4, delta=0.7, lambda=0.8, verbose=TRUE)
dp_ty_4 <- dpmfHG.imperfect.item(ty=4, N=100, barN=4, Tx=0, b=4, delta=0.7, lambda=0.8, verbose=TRUE)

sum(c(dp_ty_0^2/p_ty_0,dp_ty_1^2/p_ty_1,dp_ty_2^2/p_ty_2,dp_ty_3^2/p_ty_3,dp_ty_4^2/p_ty_4))

png("Figure/Fig 2 HG AD Based.png",width = 8,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))


# Plot results: Group level

FI_Range <- range(FI_case1_AD_Based_Item,FI_case2_AD_Based_Item,FI_case4_AD_Based_Item,
                  FI_case1_AD_Based_Group,FI_case2_AD_Based_Group,FI_case4_AD_Based_Group)

plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD based FI: All Cases under HG")
# lines(Tx_values, FI_case1_AD_Based_Perfect, col = "black", lwd = 1, lty = 3)
# lines(Tx_values, FI_case1_AD_Based_Group, col = "black", lwd = 1, lty = 4)
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




dev.off()



# FI Figures : HG Distribution using Probability Mass Function (PMF) Based ------------------
# As per the paper : HG Distribution  ---------#


Tx_values <- c(0:100)
# Different cases: Item level sensitivity
FI_case1_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "item"))
FI_case2_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "item"))
FI_case4_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "item"))
# Different cases: Group level sensitivity
FI_case1_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "group"))
FI_case2_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "group"))
FI_case4_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "group"))


png("Figure/Fig 2 HG PMF Based.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))


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





dev.off()



# FI Figures : HG Distribution using Analytic Derivative (AD) and PMF Based FI ------------------
# As per the paper : HG Distribution  ---------#

Tx_values <- c(0:80)
# Different cases: Item level sensitivity
FI_case1_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "item"))
FI_case1_AD_Based_Perfect <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.perfect(Tx, N=100, b=4, barN=4, method = "AD"))
FI_case2_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "item",verbose = FALSE))
FI_case4_AD_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "item",verbose = FALSE))
# Different cases: Group level sensitivity
FI_case1_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "AD", type = "group",verbose = FALSE))
FI_case2_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "AD", type = "group",verbose = FALSE))
FI_case4_AD_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "AD", type = "group",verbose = FALSE))


Tx_values <- c(0:80)
# Different cases: Item level sensitivity
FI_case1_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "item"))
FI_case2_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "item"))
FI_case4_PMF_Based_Item <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "item"))
# Different cases: Group level sensitivity
FI_case1_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=1.0, lambda=1.0, method = "PMF", type = "group"))
FI_case2_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=1.0, method = "PMF", type = "group"))
FI_case4_PMF_Based_Group <- sapply(Tx_values, function(Tx) FIpmfHG.Tx.imperfect(Tx, N=100, b=4, barN=4, delta=0.7, lambda=0.8, method = "PMF", type = "group"))




png("Figure/Fig 2 HG AD PMF Based Item.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(2,2))


# Plot results: Item level

FI_Range <- range(FI_case1_PMF_Based_Item,FI_case2_PMF_Based_Item,FI_case4_PMF_Based_Item,
                  FI_case1_AD_Based_Item,FI_case2_AD_Based_Item,FI_case4_AD_Based_Item)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "PMF Based FI: Item level Cases under HG")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "red", pch = 19, cex = 0.3)

legend("topright",
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[delta] ),
         bquote(Imperfect ~ D[delta] ~ "," ~ G[lambda])
         ),
       col = c("black","blue", "red"),
       lty = c(1,2,3,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size


plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD Based FI: All Cases under HG")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "red", pch = 19, cex = 0.3)

legend("topright",
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[delta] ),
         bquote(Imperfect ~ D[delta] ~ "," ~ G[lambda])
       ),
       col = c("black","blue", "red"),
       lty = c(1,2,3,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size



# Plot results: Item level
# Plot results: Group level

FI_Range <- c(0,0.02)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "PMF Based FI: All Cases under HG")
lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Item, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Item, col = "red", pch = 19, cex = 0.3)

legend("topright",
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[delta] ),
         bquote(Imperfect ~ D[delta] ~ "," ~ G[lambda])
       ),
       col = c("black","blue", "red"),
       lty = c(1,2,3,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size


plot(Tx_values, FI_case1_AD_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD Based FI: All Cases under HG")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Item, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Item, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Item, col = "red", pch = 19, cex = 0.3)

legend("topright",
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[delta] ),
         bquote(Imperfect ~ D[delta] ~ "," ~ G[lambda])
       ),
       col = c("black","blue", "red"),
       lty = c(1,2,3,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size




dev.off()




png("Figure/Fig 2 HG AD PMF Based Item Alternative.png",width =12,height = 8,units = "in", res = 300)

par(mfrow=c(2,2))


# Plot results: Item level

FI_Range <- range(FI_case1_PMF_Based_Item,FI_case2_PMF_Based_Item,FI_case4_PMF_Based_Item,
                  FI_case1_AD_Based_Item,FI_case2_AD_Based_Item,FI_case4_AD_Based_Item)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD and PMF Based FI: Item level Cases under HG")
lines(Tx_values, FI_case1_AD_Based_Item, col = "black", lwd = 1, lty = 2)

lines(Tx_values, FI_case2_PMF_Based_Item, col = "blue", lwd = 1, lty = 1)
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)

lines(Tx_values, FI_case4_PMF_Based_Item, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Item, col = "red", lwd = 1, lty = 3)

legend("topright",
       legend = c(
         "Perfect Case: PMF",
         "Perfect Case: AD",
         bquote(PMF: Imperfect ~ D[delta] ),
         bquote(AD: Imperfect ~ D[delta] ),
         bquote(PMF:Imperfect ~ D[delta] ~ "," ~ G[lambda]),
         bquote(AD:Imperfect ~ D[delta] ~ "," ~ G[lambda])
       ),
       col = c("black","black","blue", "blue", "red", "red"),
       lty = c(1,2,1,2,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size




FI_Range <- range(FI_case1_PMF_Based_Item,
                  FI_case1_AD_Based_Item)

plot(Tx_values, FI_case1_PMF_Based_Item, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD and PMF Based FI: Perfect cases under HG")
lines(Tx_values, FI_case1_AD_Based_Item, col = "black", lwd = 1, lty = 2)

legend("topright",
       legend = c(
         "Perfect Case: PMF",
         "Perfect Case: AD",
         bquote(PMF: Imperfect ~ D[delta] ),
         bquote(AD: Imperfect ~ D[delta] ),
         bquote(PMF:Imperfect ~ D[delta] ~ "," ~ G[lambda]),
         bquote(AD:Imperfect ~ D[delta] ~ "," ~ G[lambda])
       ),
       col = c("black","black","blue", "blue", "red", "red"),
       lty = c(1,2,1,2,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size

FI_Range <- range(FI_case2_PMF_Based_Item,FI_case2_AD_Based_Item)

plot(Tx_values, FI_case2_PMF_Based_Item, type = "l", col = "blue", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD and PMF Based FI: Imperfect Sensitivity under HG")
lines(Tx_values, FI_case2_AD_Based_Item, col = "blue", lwd = 1, lty = 2)

legend("topright",
       legend = c(
         "Perfect Case: PMF",
         "Perfect Case: AD",
         bquote(PMF: Imperfect ~ D[delta] ),
         bquote(AD: Imperfect ~ D[delta] ),
         bquote(PMF:Imperfect ~ D[delta] ~ "," ~ G[lambda]),
         bquote(AD:Imperfect ~ D[delta] ~ "," ~ G[lambda])
       ),
       col = c("black","black","blue", "blue", "red", "red"),
       lty = c(1,2,1,2,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size

FI_Range <- range(FI_case4_PMF_Based_Item,FI_case4_AD_Based_Item)

plot(Tx_values, FI_case4_PMF_Based_Item, type = "l", col = "red", lwd = 1, ylim = FI_Range,lty=2,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD and PMF Based FI: Imperfect case under HG")
lines(Tx_values, FI_case4_AD_Based_Item, col = "red", lwd = 1, lty = 3)

legend("topright",
       legend = c(
         "Perfect Case: PMF",
         "Perfect Case: AD",
         bquote(PMF: Imperfect ~ D[delta] ),
         bquote(AD: Imperfect ~ D[delta] ),
         bquote(PMF:Imperfect ~ D[delta] ~ "," ~ G[lambda]),
         bquote(AD:Imperfect ~ D[delta] ~ "," ~ G[lambda])
       ),
       col = c("black","black","blue", "blue", "red", "red"),
       lty = c(1,2,1,2,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size


dev.off()




png("Figure/Fig 2 HG AD PMF Based Group.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(2,2))


# Plot results: Group level

FI_Range <- range(FI_case1_PMF_Based_Group,FI_case2_PMF_Based_Group,FI_case4_PMF_Based_Group,
                  FI_case1_AD_Based_Group,FI_case2_AD_Based_Group,FI_case4_AD_Based_Group)

plot(Tx_values, FI_case1_PMF_Based_Group, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "PMF Based FI: Group level Cases under HG")
lines(Tx_values, FI_case2_PMF_Based_Group, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Group, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[Delta] ),
         bquote(Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","blue", "red"),
       lty = c(1,2,3,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size


plot(Tx_values, FI_case1_AD_Based_Group, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD Based FI: Group level under HG")
lines(Tx_values, FI_case2_AD_Based_Group, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Group, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[Delta] ),
         bquote(Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","blue", "red"),
       lty = c(1,2,3,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size




FI_Range <- c(0,0.02)

plot(Tx_values, FI_case1_PMF_Based_Group, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "PMF Based FI: Group level Cases under HG")
lines(Tx_values, FI_case2_PMF_Based_Group, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_PMF_Based_Group, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_PMF_Based_Group, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_PMF_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[Delta] ),
         bquote(Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","blue", "red"),
       lty = c(1,2,3,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size


plot(Tx_values, FI_case1_AD_Based_Group, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD Based FI: Group level Cases under HG")
lines(Tx_values, FI_case2_AD_Based_Group, col = "blue", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

# Overlay tiny dots at data points
points(Tx_values, FI_case1_AD_Based_Group, col = "black", pch = 19, cex = 0.3)
points(Tx_values, FI_case2_AD_Based_Group, col = "blue", pch = 19, cex = 0.3)
points(Tx_values, FI_case4_AD_Based_Group, col = "red", pch = 19, cex = 0.3)

legend("topright",
       legend = c(
         "Perfect Case",
         bquote(Imperfect ~ D[delta] ),
         bquote(Imperfect ~ D[delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","blue", "red"),
       lty = c(1,2,3,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size




dev.off()


png("Figure/Fig 2 HG AD PMF Based Group Alternative.png",width =12,height = 8,units = "in", res = 300)

par(mfrow=c(2,2))


# Plot results: Group level

FI_Range <- range(FI_case1_PMF_Based_Group,FI_case2_PMF_Based_Group,FI_case4_PMF_Based_Group,
                  FI_case1_AD_Based_Group,FI_case2_AD_Based_Group,FI_case4_AD_Based_Group)

plot(Tx_values, FI_case1_PMF_Based_Group, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD and PMF Based FI: Group level Cases under HG")
lines(Tx_values, FI_case1_AD_Based_Group, col = "black", lwd = 1, lty = 2)

lines(Tx_values, FI_case2_PMF_Based_Group, col = "blue", lwd = 1, lty = 1)
lines(Tx_values, FI_case2_AD_Based_Group, col = "blue", lwd = 1, lty = 2)

lines(Tx_values, FI_case4_PMF_Based_Group, col = "red", lwd = 1, lty = 2)
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

legend("topright",
       legend = c(
         "Perfect Case: PMF",
         "Perfect Case: AD",
         bquote(PMF: Imperfect ~ D[Delta] ),
         bquote(AD: Imperfect ~ D[Delta] ),
         bquote(PMF:Imperfect ~ D[Delta] ~ "," ~ G[Lambda]),
         bquote(AD:Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","black","blue", "blue", "red", "red"),
       lty = c(1,2,1,2,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size




FI_Range <- range(FI_case1_PMF_Based_Group,
                  FI_case1_AD_Based_Group)

plot(Tx_values, FI_case1_PMF_Based_Group, type = "l", col = "black", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD and PMF Based FI: Perfect cases under HG")
lines(Tx_values, FI_case1_AD_Based_Group, col = "black", lwd = 1, lty = 2)

legend("topright",
       legend = c(
         "Perfect Case: PMF",
         "Perfect Case: AD",
         bquote(PMF: Imperfect ~ D[Delta] ),
         bquote(AD: Imperfect ~ D[Delta] ),
         bquote(PMF:Imperfect ~ D[Delta] ~ "," ~ G[Lambda]),
         bquote(AD:Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","black","blue", "blue", "red", "red"),
       lty = c(1,2,1,2,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size

FI_Range <- range(FI_case2_PMF_Based_Group,FI_case2_AD_Based_Group)

plot(Tx_values, FI_case2_PMF_Based_Group, type = "l", col = "blue", lwd = 1, ylim = FI_Range,lty=1,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD and PMF Based FI: Imperfect Sensitivity under HG")
lines(Tx_values, FI_case2_AD_Based_Group, col = "blue", lwd = 1, lty = 2)

legend("topright",
       legend = c(
         "Perfect Case: PMF",
         "Perfect Case: AD",
         bquote(PMF: Imperfect ~ D[Delta] ),
         bquote(AD: Imperfect ~ D[Delta] ),
         bquote(PMF:Imperfect ~ D[Delta] ~ "," ~ G[Lambda]),
         bquote(AD:Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","black","blue", "blue", "red", "red"),
       lty = c(1,2,1,2,2,3),
       pch = c(19, 19, 19, 19, 19),
       bty = "n", # Remove box around the legend
       cex = 0.8) # Adjust text size

FI_Range <- range(FI_case4_PMF_Based_Group,FI_case4_AD_Based_Group)

plot(Tx_values, FI_case4_PMF_Based_Group, type = "l", col = "red", lwd = 1, ylim = FI_Range,lty=2,
     xlab = "Tx (Number of Contaminated Items)", ylab = "Fisher Information",
     main = "AD and PMF Based FI: Imperfect case under HG")
lines(Tx_values, FI_case4_AD_Based_Group, col = "red", lwd = 1, lty = 3)

legend("topright",
       legend = c(
         "Perfect Case: PMF",
         "Perfect Case: AD",
         bquote(PMF: Imperfect ~ D[Delta] ),
         bquote(AD: Imperfect ~ D[Delta] ),
         bquote(PMF:Imperfect ~ D[Delta] ~ "," ~ G[Lambda]),
         bquote(AD:Imperfect ~ D[Delta] ~ "," ~ G[Lambda])
       ),
       col = c("black","black","blue", "blue", "red", "red"),
       lty = c(1,2,1,2,2,3),
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
