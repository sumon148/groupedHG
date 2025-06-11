# Application to the study of estimating levels of Salmonella bacteria on pig farms

# Define constants
N <- 60  # Population size
barN <- c(20, 15, 10, 5, 1)  # Different group sizes
b <- rep(1,5)  # Number of groups
lambda <- 1  # Specificity for item-level
wp <- 20  # Total weight of the pooled sample
# wn <- b*ws*barN <- 20  # Total weight of the pooled sample
# wp <- ws*barN <- 20  # Total weight of the pooled sample
ws <- 20/barN  # Total weight of the pooled sample
rho_hat <- 0.55  # Empirical parameter for PCR
C <- 7.3  # Average number of clusters per gram
Tx_values <- seq(0, 60, by = 1)  # Range of contaminated items

# Sensitivity function (Equation 33) - Computes (1 - delta)
calculate_delta_ws <- function(barN, ws, rho_hat, C) {
  1 -  exp(-ws * C * (1 - exp(-rho_hat / (ws*barN))))
}


# Calculate item-level sensitivity (delta)
item_level_sensitivities_case2 <- calculate_delta_ws(barN, ws, rho_hat, C)


# Calculation of Hellinger Information

fisher_results_HG_PMF <- lapply(seq_along(barN), function(i) {

  # Extract corresponding values
  n <- barN[i]
  b_i <- b[i]  # Ensure correct mapping of b to barN
  ws_i <- ws[i]
  # Calculate delta based on barN
  delta <- calculate_delta_ws(barN=n, ws=ws_i, rho_hat=rho_hat, C=C)  # Ensure single numeric value
  lambda <- 1
  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfHG.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "PMF-HI", type = "item"
    )
  })

  # Store results for this barN
  list(barN = n, b=b_i,ws=ws_i,delta=delta,lambda=lambda,HI = fisher_info)
})



png("Paper Figures/Figure 2.png",width = 8,height = 6,units = "in", res = 300)


# Define colors for each barN
colors <- c("black", "blue", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_HG_PMF[[1]]$HI,
                  fisher_results_HG_PMF[[2]]$HI,
                  fisher_results_HG_PMF[[3]]$HI,
                  fisher_results_HG_PMF[[4]]$HI,
                  fisher_results_HG_PMF[[5]]$HI)

plot(Tx_values[c(2:60)], fisher_results_HG_PMF[[1]]$HI[c(2:60)], type = "l",
  xlab = expression(paste("Number of contaminated items ", T[x])), ylab = "Hellinger information (HI)",
  main = "HI: Grouped-hypergeometric model",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n")
axis(side = 1, at = seq(0, 60, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.005))


# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values[c(2:60)], fisher_results_HG_PMF[[i]]$HI[c(2:60)], type = "l", col = colors[i], lwd = 2, lty = i)
}

# Create legend labels
legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x * ~ ", " ~ bar(N) == y * ", " ~ delta == z, list(x = b[i], y =  barN[i], z =  round(item_level_sensitivities_case2[i],2)))
})

# Add legend to the plot
legend(
  20,0.04, legend = legend_labels,
  col = colors, lty = 1:length(barN), cex = 1, bty = "n"
)

# Add vertical lines
abline(v = c(9, 15, 51), col = c(2, 2, 2), lwd = 2)

# Define labels corresponding to each line
line_labels <- c(expression(T[x] == 0.15 * N ~ " = 9"),
                 expression(T[x] == 0.25 * N ~ " = 15"),
                 expression(T[x] == 0.85 * N ~ " = 51"))

# Add labels near the lines
text(x = c(9, 15, 51), y = par("usr")[4]*0.5, labels = line_labels, pos = 3, col = 2, cex = 0.7,srt = 90, # Rotation angle
     adj = 0)  # Align to the left of the line


dev.off()
