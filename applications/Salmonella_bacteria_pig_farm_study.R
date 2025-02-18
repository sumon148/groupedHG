# Application to the study of estimating levels of Salmonella bacteria on pig farms

# Application to Levels of Salmonella bacteria on pig farms: Case 1  --------------


# Define constants
N <- 80  # Population size
barN <- c(20, 15, 10, 5, 1)  # Different group sizes
b <- 1  # Number of groups
lambda <- 1  # Specificity for item-level
wp <- 20  # Total weight of the pooled sample
rho_hat <- 0.55  # Empirical parameter for PCR
C <- 7.3  # Average number of clusters per gram
Tx_values <- seq(0, 70, by = 1)  # Range of contaminated items

# Sensitivity function (Equation 33) - Computes (1 - delta)
calculate_delta <- function(barN, wp, rho_hat, C) {
  1 -  exp(-wp * C / barN * (1 - exp(-rho_hat / wp)))
}


calculate_delta_alt <- function(barN, wp, rho_hat, C) {
  1 -  (exp(-wp * C * (1 - exp(-rho_hat / wp))))^{1/barN}
}

# Calculate item-level sensitivity (delta)
item_level_sensitivities <- calculate_delta(barN, wp, rho_hat, C)
item_level_sensitivities_alt <- calculate_delta_alt(barN, wp, rho_hat, C)



# Generate FI for varying group sizes : AD and PMF based using HG Distribution ----------

item_level_sensitivities
# 0.1796423 0.2320416 0.3270132 0.5470888 0.9809425

# Different cases: Item level sensitivity

fisher_results_HG_AD <- lapply(barN, function(n) {
  # Calculate delta based on barN
  delta <- calculate_delta(n, wp, rho_hat, C)
  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfHG.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b,
      delta = delta, lambda = lambda,
      method = "AD", type = "item"
    )
  }
  )

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})


fisher_results_HG_PMF <- lapply(barN, function(n) {
  # Calculate delta based on barN
  delta <- calculate_delta(n, wp, rho_hat, C)

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfHG.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b,
      delta = delta, lambda = lambda,
      method = "PMF", type = "item"
    )
  }
  )

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})






png("Figure/Fig 3 HG AD Case 1.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "red", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_HG_AD[[1]]$fisher_info,
                  fisher_results_HG_AD[[2]]$fisher_info,
                  fisher_results_HG_AD[[3]]$fisher_info,
                  fisher_results_HG_AD[[4]]$fisher_info,
                  fisher_results_HG_AD[[5]]$fisher_info)
plot(
  Tx_values, fisher_results_HG_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.05))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- range(fisher_results_HG_AD[[1]]$fisher_info,
                  fisher_results_HG_AD[[2]]$fisher_info,
                  fisher_results_HG_AD[[3]]$fisher_info,
                  fisher_results_HG_AD[[4]]$fisher_info,
                  fisher_results_HG_AD[[5]]$fisher_info)
range_FI <- c(range_FI[1],0.05)


plot(
  Tx_values, fisher_results_HG_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.05, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



dev.off()






png("Figure/Fig 3 HG PMF Case 1.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "red", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_HG_PMF[[1]]$fisher_info,
                  fisher_results_HG_PMF[[2]]$fisher_info,
                  fisher_results_HG_PMF[[3]]$fisher_info,
                  fisher_results_HG_PMF[[4]]$fisher_info,
                  fisher_results_HG_PMF[[5]]$fisher_info)
plot(
  Tx_values, fisher_results_HG_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- range(fisher_results_HG_PMF[[1]]$fisher_info,
                  fisher_results_HG_PMF[[2]]$fisher_info,
                  fisher_results_HG_PMF[[3]]$fisher_info,
                  fisher_results_HG_PMF[[4]]$fisher_info,
                  fisher_results_HG_PMF[[5]]$fisher_info)
range_FI <- c(range_FI[1],0.01)

plot(
  Tx_values, fisher_results_HG_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.05, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)




dev.off()



# Generate FI for varying group sizes : AD and PMF based using BN Distribution ----------

item_level_sensitivities
# 0.1796423 0.2320416 0.3270132 0.5470888 0.9809425

# Different cases: Item level sensitivity

fisher_results_BN_AD <- lapply(barN, function(n) {
  # Calculate delta based on barN
  delta <- calculate_delta(n, wp, rho_hat, C)
  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfBN.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b,
      delta = delta, lambda = lambda,
      method = "AD", type = "item"
    )
  }
  )

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})


fisher_results_BN_PMF <- lapply(barN, function(n) {
  # Calculate delta based on barN
  delta <- calculate_delta(n, wp, rho_hat, C)

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfBN.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b,
      delta = delta, lambda = lambda,
      method = "PMF", type = "item"
    )
  }
  )

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})






png("Figure/Fig 3 BN AD Case 1.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "red", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_BN_AD[[1]]$fisher_info,
                  fisher_results_BN_AD[[2]]$fisher_info,
                  fisher_results_BN_AD[[3]]$fisher_info,
                  fisher_results_BN_AD[[4]]$fisher_info,
                  fisher_results_BN_AD[[5]]$fisher_info)
plot(
  Tx_values, fisher_results_BN_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Binomial",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.0025))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- range(fisher_results_BN_AD[[1]]$fisher_info,
                  fisher_results_BN_AD[[2]]$fisher_info,
                  fisher_results_BN_AD[[3]]$fisher_info,
                  fisher_results_BN_AD[[4]]$fisher_info,
                  fisher_results_BN_AD[[5]]$fisher_info)
range_FI <- c(range_FI[1],0.01)

plot(
  Tx_values, fisher_results_BN_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.01, by = 0.001))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)




dev.off()






png("Figure/Fig 3 BN PMF Case 1.png",width = 12,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "red", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_BN_PMF[[1]]$fisher_info,
                  fisher_results_BN_PMF[[2]]$fisher_info,
                  fisher_results_BN_PMF[[3]]$fisher_info,
                  fisher_results_BN_PMF[[4]]$fisher_info,
                  fisher_results_BN_PMF[[5]]$fisher_info)
plot(
  Tx_values, fisher_results_BN_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.0025))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- range(fisher_results_BN_PMF[[1]]$fisher_info,
                  fisher_results_BN_PMF[[2]]$fisher_info,
                  fisher_results_BN_PMF[[3]]$fisher_info,
                  fisher_results_BN_PMF[[4]]$fisher_info,
                  fisher_results_BN_PMF[[5]]$fisher_info)
range_FI <- c(range_FI[1],0.01)

plot(
  Tx_values, fisher_results_BN_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.01, by = 0.001))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(barN, function(n) {
  substitute(bar(N) == x, list(x = n)) # Creates \bar{N} = n
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)




dev.off()




# Application to Levels of Salmonella bacteria on pig farms: Case 2  --------------

# Define constants
N <- 80  # Population size
barN <- c(20, 10, 5, 1)  # Different group sizes
b <- c(1,2,4,20)  # Number of groups
lambda <- 1  # Specificity for item-level
wp <- 20  # Total weight of the pooled sample
# wn <- b*ws*barN <- 20  # Total weight of the pooled sample
# wp <- ws*barN <- 20  # Total weight of the pooled sample
ws <- 1  # Total weight of the pooled sample
rho_hat <- 0.55  # Empirical parameter for PCR
C <- 7.3  # Average number of clusters per gram
Tx_values <- seq(0, 70, by = 1)  # Range of contaminated items

# Sensitivity function (Equation 33) - Computes (1 - delta)
calculate_delta_ws <- function(barN, ws, rho_hat, C) {
  1 -  exp(-ws * C * (1 - exp(-rho_hat / (ws*barN))))
}


# Calculate item-level sensitivity (delta)
item_level_sensitivities_case2 <- calculate_delta_ws(barN, ws, rho_hat, C)

# Generate FI for varying group sizes : AD and PMF based using HG Distribution -------------

item_level_sensitivities_case2
# 0.1796423 0.3233888 0.5325257 0.9544191

fisher_results_HG_AD <- lapply(seq_along(barN), function(i) {

  # Extract corresponding values
  n <- barN[i]
  b_i <- b[i]  # Ensure correct mapping of b to barN

  # Calculate delta based on barN
  delta <- calculate_delta_ws(n, ws, rho_hat, C)  # Ensure single numeric value

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfHG.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "AD", type = "item"
    )
  })

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})



fisher_results_HG_PMF <- lapply(seq_along(barN), function(i) {

  # Extract corresponding values
  n <- barN[i]
  b_i <- b[i]  # Ensure correct mapping of b to barN

  # Calculate delta based on barN
  delta <- calculate_delta_ws(n, ws, rho_hat, C)  # Ensure single numeric value

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfHG.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "PMF", type = "item"
    )
  })

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})


png("Figure/Fig 3 HG AD Case 2.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_HG_AD[[1]]$fisher_info,
                  fisher_results_HG_AD[[2]]$fisher_info,
                  fisher_results_HG_AD[[3]]$fisher_info,
                  fisher_results_HG_AD[[4]]$fisher_info)

plot(
  Tx_values, fisher_results_HG_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.05))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Create legend labels
legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- c(range_FI[1],0.05)

plot(
  Tx_values, fisher_results_HG_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.05, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



dev.off()



png("Figure/Fig 3 HG PMF Case 2.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_HG_PMF[[1]]$fisher_info,
                  fisher_results_HG_PMF[[2]]$fisher_info,
                  fisher_results_HG_PMF[[3]]$fisher_info,
                  fisher_results_HG_PMF[[4]]$fisher_info)

#range_FI <- c(0,1)
plot(
  Tx_values, fisher_results_HG_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.05))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}


legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



# Plot results
range_FI <- c(range_FI[1],0.05)

plot(
  Tx_values, fisher_results_HG_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.03, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_HG_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}


legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)



dev.off()





# ----------------------------------------------------#
# Generate FI for varying group sizes : AD and PMF based using BN Distribution -------------

item_level_sensitivities_case2
# 0.1796423 0.3233888 0.5325257 0.9544191

fisher_results_BN_AD <- lapply(seq_along(barN), function(i) {

  # Extract corresponding values
  n <- barN[i]
  b_i <- b[i]  # Ensure correct mapping of b to barN

  # Calculate delta based on barN
  delta <- calculate_delta_ws(n, ws, rho_hat, C)  # Ensure single numeric value

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfBN.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "AD", type = "item"
    )
  })

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})



fisher_results_BN_PMF <- lapply(seq_along(barN), function(i) {

  # Extract corresponding values
  n <- barN[i]
  b_i <- b[i]  # Ensure correct mapping of b to barN

  # Calculate delta based on barN
  delta <- calculate_delta_ws(n, ws, rho_hat, C)  # Ensure single numeric value

  # Compute Fisher Information for each Tx value
  fisher_info <- sapply(Tx_values, function(Tx) {
    FIpmfBN.Tx.imperfect(
      N = N, barN = n, Tx = Tx, b = b_i,  # Use correct `b_i`
      delta = delta, lambda = lambda,
      method = "PMF", type = "item"
    )
  })

  # Store results for this barN
  list(barN = n, fisher_info = fisher_info)
})


png("Figure/Fig 3 BN AD Case 2.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_BN_AD[[1]]$fisher_info,
                  fisher_results_BN_AD[[2]]$fisher_info,
                  fisher_results_BN_AD[[3]]$fisher_info,
                  fisher_results_BN_AD[[4]]$fisher_info)

plot(
  Tx_values, fisher_results_BN_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.025))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Create legend labels
legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- c(range_FI[1],0.05)

plot(
  Tx_values, fisher_results_BN_AD[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (AD) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.05, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_AD[[i]]$fisher_info, col = colors[i], cex = 0.5)
}

# Create legend labels
legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)

dev.off()



png("Figure/Fig 3 BN PMF Case 2.png",width = 10,height = 6,units = "in", res = 300)

par(mfrow=c(1,2))

# Define colors for each barN
colors <- c("blue", "green", "purple", "orange")  # Unique colors for each group size


# Plot results
range_FI <- range(fisher_results_BN_PMF[[1]]$fisher_info,
                  fisher_results_BN_PMF[[2]]$fisher_info,
                  fisher_results_BN_PMF[[3]]$fisher_info,
                  fisher_results_BN_PMF[[4]]$fisher_info)

#range_FI <- c(0,1)
plot(
  Tx_values, fisher_results_BN_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, range_FI[2], by = 0.01))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}


legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


# Plot results
range_FI <- c(range_FI[1],0.05)

plot(
  Tx_values, fisher_results_BN_PMF[[1]]$fisher_info, type = "l",
  xlab = expression(T[X]), ylab = "Fisher Information: Hypergeometric",
  col = colors[1], lwd = 2,ylim=range_FI,xaxt = "n",yaxt="n",
  main = expression("Item level (PMF) FI for varying group" ~ bar(N)) # \bar{N}
)
axis(side = 1, at = seq(0, 70, by = 5))
axis(side = 2, at = seq(0, 0.05, by = 0.005))

# Add lines for each barN with distinct colors
for (i in 2:length(barN)) {
  lines(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, type = "l", col = colors[i], lwd = 2, lty = i)
}


# Add points for each barN with matching colors
for (i in 1:length(barN)) {
  points(Tx_values, fisher_results_BN_PMF[[i]]$fisher_info, col = colors[i], cex = 0.5)
}


legend_labels <- lapply(seq_along(barN), function(i) {
  substitute(b == x ~ ", " ~ bar(N) == y, list(x = b[i], y =  barN[i]))
})

# Add legend to the plot
legend(
  "topright", legend = legend_labels,
  col = colors, lty = 1:length(barN), pch = rep(19, length(barN)), cex = 1, bty = "n"
)


dev.off()



# Test Cases by Mahadi : -------------------

# Case 1:##################################################
b = 1;
barN = 20;
N = 80;
delta = 0.1796423;
lambda = 1;
Tx.values = c(0,1, 3,5,10,15,20)
ty.values <- c(0:b)

# Compute the PMF for all (ty, Tx) pairs
PRtydeltalambda.HG <- outer(ty.values, Tx.values, Vectorize(function(ty, Tx) {
  pmfHG.imperfect.item(ty, N, barN, Tx, b, delta=delta, lambda=lambda)
}))

# Convert to a data frame for readability
colnames(PRtydeltalambda.HG) <- paste0("Tx=", Tx.values)
rownames(PRtydeltalambda.HG) <- paste0("ty=", ty.values)
PRtydeltalambda.HG
colSums(PRtydeltalambda.HG)

# Case 2:##################################################
b = 1;
barN = 20;
N = b*barN;
delta = 0.1796423;
lambda = 1;
Tx.values = c(0,1, 3,5,10,15,20)
ty.values <- c(0:b)

# Compute the PMF for all (ty, Tx) pairs
PRtydeltalambda.HG <- outer(ty.values, Tx.values, Vectorize(function(ty, Tx) {
  pmfHG.imperfect.item(ty, N, barN, Tx, b, delta=delta, lambda=lambda)
}))

# Convert to a data frame for readability
colnames(PRtydeltalambda.HG) <- paste0("Tx=", Tx.values)
rownames(PRtydeltalambda.HG) <- paste0("ty=", ty.values)
PRtydeltalambda.HG
colSums(PRtydeltalambda.HG)


# Case 3:##################################################
b = 5;
barN = 20;
N = b*barN;
delta = 0.1796423;
lambda = 0.9;
Tx.values = c(0,1, 3,5,10,15,20)
ty.values <- c(0:b)

# Compute the PMF for all (ty, Tx) pairs
PRtydeltalambda.HG <- outer(ty.values, Tx.values, Vectorize(function(ty, Tx) {
  pmfHG.imperfect.item(ty, N, barN, Tx, b, delta=delta, lambda=lambda)
}))

# Convert to a data frame for readability
colnames(PRtydeltalambda.HG) <- paste0("Tx=", Tx.values)
rownames(PRtydeltalambda.HG) <- paste0("ty=", ty.values)
PRtydeltalambda.HG
colSums(PRtydeltalambda.HG)



# Case 4:##################################################
b = 5;
barN = 20;
N = b*barN;
delta = 0.1796423;
lambda = 0.7;
Tx.values = c(0,1, 3,5,10,15,20)
ty.values <- c(0:b)

# Compute the PMF for all (ty, Tx) pairs
PRtydeltalambda.HG <- outer(ty.values, Tx.values, Vectorize(function(ty, Tx) {
  pmfHG.imperfect.item(ty, N, barN, Tx, b, delta=delta, lambda=lambda)
}))

# Convert to a data frame for readability
colnames(PRtydeltalambda.HG) <- paste0("Tx=", Tx.values)
rownames(PRtydeltalambda.HG) <- paste0("ty=", ty.values)
PRtydeltalambda.HG
colSums(PRtydeltalambda.HG)


# Case 5:##################################################
b = 1;
barN = 20;
N = 80;
delta = 1;
lambda = 1;
Tx.values = c(0,1, 3,5,10,15,20)
ty.values <- c(0:b)

# Compute the PMF for all (ty, Tx) pairs
PRtydeltalambda.HG <- outer(ty.values, Tx.values, Vectorize(function(ty, Tx) {
  pmfHG.imperfect.item(ty, N, barN, Tx, b, delta=delta, lambda=lambda)
}))

# Convert to a data frame for readability
colnames(PRtydeltalambda.HG) <- paste0("Tx=", Tx.values)
rownames(PRtydeltalambda.HG) <- paste0("ty=", ty.values)
PRtydeltalambda.HG
colSums(PRtydeltalambda.HG)

# Case 6:##################################################
b = 1;
barN = 20;
N = b*barN;
delta = 1;
lambda = 0.7;
Tx.values = c(0,1, 3,5,10,15,20)
ty.values <- c(0:b)

# Compute the PMF for all (ty, Tx) pairs
PRtydeltalambda.HG <- outer(ty.values, Tx.values, Vectorize(function(ty, Tx) {
  pmfHG.imperfect.item(ty, N, barN, Tx, b, delta=delta, lambda=lambda)
}))

# Convert to a data frame for readability
colnames(PRtydeltalambda.HG) <- paste0("Tx=", Tx.values)
rownames(PRtydeltalambda.HG) <- paste0("ty=", ty.values)
PRtydeltalambda.HG
colSums(PRtydeltalambda.HG)






# Check sum of t+y by Tx: Case 1 ----------------

# Define constants
N <- 80  # Population size
barN <- c(20, 15, 10, 5, 1)  # Different group sizes
b <- 1  # Number of groups
lambda <- 1  # Specificity for item-level
wp <- 20  # Total weight of the pooled sample
rho_hat <- 0.55  # Empirical parameter for PCR
C <- 7.3  # Average number of clusters per gram
Tx_values <- seq(0, 70, by = 1)  # Range of contaminated items

# Sensitivity function (Equation 33) - Computes (1 - delta)
calculate_delta <- function(barN, wp, rho_hat, C) {
  1 -  exp(-wp * C / barN * (1 - exp(-rho_hat / wp)))
}


calculate_delta_alt <- function(barN, wp, rho_hat, C) {
  1 -  (exp(-wp * C * (1 - exp(-rho_hat / wp))))^{1/barN}
}

# Calculate item-level sensitivity (delta)
item_level_sensitivities <- calculate_delta(barN, wp, rho_hat, C)

FIpmfHG.Tx.imperfect.detail(N=80,barN=20,Tx=0,b=1,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=20,Tx=1,b=1,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=20,Tx=2,b=1,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=20,Tx=3,b=1,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=20,Tx=4,b=1,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)

Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_salmonella_case1_AD <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N=80,barN=20,b=1,delta=item_level_sensitivities[1],lambda = 1,method = "AD",type = "item",verbose = FALSE)

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_salmonella_case1_AD <- bind_rows(FI_results_salmonella_case1_AD)

FI_results_salmonella_case1_AD_delta_1 <- data.frame(cbind(Tx=FI_results_salmonella_case1_AD$Tx,ty=FI_results_salmonella_case1_AD$ty,p_ty=FI_results_salmonella_case1_AD$p_ty,
                                        dp_ty=FI_results_salmonella_case1_AD$dp_ty,FI_Tx=FI_results_salmonella_case1_AD$FI_Tx))
# View(FI_results_salmonella_case1_AD_delta_1)

tapply(FI_results_salmonella_case1_AD_delta_1$p_ty,FI_results_salmonella_case1_AD_delta_1$Tx,sum)


Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_salmonella_case1_AD <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N=80,barN=20,b=1,delta=item_level_sensitivities[2],lambda = 1,method = "AD",type = "item",verbose = FALSE)

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_salmonella_case1_AD <- bind_rows(FI_results_salmonella_case1_AD)

FI_results_salmonella_case1_AD_delta_2 <- data.frame(cbind(Tx=FI_results_salmonella_case1_AD$Tx,ty=FI_results_salmonella_case1_AD$ty,p_ty=FI_results_salmonella_case1_AD$p_ty,
                                                dp_ty=FI_results_salmonella_case1_AD$dp_ty,FI_Tx=FI_results_salmonella_case1_AD$FI_Tx))
# View(FI_results_salmonella_case1_AD_delta_2)

tapply(FI_results_salmonella_case1_AD_delta_2$p_ty,FI_results_salmonella_case1_AD_delta_2$Tx,sum)



Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_salmonella_case1_AD <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N=80,barN=20,b=1,delta=item_level_sensitivities[3],lambda = 1,method = "AD",type = "item",verbose = FALSE)

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_salmonella_case1_AD <- bind_rows(FI_results_salmonella_case1_AD)

FI_results_salmonella_case1_AD_delta_3 <- data.frame(cbind(Tx=FI_results_salmonella_case1_AD$Tx,ty=FI_results_salmonella_case1_AD$ty,p_ty=FI_results_salmonella_case1_AD$p_ty,
                                                           dp_ty=FI_results_salmonella_case1_AD$dp_ty,FI_Tx=FI_results_salmonella_case1_AD$FI_Tx))
# View(FI_results_salmonella_case1_AD_delta_3)

tapply(FI_results_salmonella_case1_AD_delta_3$p_ty,FI_results_salmonella_case1_AD_delta_3$Tx,sum)


Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_salmonella_case1_AD <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N=80,barN=20,b=1,delta=item_level_sensitivities[4],lambda = 1,method = "AD",type = "item",verbose = FALSE)

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_salmonella_case1_AD <- bind_rows(FI_results_salmonella_case1_AD)

FI_results_salmonella_case1_AD_delta_4 <- data.frame(cbind(Tx=FI_results_salmonella_case1_AD$Tx,ty=FI_results_salmonella_case1_AD$ty,p_ty=FI_results_salmonella_case1_AD$p_ty,
                                                           dp_ty=FI_results_salmonella_case1_AD$dp_ty,FI_Tx=FI_results_salmonella_case1_AD$FI_Tx))
# View(FI_results_salmonella_case1_AD_delta_4)

tapply(FI_results_salmonella_case1_AD_delta_4$p_ty,FI_results_salmonella_case1_AD_delta_4$Tx,sum)



Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_salmonella_case1_AD <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N=80,barN=20,b=1,delta=item_level_sensitivities[5],lambda = 1,method = "AD",type = "item",verbose = FALSE)

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_salmonella_case1_AD <- bind_rows(FI_results_salmonella_case1_AD)

FI_results_salmonella_case1_AD_delta_5 <- data.frame(cbind(Tx=FI_results_salmonella_case1_AD$Tx,ty=FI_results_salmonella_case1_AD$ty,p_ty=FI_results_salmonella_case1_AD$p_ty,
                                                           dp_ty=FI_results_salmonella_case1_AD$dp_ty,FI_Tx=FI_results_salmonella_case1_AD$FI_Tx))
# View(FI_results_salmonella_case1_AD_delta_5)

tapply(FI_results_salmonella_case1_AD_delta_5$p_ty,FI_results_salmonella_case1_AD_delta_5$Tx,sum)

plot(FI_results_salmonella_case1_AD_delta_5$Tx,FI_results_salmonella_case1_AD_delta_5$FI_Tx)
plot(FI_results_salmonella_case1_AD_delta_5$Tx,FI_results_salmonella_case1_AD_delta_4$FI_Tx)
plot(FI_results_salmonella_case1_AD_delta_5$Tx,FI_results_salmonella_case1_AD_delta_3$FI_Tx)
plot(FI_results_salmonella_case1_AD_delta_5$Tx,FI_results_salmonella_case1_AD_delta_2$FI_Tx)
plot(FI_results_salmonella_case1_AD_delta_5$Tx,FI_results_salmonella_case1_AD_delta_1$FI_Tx)


# Check sum of ty by Tx: Case 2 ----------------

# Define constants
N <- 80  # Population size
barN <- c(20, 10, 5, 1)  # Different group sizes
b <- c(1,2,4,20)  # Number of groups
lambda <- 1  # Specificity for item-level
wp <- 20  # Total weight of the pooled sample
# wn <- b*ws*barN <- 20  # Total weight of the pooled sample
# wp <- ws*barN <- 20  # Total weight of the pooled sample
ws <- 1  # Total weight of the pooled sample
rho_hat <- 0.55  # Empirical parameter for PCR
C <- 7.3  # Average number of clusters per gram
Tx_values <- seq(0, 70, by = 1)  # Range of contaminated items

# Sensitivity function (Equation 33) - Computes (1 - delta)
calculate_delta_ws <- function(barN, ws, rho_hat, C) {
  1 -  exp(-ws * C * (1 - exp(-rho_hat / (ws*barN))))
}


# Calculate item-level sensitivity (delta)
item_level_sensitivities_case2 <- calculate_delta_ws(barN, ws, rho_hat, C)

FIpmfHG.Tx.imperfect.detail(N=80,barN=20,Tx=0,b=1,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=10,Tx=0,b=2,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=5,Tx=0,b=4,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=1,Tx=0,b=20,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)

FIpmfHG.Tx.imperfect.detail(N=80,barN=20,Tx=1,b=1,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=10,Tx=1,b=2,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=5,Tx=1,b=4,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=1,Tx=1,b=20,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)

FIpmfHG.Tx.imperfect.detail(N=80,barN=20,Tx=20,b=1,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=10,Tx=20,b=2,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=5,Tx=20,b=4,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)
FIpmfHG.Tx.imperfect.detail(N=80,barN=1,Tx=20,b=20,delta=0.1796423,lambda = 1,method = "AD",type = "item",verbose = FALSE)


Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_salmonella_case2_AD <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N=80,barN=barN[1],b=b[1],delta=item_level_sensitivities_case2[1],lambda = 1,method = "AD",type = "item",verbose = FALSE)

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_salmonella_case2_AD <- bind_rows(FI_results_salmonella_case2_AD)

FI_results_salmonella_case2_AD_delta_1 <- data.frame(cbind(Tx=FI_results_salmonella_case2_AD$Tx,ty=FI_results_salmonella_case2_AD$ty,p_ty=FI_results_salmonella_case2_AD$p_ty,
                                                           dp_ty=FI_results_salmonella_case2_AD$dp_ty,FI_Tx=FI_results_salmonella_case2_AD$FI_Tx))
# View(FI_results_salmonella_case2_AD_delta_1)

tapply(FI_results_salmonella_case2_AD_delta_1$p_ty,FI_results_salmonella_case2_AD_delta_1$Tx,sum)



Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_salmonella_case2_AD <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N=80,barN=barN[2],b=b[2],delta=item_level_sensitivities_case2[2],lambda = 1,method = "AD",type = "item",verbose = FALSE)

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_salmonella_case2_AD <- bind_rows(FI_results_salmonella_case2_AD)

FI_results_salmonella_case2_AD_delta_2 <- data.frame(cbind(Tx=FI_results_salmonella_case2_AD$Tx,ty=FI_results_salmonella_case2_AD$ty,p_ty=FI_results_salmonella_case2_AD$p_ty,
                                                           dp_ty=FI_results_salmonella_case2_AD$dp_ty,FI_Tx=FI_results_salmonella_case2_AD$FI_Tx))
# View(FI_results_salmonella_case2_AD_delta_2)

tapply(FI_results_salmonella_case2_AD_delta_2$p_ty,FI_results_salmonella_case2_AD_delta_2$Tx,sum)



Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_salmonella_case2_AD <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N=80,barN=barN[3],b=b[3],delta=item_level_sensitivities_case2[3],lambda = 1,method = "AD",type = "item",verbose = FALSE)

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_salmonella_case2_AD <- bind_rows(FI_results_salmonella_case2_AD)

FI_results_salmonella_case2_AD_delta_3 <- data.frame(cbind(Tx=FI_results_salmonella_case2_AD$Tx,ty=FI_results_salmonella_case2_AD$ty,p_ty=FI_results_salmonella_case2_AD$p_ty,
                                                           dp_ty=FI_results_salmonella_case2_AD$dp_ty,FI_Tx=FI_results_salmonella_case2_AD$FI_Tx))
# View(FI_results_salmonella_case2_AD_delta_3)

tapply(FI_results_salmonella_case2_AD_delta_3$p_ty,FI_results_salmonella_case2_AD_delta_3$Tx,sum)



Tx_values <- c(0:80)
# Apply function over all Tx values
FI_results_salmonella_case2_AD <- lapply(Tx_values, function(Tx) {
  result <- FIpmfHG.Tx.imperfect.detail(Tx, N=80,barN=barN[4],b=b[4],delta=item_level_sensitivities_case2[4],lambda = 1,method = "AD",type = "item",verbose = FALSE)

  # Extract FI_Tx and ty.data, and add Tx column
  data.frame(Tx = Tx, FI_Tx = result$FI_Tx, result$ty.data)
})

# Combine into a single dataframe
FI_results_salmonella_case2_AD <- bind_rows(FI_results_salmonella_case2_AD)

FI_results_salmonella_case2_AD_delta_4 <- data.frame(cbind(Tx=FI_results_salmonella_case2_AD$Tx,ty=FI_results_salmonella_case2_AD$ty,p_ty=FI_results_salmonella_case2_AD$p_ty,
                                                           dp_ty=FI_results_salmonella_case2_AD$dp_ty,FI_Tx=FI_results_salmonella_case2_AD$FI_Tx))
# View(FI_results_salmonella_case2_AD_delta_4)

tapply(FI_results_salmonella_case2_AD_delta_4$p_ty,FI_results_salmonella_case2_AD_delta_4$Tx,sum)
