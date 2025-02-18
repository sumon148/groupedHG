#' Factorial Custom Function
#'
#' Helper function to calculate the factorial of a non-negative integer.
#' @param n A non-negative integer.
#' @return The factorial of \code{n}.
#' @examples
#' factorial_custom(5) # Returns 120
#' factorial_custom(0) # Returns 1
#' @export
factorial_custom <- function(n) {
  if (n == 0) return(1)
  return(prod(1:n))
}



#' Binomial Coefficient
#'
#' Function to calculate the binomial coefficient (\code{n choose k}).
#' @param n The number of items.
#' @param k The number of items to choose.
#' @return The binomial coefficient (\code{n choose k}). Returns 0 if \code{k > n} or \code{k < 0}.
#' @examples
#' choose_custom(5, 2) # Returns 10
#' choose_custom(5, 6) # Returns 0
#' @export
choose_custom <- function(n, k) {
  if (k > n || k < 0) return(0)
  return(factorial_custom(n) / (factorial_custom(k) * factorial_custom(n - k)))
}


#' Binomial Coefficient: log-gamma  form
#'
#' Function to calculate the binomial coefficient (\code{n choose k}) in log-gamma  form.
#' @param n The number of items.
#' @param k The number of items to choose.
#' @return The binomial coefficient (\code{n choose k}). Returns 0 if \code{k > n} or \code{k < 0}.
#' @examples
#' choose_custom_log(5, 2) # Returns 10
#' choose_custom_log(5, 6) # Returns 0
#' @export
choose_custom_log <- function(n, k) {
  if (k > n || k < 0) return(0)  # Handle invalid cases
  exp(lgamma(n + 1) - (lgamma(k + 1) + lgamma(n - k + 1)))
}



#' Safe Digamma Function
#'
#' A safe version of the digamma function that returns 0 for invalid arguments (\code{x <= 0}).
#' @param x A numeric value.
#' @return The digamma of \code{x} if \code{x > 0}, otherwise 0.
#' @examples
#' safe_digamma(5) # Returns digamma(5)
#' safe_digamma(-1) # Returns 0
#' @export
safe_digamma <- function(x) {
  if (x <= 0) {
    return(0)  # Return 0 for invalid arguments
  }
  return(digamma(x))
}

#' PMF for Grouped Hypergeometric Distribution under Perfect Test
#'
#' Function to calculate the PMF for a grouped hypergeometric distribution under a perfect test.
#'
#' @param ty Number of groups where the target outcome is observed.
#' @param N Total population size.
#' @param barN Average number of individuals per group.
#' @param Tx Number of target individuals in the population.
#' @param b Total number of groups.
#' @return The probability mass function (PMF) value for the given parameters.
#' @examples
#' pmfHG.perfect(ty = 2, N = 100, barN = 10, Tx = 20, b = 5)
#' @export
pmfHG.perfect <- function(ty, N, barN, Tx, b) {
  # Initialize the outer binomial coefficient
  outer_binomial <- choose_custom_log(b, ty)  # Binomial coefficient for selecting t_y groups

  # Initialize the summation for inclusion-exclusion
  sum_inclusion_exclusion <- 0

  # Loop through each subset of groups to apply inclusion-exclusion
  for (j in 0:ty) {
    # Calculate the binomial coefficient for selecting j groups from ty
    inner_binomial <- choose_custom_log(ty, j)

    # Calculate the remaining population size after excluding certain groups
    exclusion_population <- N - barN * (b - j)

    # Check if the exclusion population is large enough to support Tx


    if (exclusion_population < 0) {
      next  # Skip this iteration as the exclusion population is invalid
    } else if (exclusion_population == Tx && Tx == 0) {
      hypergeo_prob <- 1  # Set to 1 when both are zero
    } else if (exclusion_population >= Tx) {
      # Compute the hypergeometric probability
      hypergeo_prob <- choose_custom_log(exclusion_population, Tx) / choose_custom_log(N, Tx)
    }
      # Calculate the inclusion-exclusion term
      term <- (-1)^(ty - j) * inner_binomial * hypergeo_prob

      # Add the term to the summation
      sum_inclusion_exclusion <- sum_inclusion_exclusion + term
    }


  # Compute the final PMF value
  pmf <- outer_binomial * sum_inclusion_exclusion

  # Apply a threshold for small probabilities to prevent numerical instability
  if (abs(pmf) < 1e-5) pmf <- 0

  # Return the computed PMF value
  return(pmf)
}


#' Derivative of PMF w.r.t Tx (Perfect Case)
#'
#' Computes the derivative of the PMF for the perfect case as per Equation 26.
#'
#' @param ty Integer. Number of groups with observed target outcomes.
#' @param N Integer. Total population size.
#' @param barN Numeric. Average group size.
#' @param Tx Integer. Number of target individuals in the population.
#' @param b Integer. Total number of groups.
#' @param verbose Logical. If TRUE, prints conditions that are not fulfilled.
#' @return Numeric. The derivative of the PMF for the perfect case.
#' @export
dpmfHG.perfect <- function(ty, N, barN, Tx, b, verbose=TRUE) {
  # Calculate the PMF for the given inputs
  p_ty_given_Tx <- pmfHG.perfect(ty, N, barN, Tx, b)

  if (ty > b || ty > Tx && p_ty_given_Tx==0) {  # ty must be <= min(b, Tx) in perfect case
    return(0)
  }

  # If the PMF is zero, the derivative is also zero
  if (p_ty_given_Tx > 0) {
    # Binomial coefficient for selecting t_y groups out of b
    outer_binomial <- choose_custom_log(b, ty)
    # Initialize the sum for the inclusion-exclusion principle
    sum_inclusion_exclusion <- 0

    # Iterate over subsets (j = 0 to ty) to compute inclusion-exclusion terms
    for (j in 0:ty) {
      # Binomial coefficient for selecting j groups out of ty
      inner_binomial <- choose_custom_log(ty, j)
      # Size of the population excluding j groups
      group_size <- barN * (b - j)
      exclusion_population <- N - group_size

      # Handle special cases
      if (exclusion_population < 0) {
        next  # Skip this iteration as the exclusion population is invalid
      } else if (exclusion_population == Tx && Tx == 0) {
        der_hypergeo_prob <- 1  # Set to 1 when both are zero
      } else if (exclusion_population >= Tx) {
        # Calculate the hypergeometric probability for the excluded population
        hypergeo_prob <- choose_custom_log(N-Tx, group_size) / choose_custom_log(N, group_size)
        # Compute the derivative term using the digamma function
        digamma_term <- safe_digamma(N - Tx - group_size + 1) - safe_digamma(N - Tx + 1) # paper
        # digamma_term <- safe_digamma(N - Tx + 1) - safe_digamma(N - Tx - group_size + 1) - safe_digamma(Tx + 1) # ChatGPT
        der_hypergeo_prob <- hypergeo_prob * digamma_term
      } else {
        if (verbose){
          # Debugging output for intermediate values (optional)
          cat("  j =", j,
              "| inner_binomial =", inner_binomial,
              "| hypergeo_prob =", hypergeo_prob,
              "| digamma_term =", digamma_term,
              "| term =", term, "\n")
        }
        next
      }

      # Inclusion-exclusion term with alternating signs
      term <- (-1)^(ty - j) * inner_binomial * der_hypergeo_prob

      # Accumulate the term in the sum
      sum_inclusion_exclusion <- sum_inclusion_exclusion + term
    }


    # Final derivative is the product of the outer binomial and the sum of terms
    der_pmf <- outer_binomial * sum_inclusion_exclusion
  } else {
    # If the PMF is zero, the derivative is set to zero
    der_pmf <- 0
  }

  # Return the computed derivative
  return(der_pmf)
}


#' Fisher Information for \eqn{T_x} (Perfect Case)
#'
#' Calculates Fisher Information for the perfect case using different methods.
#'
#' @param N Integer. Total population size.
#' @param barN Numeric. Average group size.
#' @param Tx Integer. Number of target individuals in the population.
#' @param b Integer. Total number of groups.
#' @param method Character. Method to calculate the derivative. Options are:
#'   \itemize{
#'     \item "AD" (Analytic Derivative)
#'     \item "ND" (Numerical Derivative)
#'     \item "PMF" (PMF-based approximation)
#'   }
#' @return Numeric. Fisher Information for the given \eqn{T_x}.
#' @export
FIpmfHG.Tx.perfect <- function(N, barN, Tx, b, method = c("AD", "ND", "PMF")) {

  FI_Tx <- 0  # Initialize Fisher Information

  # Ensure the method is valid
  method <- match.arg(method)

  for (ty in 0:b) {
    # Skip invalid cases where ty > Tx
    if (ty > b || ty > Tx) next

    # PMF for ty given Tx
    P_ty <- pmfHG.perfect(ty, N, barN, Tx, b)

    # Skip terms with P(t_y | T_x) = 0 to avoid division by zero
    if (P_ty == 0) next

    if (method == "AD") {
      # Method: Analytic Derivative (AD)
      dP_Tx <- dpmfHG.perfect(ty, N, barN, Tx, b)
      if (!is.na(dP_Tx) && dP_Tx != 0) {  # Check for non-NA and non-zero derivative
        FI_Tx <- FI_Tx + (dP_Tx^2 / P_ty)
      }

    } else if (method == "PMF") {
      # Method: PMF-based approximation (PMF)
      P_ty_Tx_plus1 <- pmfHG.perfect(ty, N, barN, Tx + 1, b)
      if (P_ty_Tx_plus1 > 0) {
        FI_Tx <- FI_Tx + 4 * (sqrt(P_ty_Tx_plus1) - sqrt(P_ty))^2
      }

    } else if (method == "ND") {
      # Method: Numerical Derivative (ND)
      h <- 1  # Step size for forward difference
      P_plus <- pmfHG.perfect(ty, N, barN, Tx + h, b)
      if (P_plus > 0) {
        FI_Tx <- FI_Tx + ((P_plus - P_ty) / h)^2 / P_ty
      }
    }
  }

  # Return the computed Fisher Information
  return(FI_Tx)
}


#' Compute \eqn{q_{{\Delta\Lambda} k}} for group-Level Sensitivity and Specificity
#'
#' Calculates \eqn{q_{{\Delta\Lambda} k}}, the probability that none are identified as contaminated for the grouped case.
#'
#' @param k Integer. Number of groups considered in the summation.
#' @param N Integer. Total population size.
#' @param Tx Integer. Number of successes (positive cases) to be distributed.
#' @param barN Numeric. Average group size.
#' @param delta Numeric. False negative rate (0 ≤ delta ≤ 1).
#' @param lambda Numeric. False positive rate (0 ≤ lambda ≤ 1).
#' @return Numeric. The probability \eqn{q_{\Delta\Lambda k}}.
#' @name qDLkGroup
#' @rdname qDLkGroup
#' @export
qDLkGroup <- function(k, N, Tx, barN, delta, lambda) {
  # Input validation
  if (N <= 0 || Tx < 0 || barN <= 0 || delta < 0 || delta > 1 || lambda < 0 || lambda > 1) {
    stop("Invalid input parameters.")
  }

  # Precompute constant values
  total_choose <- choose_custom_log(N, Tx)

  # Initialize q_Delta_Lambda_k to 0
  q_Delta_Lambda_k <- 0

  # Outer loop: iterate over ell from 0 to k
  for (ell in 0:k) {
    # Initialize the inner summation for the current value of ell
    inner_sum <- 0

    # Inner loop: iterate over i from 0 to ell
    for (i in 0:ell) {
      # Calculate the exclusion population size
      exclusion_population <- N - barN * (k - (ell - i))


    # Check if the exclusion population is valid for hypergeometric probability
      if (exclusion_population >= Tx) {
        # Compute the term for this combination of ell and i
        term <- (-1)^i * choose_custom_log(ell, i) *
          (choose_custom_log(exclusion_population, Tx) / total_choose)

        # Add the term to the inner summation
        inner_sum <- inner_sum + term
      }
    }

    # Add the contribution from the current ell to q_Delta_Lambda_k
    q_Delta_Lambda_k <- q_Delta_Lambda_k + choose_custom_log(k, ell) * inner_sum *
      (1 - delta)^ell * lambda^(k - ell)
  }

  # Return the computed value of q_Delta_Lambda_k
  return(q_Delta_Lambda_k)
}


#' Derivative of \eqn{q_{{\Delta\Lambda} k}} with respect to \eqn{T_x}
#'
#' Calculates the derivative of \eqn{q_{{\Delta \Lambda} k}} with respect to \eqn{T_x} for the grouped case.
#'
#' @param k Integer. Number of groups considered in the summation.
#' @param N Integer. Total population size.
#' @param Tx Integer. Number of successes (positive cases) to be distributed.
#' @param barN Numeric. Average group size.
#' @param delta Numeric. False negative rate (0 ≤ delta ≤ 1).
#' @param lambda Numeric. False positive rate (0 ≤ lambda ≤ 1).
#' @param verbose Logical. If TRUE, prints conditions that are not fulfilled.
#' @return Numeric. The derivative of \eqn{q_{\Delta\Lambda k}} with respect to \eqn{T_x}.
#' @name dqDLkGroup
#' @rdname dqDLkGroup
#' @export
dqDLkGroup <- function(k, N, Tx, barN, delta, lambda,verbose=TRUE) {
  # Step 1: Function Purpose
  # This function calculates the derivative of q_{Delta, Lambda, k} with respect to T_x for the grouped case.
  # The calculation follows Equation 26, considering a grouped population model with false negative and false positive rates.

  # Step 2: Input Parameters Description
  # k: Number of groups considered in the summation.
  # N: Total population size.
  # Tx: Number of successes (positive cases) to be distributed.
  # barN: Average size of each group.
  # delta: False negative rate (0 ≤ delta ≤ 1).
  # lambda: False positive rate (0 ≤ lambda ≤ 1).

  # Step 3: Input Validation
  # Validate input parameters to ensure they are within expected ranges.
  if (N <= 0 || Tx < 0 || barN <= 0 || delta < 0 || delta > 1 || lambda < 0 || lambda > 1) {
    stop("Invalid input parameters.")
  }

  # Step 4: Initialize Output
  # Initialize the derivative value as zero to accumulate results from the summations.
  dq_Delta_Lambda_k <- 0

  # Step 5: Outer Summation Over ell
  # Iterate over possible values of ell (from 0 to k).
  for (ell in 0:k) {
    inner_sum <- 0  # Initialize the inner summation for each ell.

    # Step 6: Inner Summation Over i
    # Iterate over possible values of i (from 0 to ell) for each ell.
    for (i in 0:ell) {
      # Step 7: Exclusion Population Calculation
      # Calculate the exclusion population based on current ell and i.
      exclusion_population <- N - barN * (k - (ell - i))

      # Step 8: Validity Check for Hypergeometric Term
      # Skip terms where the exclusion population is less than Tx, as they are invalid.
      if (exclusion_population >= Tx) {
        # Step 9: Calculate the Hypergeometric Term
        # Compute the hypergeometric probability term for the given parameters.
        hypergeo_term <- choose_custom_log(exclusion_population, Tx) / choose_custom_log(N, Tx)

        # Step 10: Derivative of Hypergeometric Term
        # Compute the derivative of the hypergeometric term using the digamma function.
        d_hypergeo_Tx <- hypergeo_term * (
          safe_digamma(N - Tx - (k - (ell - i)) * barN + 1) - safe_digamma(N - Tx + 1)
        )

        # Step 11: Contribution from Inner Summation
        # Calculate the contribution of the current term in the inner summation.
        term <- (-1)^i * choose_custom_log(ell, i) * d_hypergeo_Tx
        inner_sum <- inner_sum + term  # Accumulate the result.
      } else {
        if (verbose) {
          # Step 12: Skip Invalid Terms
          # Print a message for skipped terms with invalid exclusion populations.
          cat("Skipping invalid term: ell =", ell, "i =", i,
              "| exclusion_population =", exclusion_population,
              "| Tx =", Tx, "\n")
        }
      }
    }

    # Step 13: Combine Inner Sum with Outer Sum
    # Combine the result of the inner summation with the outer summation.
    dq_Delta_Lambda_k <- dq_Delta_Lambda_k + choose_custom_log(k, ell) * inner_sum *
      (1 - delta)^ell * lambda^(k - ell)
  }

  # Step 14: Return Final Result
  # Return the calculated derivative value.
  return(dq_Delta_Lambda_k)
}


#' PMF for Imperfect Group-Level Sensitivity and Specificity
#'
#' Computes the probability mass function (PMF) for group-level sensitivity and specificity,
#' accounting for false negative rate (\eqn{\Delta}) and false positive rate (\eqn{\Lambda}).
#'
#' @param ty Integer. Number of identified contaminated groups (out of \eqn{b} total groups).
#' @param N Integer. Total population size.
#' @param barN Numeric. Average size of each group.
#' @param Tx Integer. Number of successes (positive cases) to be distributed in the population.
#' @param b Integer. Total number of groups.
#' @param delta Numeric. Group-level false negative rate (\eqn{\Delta}), the probability that a positive case is classified as negative (0 ≤ delta ≤ 1).
#' @param lambda Numeric. Group-level false positive rate (\eqn{\Lambda}), the probability that a negative case is classified as positive (0 ≤ lambda ≤ 1).
#' @param verbose Logical. If TRUE, prints conditions that are not fulfilled.
#' @return Numeric. The computed PMF value adjusted for imperfect group-level sensitivity and specificity.
#' @examples
#' # Example usage:
#' pmfHG.imperfect.group(ty = 3, N = 100, barN = 10, Tx = 20, b = 10, delta = 0.1, lambda = 0.05)
#' @export
pmfHG.imperfect.group <- function(ty, N, barN, Tx, b, delta, lambda, verbose = TRUE) {

  pmf <- 0 # % Initialize the probability to zero

  if (ty > b || ty > Tx && lambda==1) {  # ty must be <= min(b, Tx) when lambda ==1
    return(c(pmf))
  }

  # Binomial coefficient for selecting t_y groups from b total groups
  outer_binomial <- choose_custom_log(b, ty)

  # Initialize the summation for inclusion-exclusion
  sum_inclusion_exclusion <- 0

  # Loop over possible values of j (number of excluded groups from the b total groups)
  for (j in 0:ty) {
    # Calculate the population left after excluding certain groups
    exclusion_population <- N - barN * (b - j)  # Adjust population size based on the excluded groups

    # Check if the exclusion population is valid (should be >= Tx to proceed)
    if (exclusion_population >= Tx) {
      # Compute q(Delta, Lambda, k) which considers the sensitivity (Delta) and specificity (Lambda)
      q_Delta_Lambda <- qDLkGroup(b - j, N, Tx, barN, delta, lambda)

      # Calculate the term for inclusion-exclusion and adjust the sign based on j
      term <- (-1)^(ty - j) * choose_custom_log(ty, j) * q_Delta_Lambda

      # Add the term to the summation for inclusion-exclusion
      sum_inclusion_exclusion <- sum_inclusion_exclusion + term
    } else {

      if (verbose) {
      # Debug message for invalid terms where the exclusion population is not sufficient
      cat("Skipping invalid term: j =", j,
          "| exclusion_population =", exclusion_population,
          "| Tx =", Tx, "\n")
    }
  }
}

  # Multiply the result of the inclusion-exclusion summation by the outer binomial coefficient to get the final PMF
  pmf <- outer_binomial * sum_inclusion_exclusion

  # Apply a threshold for small probabilities (numerical stability)
  if (abs(pmf) < 1e-5) pmf <- 0

  # Return the final computed PMF
  return(pmf)
}




#' Derivative of PMF w.r.t Tx for Imperfect Group-Level Sensitivity and Specificity
#'
#' Computes derivative of the probability mass function (PMF) for group-level sensitivity and specificity,
#' w.r.t. Tx accounting for false negative rate (\eqn{\delta}) and false positive rate (\eqn{\lambda}).
#'
#' @param ty Integer. Number of identified contaminated groups (out of \eqn{b} total groups).
#' @param N Integer. Total population size.
#' @param barN Numeric. Average size of each group.
#' @param Tx Integer. Number of successes (positive cases) to be distributed in the population.
#' @param b Integer. Total number of groups.
#' @param delta Numeric. Group-level false negative rate (\eqn{\delta}), the probability that a positive case is classified as negative (0 ≤ delta ≤ 1).
#' @param lambda Numeric. Group-level false positive rate (\eqn{\lambda}), the probability that a negative case is classified as positive (0 ≤ lambda ≤ 1).
#' @param verbose Logical. If TRUE, prints conditions that are not fulfilled.
#' @return Numeric. The derivative of the computed PMF value w.r.t Tx adjusted for imperfect group-level sensitivity and specificity.
#' @examples
#' # Example usage:
#' pmfHG.imperfect.group(ty = 3, N = 100, barN = 10, Tx = 20, b = 10, delta = 0.1, lambda = 0.05)
#' @export
dpmfHG.imperfect.group <- function(ty, N, barN, Tx, b, delta, lambda, verbose = TRUE) {
  # Derivative of PMF with respect to Tx (following Equation 26)
  # This function calculates the derivative of the PMF for group-level sensitivity and specificity
  # based on the provided parameters. It handles the inclusion-exclusion principle and adjusts for
  # sensitivity (Delta) and specificity (Lambda) at the group level.

  # ty: Number of identified contaminated groups (out of total b groups)
  # N: Total population size (base size for hypergeometric distribution)
  # barN: Average size of each group
  # Tx: Number of successes (contaminated cases to be distributed)
  # b: Total number of groups
  # Delta: False negative rate (probability of a true positive being missed)
  # Lambda: False positive rate (probability of a true negative being falsely identified as positive)

  # Calculate the PMF first using the existing function
  pmf <- pmfHG.imperfect.group(ty, N, barN, Tx, b, delta, lambda,verbose)

  # Use a stricter threshold to avoid calculating derivatives for zero PMF
  if (pmf > 0) {
    # Proceed with derivative calculation only if PMF is non-zero

    # Binomial coefficient for selecting ty groups from the b total groups
    outer_binomial <- choose_custom_log(b, ty)

    # Initialize the summation for the inclusion-exclusion derivative
    sum_inclusion_exclusion_derivative <- 0

    # Loop over possible values of j (number of excluded groups)
    for (j in 0:ty) {
      # Compute the derivative of q(Delta, Lambda) with respect to Tx for the current value of j
      dq_Delta_Lambda_derivative <- dqDLkGroup(b - j, N, Tx, barN, delta, lambda,verbose)

      # Update the sum for the inclusion-exclusion principle, applying the appropriate sign
      sum_inclusion_exclusion_derivative <- sum_inclusion_exclusion_derivative +
        (-1)^(ty - j) * choose_custom_log(ty, j) * dq_Delta_Lambda_derivative
    }

    # Final derivative of the PMF: Multiply the outer binomial coefficient by the summation
    pmf_derivative <- outer_binomial * sum_inclusion_exclusion_derivative
  } else {
    if(verbose) {
    # Debug message: If PMF is effectively zero, skip derivative calculation
    cat("PMF is effectively zero for ty =", ty, ", Tx =", Tx,
        ". Returning derivative as 0.\n")
    }
    pmf_derivative <- 0  # If PMF is zero, the derivative is also zero

  }

  # Return the derivative of the PMF
  return(pmf_derivative)
}



#' Compute \eqn{q_{\delta \lambda k} } for Item-Level Sensitivity and Specificity
#'
#' Calculates the probability \eqn{q_{\delta \lambda k} } that none of the items
#' in the groups are identified as contaminated, considering item-level false
#' negative and false positive rates.
#'
#' @param k Integer. Number of groups considered in the summation.
#' @param N Integer. Total population size.
#' @param Tx Integer. Number of successes (positive cases) to be distributed.
#' @param barN Numeric. Average size of each group.
#' @param delta Numeric. False negative rate (item level).
#' @param lambda Numeric. False positive rate (item level).
#' @return Numeric. The computed value of \eqn{q_{\delta\lambda k} }.
#' @param verbose Logical. If TRUE, prints conditions that are not fulfilled.
#' @name qdlkItem
#' @rdname qdlkItem
#' @examples
#' qdlkqdlkItem(k = 3, N = 100, Tx = 5, barN = 10, delta = 0.1, lambda = 0.05)
#' @export
qdlkItem <- function(k, N, Tx, barN, delta, lambda, verbose=TRUE) {
  # Function to compute qδλk for item-level sensitivity and specificity (Equation 13)
  # qdlk: For Item case
  # qΔΛk : Probability that none are identified as contaminated
  # Function to calculate qΔΛk (Equation 10)
  # k: Number of groups considered in the summation.
  # N: Total population size.
  # Tx: Number of successes (positive cases) to be distributed.
  # barN: Average size of each group.
  # delta: False negative rate (item level).
  # lambda: False positive rate (item level).

  # Initialize q_delta_lambda_item_k to 0
  q_delta_lambda_item_k <- 0

  # Iterate over possible values of X_hat_k (number of contaminated items across k groups)
  # This sum accounts for all possible contaminations within k groups
  for (X_hat_k in 0:(k * barN)) {

    # Calculate exclusion population (those not in the groups)
    exclusion_population <- N - k * barN

    # Required number of items to be chosen from the exclusion population
    required_items <- Tx - X_hat_k

    # Check if the hypergeometric term is valid
    if (exclusion_population >= required_items && required_items >= 0) {

      # Hypergeometric term: choosing X_hat_k contaminated items from the k*barN possible items
      # and the remaining required items from the exclusion population
      hypergeo_term <- choose_custom_log(k * barN, X_hat_k) *
        (choose_custom_log(exclusion_population, required_items) / choose_custom_log(N, Tx))

      # Probability term: adjusting for false negative rate (delta) and false positive rate (lambda)
      probability_term <- (1 - delta)^X_hat_k * lambda^(k * barN - X_hat_k)

      # Add the product of the hypergeometric term and probability term to the sum
      q_delta_lambda_item_k <- q_delta_lambda_item_k + hypergeo_term * probability_term
    } else {
      if (verbose) {
        # Debug message for invalid terms
        cat("Skipping invalid term: X_hat_k =", X_hat_k,
            "| exclusion_population =", exclusion_population,
            "| required_items =", required_items,
            "| Tx =", Tx, "\n")
      }
    }
  }

  # Return the computed value of q_delta_lambda_item_k
  return(q_delta_lambda_item_k)
}

#' Compute Derivative of \eqn{q_{\delta\lambda k} } with respect to \eqn{T_x}
#'
#' Calculates the derivative of \eqn{q_{\delta\lambda k} } with respect to the
#' number of successes \eqn{T_x} for item-level sensitivity and specificity.
#'
#' @param k Integer. Number of groups considered in the summation.
#' @param N Integer. Total population size.
#' @param Tx Integer. Number of successes (positive cases) to be distributed.
#' @param barN Numeric. Average size of each group.
#' @param delta Numeric. False negative rate (item level).
#' @param lambda Numeric. False positive rate (item level).
#' @return Numeric. The derivative of \eqn{q_{\delta\lambda k} } with respect to \eqn{T_x}.
#' @name dqdlkItem
#' @rdname dqdlkItem
#' @examples
#' dqdlk(k = 3, N = 100, Tx = 5, barN = 10, delta = 0.1, lambda = 0.05)
#' @export
dqdlkItem <- function(k, N, Tx, barN, delta, lambda) {
  # Function to compute the derivative of q_delta_lambda_item_k
  # This is the derivative of the function q_delta_lambda_item_k (from Equation 27)
  # with respect to Tx for item-level sensitivity and specificity.
  # k: Number of groups considered in the summation.
  # N: Total population size.
  # Tx: Number of successes (positive cases) to be distributed.
  # barN: Average size of each group.
  # delta: False negative rate (item level).
  # lambda: False positive rate (item level).

  dq_delta_lambda_item_k <- 0  # Initialize derivative sum

  # Loop over possible contamination levels (X_hat_k), where X_hat_k is the number of contaminated items
  for (X_hat_k in 0:(k * barN)) {
    exclusion_population <- N - k * barN  # Exclusion population: those outside the k*barN contaminated groups
    required_items <- Tx - X_hat_k  # Number of required non-contaminated items to be selected from the exclusion population

    # Check if the hypergeometric term is valid (exclusion population must be sufficient)
    if (exclusion_population >= required_items && required_items >= 0) {

      # Hypergeometric term:
      # This calculates the probability of selecting X_hat_k contaminated items from the k*barN items
      # and the required number of non-contaminated items from the exclusion population
      # Based on PMF form
      # hypergeo_term <- choose_custom_log(k * barN, X_hat_k) * (choose_custom_log(exclusion_population, required_items) / choose_custom_log(N, Tx))
      # Based on alternative probability form
      hypergeo_term <- choose_custom_log(Tx, X_hat_k) * (choose_custom_log(N-Tx, k * barN - X_hat_k) / choose_custom_log(N, k * barN))

      # Derivative of the hypergeometric term with respect to Tx
      # The derivative is calculated using the digamma function, which is the derivative of the logarithm of the gamma function
      # This term accounts for how the probability changes as the number of contaminated items Tx changes
      # Old
     # d_hypergeo_Tx <- hypergeo_term * (
     #   safe_digamma(Tx - X_hat_k + 1) - safe_digamma(Tx + 1) +
     #     safe_digamma(N - Tx - k * barN + X_hat_k) - safe_digamma(N - Tx + 1)
     # )

     # New
     d_hypergeo_Tx <- hypergeo_term * (
       safe_digamma(Tx + 1) - safe_digamma(N - Tx + 1) + safe_digamma(N - Tx - k * barN + X_hat_k + 1)  - safe_digamma(Tx - X_hat_k + 1)
      )

      # Probability term:
      # This term accounts for the effect of false negatives (delta) and false positives (lambda) on the detection process
      probability_term <- (1 - delta)^X_hat_k * lambda^(k * barN - X_hat_k)

      # Combine the derivative of the hypergeometric term with the probability term
      dq_delta_lambda_item_k <- dq_delta_lambda_item_k + d_hypergeo_Tx * probability_term
    }
  }

  # Return the computed derivative
  return(dq_delta_lambda_item_k)
}

#' Probability Mass Function (PMF) for Imperfect Item-Level Sensitivity and Specificity
#'
#' Calculates the probability mass function (PMF) with item-level sensitivity and
#' specificity, based on the number of identified contaminated groups.
#'
#' @param ty Integer. Number of identified contaminated groups (out of total \eqn{b} groups).
#' @param N Integer. Total population size.
#' @param barN Numeric. Average size of each group.
#' @param Tx Integer. Number of successes (contaminated cases to be distributed).
#' @param b Integer. Total number of groups.
#' @param delta Numeric. False negative rate (item level probability of a true positive being missed).
#' @param lambda Numeric. False positive rate (item level probability of a true negative being falsely identified as positive).
#' @return Numeric. The computed probability mass function value.
#' @param verbose Logical. If TRUE, prints conditions that are not fulfilled.
#' @examples
#' pmfHG.imperfect.item(ty = 2, N = 100, barN = 10, Tx = 5, b = 5, delta = 0.1, lambda = 0.05)
#' @export
pmfHG.imperfect.item <- function(ty, N, barN, Tx, b, delta, lambda, verbose=TRUE) {
  # Function to calculate PMF (Probability Mass Function) with item-level sensitivity and specificity.
  # This corresponds to Equation 13 and computes the probability for a given number of contaminated items (ty).

  # ty: Number of identified contaminated groups (out of total b groups)
  # N: Total population size (base size for hypergeometric distribution)
  # barN: Average size of each group
  # Tx: Number of successes (contaminated cases to be distributed)
  # b: Total number of groups
  # delta: False negative rate (item level probability of a true positive being missed)
  # lambda: False positive rate (item level probability of a true negative being falsely identified as positive)

  pmf <- 0  # Initialize probability to zero

  if (ty > b || ty > Tx && lambda==1) {  # ty must be <= min(b, Tx)
    return(pmf)
  }

    # Binomial coefficient for selecting `ty` contaminated groups from `b` total groups
    outer_binomial <- choose_custom_log(b, ty)

    # Initialize the summation for inclusion-exclusion terms
    sum_inclusion_exclusion <- 0

    # Loop over possible values of `j` to compute inclusion-exclusion terms for each possible contamination level
    for (j in 0:ty) {
      # Adjust the population for excluded groups: those not included in the current set of contaminated groups
      exclusion_population <- N - barN * (b - j)

      # Check if the exclusion population is valid (i.e., there are enough items for the hypergeometric calculation)

      # Handle special cases
      if (exclusion_population < 0) {
        next  # Skip this iteration as the exclusion population is invalid
      } else if (exclusion_population == Tx && Tx == 0) {
        q_delta_lambda_item <- 1  # Set to 1 when both are zero
      } else if (exclusion_population >= Tx) {
        # Compute dq(Delta, Lambda, k) which considers the sensitivity (Delta) and specificity (Lambda)
        q_delta_lambda_item <- qdlkItem(b - j, N, Tx, barN, delta, lambda, verbose)
      } else {
        if (verbose) {
          # Debug message for invalid terms where the exclusion population is not sufficient
          cat("Skipping invalid term: j =", j,
              "| exclusion_population =", exclusion_population,
              "| Tx =", Tx, "\n")
        }
        next
      }

      # Compute the term for inclusion-exclusion
      # The term adjusts for how many items are included or excluded from the contamination set
      term <- (-1)^(ty - j) * choose_custom_log(ty, j) * q_delta_lambda_item

      # Add the term to the summation
      sum_inclusion_exclusion <- sum_inclusion_exclusion + term

      }


    # Multiply the binomial coefficient with the sum of inclusion-exclusion terms to compute the PMF
    pmf <- outer_binomial * sum_inclusion_exclusion

    # Apply a threshold for very small probabilities, setting them to zero to avoid numerical errors
    if (abs(pmf) < 1e-5) pmf <- 0



  # Return the computed PMF value
  return(pmf)
}




#' Derivative of PMF for Imperfect Item-Level Detection
#'
#' This function calculates the derivative of the probability mass function (PMF) for
#' imperfect item-level detection based on the provided parameters. It adjusts for
#' item-level sensitivity and specificity using the inclusion-exclusion principle.
#'
#' @param ty Integer. Number of identified contaminated groups (out of total `b` groups).
#' @param N Integer. Total population size (base size for hypergeometric distribution).
#' @param barN Numeric. Average size of each group.
#' @param Tx Integer. Number of successes (contaminated cases to be distributed).
#' @param b Integer. Total number of groups.
#' @param delta Numeric. False negative rate (item-level probability of a true positive being missed).
#' @param lambda Numeric. False positive rate (item-level probability of a true negative being falsely identified as positive).
#' @return Numeric. The derivative of the PMF with respect to `Tx`.
#' @examples
#' dpmfHG.imperfect.item(ty = 3, N = 100, barN = 10, Tx = 5, b = 10, delta = 0.1, lambda = 0.05)
#' @export
dpmfHG.imperfect.item <- function(ty, N, barN, Tx, b, delta, lambda, verbose=TRUE) {

  pmf <- pmfHG.imperfect.item(ty, N, barN, Tx, b, delta, lambda, verbose )

  if (ty > b || (ty > Tx && lambda == 1) || (pmf == 0 && lambda == 1)) {
    return(0)
  }

  #  if (pmf > 0) {
  outer_binomial <- choose_custom_log(b, ty)

  sum_inclusion_exclusion_derivative <- 0

  for (j in 0:ty) {
    exclusion_population <- N - barN * (b - j)

    # Handle special cases
    if (exclusion_population < 0) {
      next  # Skip this iteration as the exclusion population is invalid
    } else if (exclusion_population == Tx && Tx == 0) {
      dq_delta_lambda_item <- 1  # Set to 1 when both are zero
    } else if (exclusion_population >= Tx) {
      # Compute dq(Delta, Lambda, k) which considers the sensitivity (Delta) and specificity (Lambda)
      dq_delta_lambda_item <- dqdlkItem(b - j, N, Tx, barN, delta, lambda)
    } else {
      if (verbose) {
        # Debug message for invalid terms where the exclusion population is not sufficient
        cat("Skipping invalid term: j =", j,
            "| exclusion_population =", exclusion_population,
            "| Tx =", Tx, "\n")
      }
      next
    }



#    if (exclusion_population >= Tx) {
#      dq_delta_lambda_item <- dqdlkItem(b - j, N, Tx, barN, delta, lambda)
      term_derivative <- (-1)^(ty - j) * choose_custom_log(ty, j) * dq_delta_lambda_item
      sum_inclusion_exclusion_derivative <- sum_inclusion_exclusion_derivative + term_derivative
 #   }
  }

  pmf_derivative <- outer_binomial * sum_inclusion_exclusion_derivative
  #  } else {
  #    pmf_derivative <- 0


  return(pmf_derivative)
}

#' Fisher Information for Imperfect Detection
#'
#' This function calculates the Fisher Information for imperfect detection using different
#' methods (analytic, numerical, or PMF-based approximation) and detection types (group or item level).
#'
#' @param N Integer. Total population size (base size for hypergeometric distribution).
#' @param barN Numeric. Average size of each group.
#' @param Tx Integer. Number of successes (contaminated cases to be distributed).
#' @param b Integer. Total number of groups.
#' @param delta Numeric. False negative rate.
#' @param lambda Numeric. False positive rate.
#' @param method Character. Method for computing Fisher Information. One of `"AD"` (analytic derivative), `"ND"` (numerical derivative), or `"PMF"` (PMF-based approximation).
#' @param type Character. Detection type. One of `"group"` or `"item"`.
#' @param verbose Logical. If TRUE, prints conditions that are not fulfilled.
#' @return Numeric. The calculated Fisher Information.
#' @examples
#' FIpmfHG.Tx.imperfect(N = 100, barN = 10, Tx = 5, b = 10, delta = 0.1, lambda = 0.05, method = "AD", type = "item")
#' @export
FIpmfHG.Tx.imperfect <- function(N, barN, Tx, b, delta, lambda, method = c("AD", "ND", "PMF"), type = c("group", "item"),verbose=TRUE) {
  FI_Tx <- 0
  method <- match.arg(method)
  type <- match.arg(type)

  for (ty in 0:b) {
    if (ty > b || ty > Tx && lambda == 1) next

    P_ty <- if (type == "group") {
      pmfHG.imperfect.group(ty, N, barN, Tx, b, delta, lambda,verbose)
    } else {
      pmfHG.imperfect.item(ty, N, barN, Tx, b, delta, lambda,verbose)
    }

    if (P_ty == 0) next

    if (method == "AD") {
      dP_Tx <- if (type == "group") {
        dpmfHG.imperfect.group(ty, N, barN, Tx, b, delta, lambda,verbose)
      } else {
        dpmfHG.imperfect.item(ty, N, barN, Tx, b, delta, lambda, verbose)
      }

      FI_Tx <- FI_Tx + (dP_Tx^2 / P_ty)
    }

    if (method == "PMF") {
      P_ty_Tx_plus1 <- if (type == "group") {
        pmfHG.imperfect.group(ty, N, barN, Tx + 1, b, delta, lambda, verbose)
      } else {
        pmfHG.imperfect.item(ty, N, barN, Tx + 1, b, delta, lambda, verbose)
      }

      if (P_ty_Tx_plus1 > 0 && P_ty > 0) {
        FI_Tx <- FI_Tx + 4 * (sqrt(P_ty_Tx_plus1) - sqrt(P_ty))^2
      }
    }

    if (method == "ND") {
      P_plus <- if (type == "group") {
        pmfHG.imperfect.group(ty, N, barN, Tx + 1, b, delta, lambda, verbose)
      } else {
        pmfHG.imperfect.item(ty, N, barN, Tx + 1, b, delta, lambda, verbose)
      }

      if (P_plus > 0 && P_ty > 0) {
        FI_Tx <- FI_Tx + ((P_plus - P_ty)^2 / P_ty)
      }
    }
  }

  return(FI_Tx)
}

#' Fisher Information for Imperfect Detection
#'
#' This function calculates the Fisher Information for imperfect detection using different
#' methods (analytic, numerical, or PMF-based approximation) and detection types (group or item level).
#' This function also provides probability, derivative of ty w.r.t Tx and finally Fisher information for a given Tx
#' @param N Integer. Total population size (base size for hypergeometric distribution).
#' @param barN Numeric. Average size of each group.
#' @param Tx Integer. Number of successes (contaminated cases to be distributed).
#' @param b Integer. Total number of groups.
#' @param delta Numeric. False negative rate.
#' @param lambda Numeric. False positive rate.
#' @param method Character. Method for computing Fisher Information. One of `"AD"` (analytic derivative), `"ND"` (numerical derivative), or `"PMF"` (PMF-based approximation).
#' @param type Character. Detection type. One of `"group"` or `"item"`.
#' @param verbose Logical. If TRUE, prints conditions that are not fulfilled.
#' @return Numeric. The calculated Fisher Information.
#' @examples
#' FIpmfHG.Tx.imperfect(N = 100, barN = 10, Tx = 5, b = 10, delta = 0.1, lambda = 0.05, method = "AD", type = "item")
#' @export
FIpmfHG.Tx.imperfect.detail <- function(N, barN, Tx, b, delta, lambda, method = c("AD", "ND", "PMF"), type = c("group", "item"), verbose = TRUE) {
  FI_Tx <- 0
  method <- match.arg(method)
  type <- match.arg(type)
  P_ty_vector <- numeric(b + 1)  # Initialize vector with zeros
  dP_ty_vector <- NULL

  for (ty in 0:b) {
    if (ty > b || (ty > Tx && lambda == 1)) next  # Skip iterations where ty > b or ty > Tx when lambda == 1

    P_ty <- if (ty > Tx && lambda == 1) {  # Ensure P_ty = 0 when ty > Tx
      0
    } else if (type == "group") {
      pmfHG.imperfect.group(ty, N, barN, Tx, b, delta, lambda, verbose)
    } else {
      pmfHG.imperfect.item(ty, N, barN, Tx, b, delta, lambda, verbose)
    }

    P_ty_vector[ty + 1] <- P_ty  # Store probability

    if (P_ty == 0) next  # Skip further calculations if P_ty is 0

    if (method == "AD") {
      dP_Tx <- if (type == "group") {
        dpmfHG.imperfect.group(ty, N, barN, Tx, b, delta, lambda, verbose)
      } else {
        dpmfHG.imperfect.item(ty, N, barN, Tx, b, delta, lambda, verbose)
      }

      dP_ty_vector[ty + 1] <- dP_Tx
      FI_Tx <- FI_Tx + (dP_Tx^2 / P_ty)
    }

    if (method == "PMF") {
      P_ty_Tx_plus1 <- if (type == "group") {
        pmfHG.imperfect.group(ty, N, barN, Tx + 1, b, delta, lambda, verbose)
      } else {
        pmfHG.imperfect.item(ty, N, barN, Tx + 1, b, delta, lambda, verbose)
      }

      if (P_ty_Tx_plus1 > 0 && P_ty > 0) {
        FI_Tx <- FI_Tx + 4 * (sqrt(P_ty_Tx_plus1) - sqrt(P_ty))^2
      }
    }

    if (method == "ND") {
      P_plus <- if (type == "group") {
        pmfHG.imperfect.group(ty, N, barN, Tx + 1, b, delta, lambda, verbose)
      } else {
        pmfHG.imperfect.item(ty, N, barN, Tx + 1, b, delta, lambda, verbose)
      }

      if (P_plus > 0 && P_ty > 0) {
        FI_Tx <- FI_Tx + ((P_plus - P_ty)^2 / P_ty)
      }
    }
  }

  # Create the data frame ensuring p_ty = 0 for ty > Tx
  ty.data <- data.frame(ty = 0:b, p_ty = P_ty_vector)

  if (method == "AD") {
    ty.data$dp_ty <- c(dP_ty_vector, rep(NA, length(ty.data$ty) - length(dP_ty_vector)))
  }

  list(FI_Tx = FI_Tx, ty.data = ty.data)
}


#' Calculate Expectation and Variance with Group-Level Sensitivity and Specificity
#'
#' Computes the expectation and variance of positive test groups using
#' hypergeometric probabilities at the group level.
#'
#' @param N Total number of items.
#' @param Tx Number of contaminated items.
#' @param barN Number of items per group.
#' @param b Number of groups.
#' @param delta Sensitivity of the test.
#' @param lambda Specificity of the test.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{V}{Probability that a single group tests positive.}
#'   \item{W}{Covariance between two groups.}
#'   \item{expectation}{Expected number of positive groups.}
#'   \item{variance}{Variance of the number of positive groups.}
#' }
#'
#' @export
ExpVarHG.imperfect.group <- function(N, Tx, barN, b, delta, lambda) {
# Step 1: Calculate hypergeometric probabilities ----------------#
# Probability that a single group does not contain any contaminated items
P_single_group <- choose_custom_log(N - barN, Tx) / choose_custom_log(N, Tx)

# Probability that two groups do not contain any contaminated items
P_two_groups <- choose_custom_log(N - 2 * barN, Tx) / choose_custom_log(N, Tx)

# Step 2: Calculate V (probability a single group tests positive) ----------------#
# Combines true positive detection and false positive errors
V <- delta * (1 - P_single_group) + (1 - lambda) * P_single_group

# Step 3: Calculate W (covariance between two groups) ----------------#
# Dependency between two groups being positive
W <- (P_two_groups - P_single_group^2) * (1 - lambda - delta)^2

# Step 4: Calculate expectation ----------------#
expectation <- b * V  # Total number of positive groups

# Step 5: Calculate variance ----------------#
variance <- b * V * (1 - V) * (1 + (b - 1) * W / (V * (1 - V)))

# Return calculated values
return(list(V = V, W = W, expectation = expectation, variance = variance))
}

#' Calculate Expectation and Variance with Item-Level Sensitivity and Specificity
#'
#' Computes the expectation and variance of positive test groups using
#' hypergeometric probabilities at the item level.
#'
#' @param N Total number of items.
#' @param Tx Number of contaminated items.
#' @param barN Number of items per group.
#' @param b Number of groups.
#' @param delta Sensitivity of the test.
#' @param lambda Specificity of the test.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{V}{Probability that a single group tests positive.}
#'   \item{W}{Covariance between two groups.}
#'   \item{expectation}{Expected number of positive groups.}
#'   \item{variance}{Variance of the number of positive groups.}
#' }
#'
#' @export
ExpVarHG.imperfect.item <- function(N, Tx, barN, b, delta, lambda) {
  # Helper function to compute Phi_Xk(s) ----------------#
  compute_phi <- function(k, N, Tx, barN, s) {
    phi <- 0
    for (X_k in 0:(k * barN)) {
      # Hypergeometric term: exact contaminated item count
      hypergeo_term <- choose_custom_log(k * barN, X_k) *
        (choose_custom_log(N - k * barN, Tx - X_k) / choose_custom_log(N, Tx))

      # Weighted by sensitivity and specificity
      phi <- phi + hypergeo_term * s^X_k
    }
    return(phi)
  }

  # Step 2: Compute Phi_X1 and Phi_X2 ----------------#
  phi_X1 <- compute_phi(1, N, Tx, barN, (1 - delta) / lambda)
  phi_X2 <- compute_phi(2, N, Tx, barN, (1 - delta) / lambda)

  # Step 3: Calculate V ----------------#
  V <- 1 - lambda^barN * phi_X1

  # Step 4: Calculate W ----------------#
  W <- lambda^(2 * barN) * phi_X2 - (lambda^barN * phi_X1)^2

  # Step 5: Calculate expectation ----------------#
  expectation <- b * V

  # Step 6: Calculate variance ----------------#
  variance <- b * V * (1 - V) * (1 + (b - 1) * W / (V * (1 - V)))

  # Return calculated values
  return(list(V = V, W = W, expectation = expectation, variance = variance))
}

#' Unified Expectation and Variance Calculation for Imperfect Tests
#'
#' Computes the expectation and variance of positive test groups based on
#' either group-level or item-level sensitivity and specificity.
#'
#' @param N Total number of items.
#' @param Tx Number of contaminated items.
#' @param barN Number of items per group.
#' @param b Number of groups.
#' @param delta Sensitivity of the test.
#' @param lambda Specificity of the test.
#' @param type Character string specifying whether to use "group" or "item" calculations.
#' Default is \code{"group"}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{V}{Probability that a single group tests positive.}
#'   \item{W}{Covariance between two groups.}
#'   \item{expectation}{Expected number of positive groups.}
#'   \item{variance}{Variance of the number of positive groups.}
#' }
#'
#' @export
ExpVarHG.imperfect <- function(N, Tx, barN, b, delta, lambda, type = c("group", "item")) {
  # Validate input type
  type <- match.arg(type)

  if (type == "group") {
    # Group-level calculation
    ExpVarHG <- ExpVarHG.imperfect.group(N, Tx, barN, b, delta, lambda)
  } else {
    # Item-level calculation
    ExpVarHG <- ExpVarHG.imperfect.item(N, Tx, barN, b, delta, lambda)
  }

  # Return calculated results
  return(list(
    V = ExpVarHG$V,
    W = ExpVarHG$W,
    expectation = ExpVarHG$expectation,
    variance = ExpVarHG$variance
  ))
}
