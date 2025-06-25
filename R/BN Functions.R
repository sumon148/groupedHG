#' Probability Mass Function for Grouped Binomial Distribution under Perfect Test
#'
#' @description
#' Computes the probability mass function (PMF) for a grouped Binomial distribution
#' under the assumption of a perfect test (no false positives or false negatives).
#'
#' @param ty Integer. Number of groups where the target outcome is observed.
#' @param Tx Integer. Number of target individuals in the population.
#' @param b Integer. Total number of groups.
#' @param barN Integer. Average number of individuals per group.
#' @param N Integer. Total population size.
#'
#' @return Numeric. The probability mass function (PMF) value for the given parameters.
#'
#' @details
#' The function calculates the probability of observing `ty` positive groups in `b` total groups,
#' using binomial probabilities for detecting target individuals in a group.
#'
#' @examples
#' pmf_bn_perfect(ty = 3, Tx = 10, b = 5, barN = 2, N = 20)
#'
#' @export
pmf_bn_perfect <- function(ty, Tx, b, barN, N) {
  # Function to calculate the PMF for a grouped Binomial distribution under Perfect Test (Equation 20)
  # Parameters:
  # ty: Number of groups where the target outcome is observed
  # N: Total population size
  # barN: Average number of individuals per group
  # Tx: Number of target individuals in the population
  # b: Total number of groups

  if (ty > b || ty < 0 || Tx < 0 || N < barN) stop("Invalid input parameters.")

  # Probabilities of detection and non-detection
  prob_non_detect <- choose_custom(N-barN,Tx) / choose_custom(N,Tx)
  prob_detect <- 1 - prob_non_detect

  # Binomial coefficient and product terms
  outer_binomial <- choose(b, ty)
  pmf <- outer_binomial * prob_detect^ty * prob_non_detect^(b - ty)

  # Apply threshold for small probabilities
  if (abs(pmf) < 1e-12) pmf <- 0

  return(pmf)
}


#' Probability Mass Function for Grouped Binomial Distribution under group-level imperfect Test
#'
#' @description
#' Computes the probability mass function (PMF) for a grouped Binomial distribution
#' under the assumption of group-level imperfect test (false positives or false negatives).
#' @param ty Numeric. The number of groups where the target outcome is observed.
#' @param Tx Numeric. The number of target individuals in the population.
#' @param b Numeric. The total number of groups.
#' @param barN Numeric. The average number of individuals per group.
#' @param N Numeric. The total population size.
#' @param delta Numeric. The group-level false negative rate, which is the probability that a positive case is classified as negative.
#' @param lambda Numeric. The group-level false positive rate, which is the probability that a negative case is classified as positive.
#' @return Numeric. The probability mass function (PMF) for the given group-level parameters.
#' @examples
#' pmf_bn_group_imperfect(ty = 5, Tx = 20, b = 10, barN = 50, N = 1000, delta = 0.1, lambda = 0.05)
#' @export
pmf_bn_group_imperfect <- function(ty, Tx, b, barN, N, delta, lambda) {
  # Function to calculate group-level PMF (Equation 21)

  # Parameters:
  # ty: Number of groups where the target outcome is observed
  # N: Total population size
  # barN: Average number of individuals per group
  # Tx: Number of target individuals in the population
  # b: Total number of groups
  # Delta: Group level False negative rate (probability that a positive case is classified as negative)
  # Lambda: Group level False positive rate (probability that a negative case is classified as positive)


  pmf <- 0 # % Initialize the probability to zero

  if (Tx > N) return(pmf)



  # Probabilities of detection and non-detection
  prob_non_detect <- lambda * choose_custom(N-barN,Tx) / choose_custom(N,Tx) + (1 - delta) * (1 - choose_custom(N-barN,Tx) / choose_custom(N,Tx))
  prob_detect <- 1 - prob_non_detect

  # Binomial coefficient and product terms
  outer_binomial <- choose(b, ty)
  pmf <- outer_binomial * prob_detect^ty * prob_non_detect^(b - ty)

  # Apply threshold for small probabilities
  if (abs(pmf) < 1e-12) pmf <- 0


  # Check for invalid values
  if (is.na(pmf) || is.nan(pmf) || is.infinite(pmf)) {
    return(0)  # Prevent NaN issues
  }

  return(pmf)


}

#' Probability Mass Function for Grouped Binomial Distribution under item-level imperfect Test
#'
#' @description
#' Computes the probability mass function (PMF) for a grouped Binomial distribution
#' under the assumption of item-level imperfect test (false positives or false negatives).
#' @param ty Numeric. The number of groups where the target outcome is observed.
#' @param Tx Numeric. The number of target individuals in the population.
#' @param b Numeric. The total number of groups.
#' @param barN Numeric. The average number of individuals per group.
#' @param N Numeric. The total population size.
#' @param delta Numeric. The group-level false negative rate, which is the probability that a positive case is classified as negative.
#' @param lambda Numeric. The group-level false positive rate, which is the probability that a negative case is classified as positive.
#' @return Numeric. The probability mass function (PMF) for the given group-level parameters.
#' @examples
#' pmf_bn_item_imperfect(ty = 5, Tx = 20, b = 10, barN = 50, N = 1000, delta = 0.1, lambda = 0.05)
#' @export
pmf_bn_item_imperfect <- function(ty, Tx, b, barN, N, delta, lambda) {
  # Function to calculate item-level PMF using Equation (22)
  # Parameters:
  # ty: Number of groups where the target outcome is observed
  # N: Total population size
  # barN: Average number of individuals per group
  # Tx: Number of target individuals in the population
  # b: Total number of groups
  # Delta: Item level False negative rate (probability that a positive case is classified as negative)
  # Lambda: Item level False positive rate (probability that a negative case is classified as positive)

  # Compute detection and non-detection probabilities
  prob_non_detect <- ( lambda*(1-Tx / N) + (1-delta) * (Tx / N) )^barN
  prob_detect <- 1 - prob_non_detect

  # Binomial coefficient and product of probabilities
  outer_binomial <- choose(b, ty)
  pmf <- outer_binomial * prob_detect^ty * prob_non_detect^(b - ty)

  # Apply threshold for small probabilities
  if (abs(pmf) < 1e-12) pmf <- 0

  return(pmf)
}


#' Fisher Information for Grouped Binomial (Perfect Case)
#'
#' This function calculates the Fisher information for a grouped binomial distribution in the perfect case
#' (when detection probability is assumed to be perfect), based on the number of target individuals in the population.
#' The formula for Fisher Information is derived for this scenario, and the function returns the corresponding value.
#'
#' @param Tx Numeric. The number of target individuals in the population.
#' @param b Numeric. The total number of groups.
#' @param barN Numeric. The average number of individuals per group.
#' @param N Numeric. The total population size.
#'
#' @return Numeric. The Fisher Information for the grouped binomial distribution in the perfect case.
#' If the parameters are invalid or special conditions occur (such as Tx = 0 or numerical instability),
#' the function returns 0 and prints a message indicating the issue.
#'
#' @details
#' The Fisher Information is calculated using a formula based on the parameters of the grouped binomial model.
#' The function handles special cases, such as when the number of target individuals (\code{Tx}) is 0, by setting the Fisher information to 0.
#' It also checks for numerical instability by verifying if the denominator is close to zero and, in such cases, sets the Fisher Information to 0.
#'
#' @examples
#' # Example calculation for valid parameters
#' FIpmfBN.Tx.perfect(Tx = 10, b = 5, barN = 50, N = 1000)
#'
#' # Example when Tx = 0
#' FIpmfBN.Tx.perfect(Tx = 0, b = 5, barN = 50, N = 1000)
#'
#' @export
FIpmfBN.Tx.perfect <- function(Tx, b, barN, N) {
  # Fisher Information: Grouped Binomial (Perfect Case)

  # Parameters:
  # Tx: Number of target individuals in the population
  # b: Total number of groups
  # barN: Average number of individuals per group
  # N: Total population size

  # Input validation
  if (Tx < 0 || Tx > N || b <= 0 || barN <= 0 || N <= 0) {
    stop("Invalid input parameters. Ensure Tx is between 0 and N, and b, barN, N are positive.")
  }

  if (Tx == 0) {
    # Handle special case where Tx = 0
    FI <- 0  # Fisher information is zero as no variability exists
    cat("Tx = 0: Fisher information set to 0 as no variability exists.\n")
  } else {
    # General case
    FI.num <- b * barN^2 * (1 - Tx / N)^(barN - 2)
    FI.denom <- (1 - (1 - Tx / N)^barN) * N^2

    # Check for numerical stability
    if (abs(FI.denom) < 1e-12) {
      FI <- 0  # Avoid division by near-zero values
      cat("Numerical instability detected: Fisher information set to 0.\n")
    } else {
      FI <- FI.num / FI.denom
    }
  }

  return(FI)
}


#' Fisher Information for Grouped Binomial (Group-level Imperfect Case)
#'
#' This function calculates the Fisher information for a grouped binomial distribution in the imperfect case,
#' where both false positive and false negative rates are taken into account at the group level.
#' @param Tx Numeric. The number of target individuals in the population.
#' @param b Numeric. The total number of groups.
#' @param barN Numeric. The average number of individuals per group.
#' @param N Numeric. The total population size.
#' @param delta Numeric. The group-level false negative rate, which is the probability that a positive case is classified as negative.
#' @param lambda Numeric. The group-level false positive rate, which is the probability that a negative case is classified as positive.
#'
#' @return Numeric. The Fisher information for the grouped binomial distribution in the imperfect case. If the parameters lead to special conditions such as \( Tx = 0 \) or numerical instability, the function returns 0 and prints a corresponding message.
#'
#' @details
#' The Fisher information in the imperfect case incorporates both false positive and false negative rates. It is computed using the given parameters and is derived from the formula provided in Equation 29.
#' Special handling is included for the case where \( Tx = 0 \) or when numerical instability occurs due to small denominator values. In these cases, the Fisher information is set to 0 and a message is printed.
#'
#' @examples
#' # Example calculation for valid parameters
#' FIpmfBN.Tx.imperfect.group(Tx = 10, b = 5, barN = 50, N = 1000, delta = 0.1, lambda = 0.05)
#'
#' # Example when Tx = 0
#' FIpmfBN.Tx.imperfect.group(Tx = 0, b = 5, barN = 50, N = 1000, delta = 1, lambda = 1)
#'
#' @export
FIpmfBN.Tx.imperfect.group <- function(Tx, b, barN, N, delta, lambda) {
  # Fisher Information: Grouped binomial: Imperfect Case (29)

  # Parameters:
  # N: Total population size
  # barN: Average number of individuals per group
  # Tx: Number of target individuals in the population
  # b: Total number of groups
  # delta: Group level False negative rate (probability that a positive case is classified as negative)
  # lambda: Group level False positive rate (probability that a negative case is classified as positive)

  if (Tx == 0 || Tx == 100 & delta == 1 & lambda == 1) {
    # Handle Tx = 0 explicitly
    FI <- 0  # Fisher information is zero if Tx = 0 since no variability exists
    cat("Tx = 0: Fisher information set to 0 as no variability exists.\n")
  } else {
    # General case
    Common.part <-  lambda * (1 - Tx / N)^{barN} + (1 - delta) * (1 - (1 - Tx / N)^barN)
    FI.denom <- (1 - Common.part) * N^2
    FI.num <- b * barN^2 * Common.part^{-1} * ((1 - delta - lambda) * (1 - Tx / N)^{barN - 1})^2

    # Check for numerical stability in denominator
    if (abs(FI.denom) < 1e-12) {
      FI <- 0  # Avoid division by near-zero values
      cat("Numerical instability detected: Fisher information set to 0.\n")
    } else {
      FI <- FI.num / FI.denom
    }
  }

  return(FI)
}


#' Fisher Information for Grouped Binomial (Imperfect Case, Item Level)
#'
#' This function calculates the Fisher information for a grouped binomial distribution in the imperfect case
#' at the item level, where both false positive and false negative rates are considered. The formula is derived
#' based on Equation 30 from the referenced document.
#'
#' @param Tx Numeric. The number of target individuals in the population.
#' @param b Numeric. The total number of groups.
#' @param barN Numeric. The average number of individuals per group.
#' @param N Numeric. The total population size.
#' @param delta Numeric. The item-level false negative rate, which is the probability that a positive case is classified as negative.
#' @param lambda Numeric. The item-level false positive rate, which is the probability that a negative case is classified as positive.
#'
#' @return Numeric. The Fisher information for the grouped binomial distribution in the imperfect case at the item level.
#' If the parameters lead to special conditions such as \( Tx = 0 \) or numerical instability, the function returns 0 and prints a corresponding message.
#'
#' @details
#' The Fisher information in the imperfect case at the item level incorporates both false positive and false negative rates.
#' It is computed using the given parameters and is derived from the formula provided in Equation 30.
#' Special handling is included for the case where \( Tx = 0 \) or when numerical instability occurs due to small denominator values.
#' In these cases, the Fisher information is set to 0 and a message is printed.
#'
#' @examples
#' # Example calculation for valid parameters
#' FIpmfBN.Tx.imperfect.item(Tx = 10, b = 5, barN = 50, N = 1000, delta = 0.1, lambda = 0.05)
#'
#' # Example when Tx = 0
#' FIpmfBN.Tx.imperfect.item(Tx = 0, b = 5, barN = 50, N = 1000, delta = 1, lambda = 1)
#'
#' @export
FIpmfBN.Tx.imperfect.item <- function(Tx, b, barN, N,delta,lambda) {
  # Fisher Information: Grouped binomial: Imperfect Case (30)

  # Parameters:
  # N: Total population size
  # barN: Average number of individuals per group
  # Tx: Number of target individuals in the population
  # b: Total number of groups
  # delta: Item level False negative rate (probability that a positive case is classified as negative)
  # lambda: Item level False positive rate (probability that a negative case is classified as positive)


  if (Tx == 0 & delta==1 & lambda==1) {
    # Handle Tx = 0 explicitly
    FI <- 0  # Fisher information is zero if Tx = 0 since no variability exists
    cat("Tx = 0: Fisher information set to 0 as no variability exists.\n")
  } else {
    # General case
    Common.part <-  lambda*(1 - Tx/N) + (1 - delta) *   Tx / N
    FI.denom <- (1 - Common.part^{barN}) * N^2
    FI.num <- b * barN^2 * Common.part^{barN-2} * (1 - delta - lambda)^{2}

    # Check for numerical stability in denominator
    if (abs(FI.denom) < 1e-12) {
      FI <- 0  # Avoid division by near-zero values
      cat("Numerical instability detected: Fisher information set to 0.\n")
    } else {
      FI <- FI.num / FI.denom
    }
  }

  return(FI)
}


#' Fisher Information for Grouped Binomial (Group or Item Level Imperfect Case)
#'
#' This function calculates the Fisher information for a grouped binomial distribution in the imperfect case.
#' The Fisher information is calculated based on different methods (PMF-based approximation, numerical derivative, or analytical derivative)
#' and for either group-level or item-level models. The function computes the Fisher information over a range of possible values
#' for the number of target outcomes observed (ty) and returns the total Fisher information for the given parameters.
#'
#' @param N Numeric. The total population size.
#' @param barN Numeric. The average number of individuals per group.
#' @param Tx Numeric. The number of target individuals in the population.
#' @param b Numeric. The total number of groups.
#' @param delta Numeric. The group or item-level false negative rate (probability that a positive case is classified as negative).
#' @param lambda Numeric. The group or item-level false positive rate (probability that a negative case is classified as positive).
#' @param method Character. Method for computing Fisher Information. One of `"AD"` (analytic derivative), `"ND"` (numerical derivative), `"PMF-SM"` (PMF-based approximation using Sanchez-Moreno et al (2009) paper where ty=0 to (b-1)), or `"PMF-HI"` (where ty=0 to b used as in Shemyakin (2023)).
#' @param type Character. The type of model for which the Fisher information is calculated. One of "group" (group-level model) or "item" (item-level model). Default is "group".
#'
#' @return Numeric. The Fisher information for the grouped binomial distribution in the imperfect case.
#' If special conditions occur, such as division by zero, the function skips those terms and continues the calculation.
#'
#' @details
#' The Fisher information is computed based on three methods:
#' \itemize{
#'   \item "PMF": PMF-based approximation method, where the Fisher information is calculated using probabilities from the PMF function.
#'   \item "ND": Numerical derivative method, where the Fisher information is computed using finite differences of the PMF function.
#'   \item "AD": Analytical derivative method, which uses previously defined Fisher information formulas for group and item-level models.
#' }
#' The Fisher information is accumulated over all possible values of the number of target outcomes observed (\code{ty}),
#' from 0 to \code{b}, where \code{b} is the total number of groups. The appropriate PMF function is used based on the \code{type} parameter
#' (group-level or item-level).
#'
#' Special cases are handled where \code{ty} exceeds \code{Tx} and when \code{Tx} equals 0 with perfect detection rates (i.e., both \code{delta} and \code{lambda} equal 1).
#' The function also skips terms where the PMF is 0 to avoid division by zero errors.
#'
#' @examples
#' # Example using PMF method for group-level model
#' info_bn_tx_imperfect(N = 1000, barN = 50, Tx = 10, b = 5, delta = 0.1, lambda = 0.05, method = "PMF", type = "group")
#'
#' # Example using numerical derivative method for item-level model
#' info_bn_tx_imperfect(N = 1000, barN = 50, Tx = 10, b = 5, delta = 0.1, lambda = 0.05, method = "ND", type = "item")
#'
#' @export
info_bn_tx_imperfect <- function(N, barN, Tx, b, delta, lambda,
                                 method = c("AD", "ND", "PMF-SM", "PMF-HI"),
                                 type = c("group", "item")) {

  # Ensure valid method and type
  method <- match.arg(method)
  type <- match.arg(type)

  # Initialize Fisher Information
  FI_Tx <- 0

  if(method == "ND" | method == "PMF-SM" | method == "PMF-HI") {
    for (ty in 0:b) {
      # Skip invalid cases where ty > Tx and perfect detection
      if (ty > Tx && delta == 1 && lambda == 1) {
        next
      }

      # PMF for ty
      if (type == "group") {
        P_ty <- pmf_bn_group_imperfect(ty, Tx, b, barN, N, delta, lambda)
      } else {
        P_ty <- pmf_bn_item_imperfect(ty, Tx, b, barN, N, delta, lambda)
      }

      # Skip terms with P(t_y | T_x) = 0 to avoid division by zero
      if (P_ty == 0) {
        next
      }

      # Method: PMF-based approximation (PMF) - Using ty=0:b
      if (method == "PMF-HI") {
        # PMF values for Tx and Tx + 1
        if (type == "group") {
          P_ty_Tx_plus1 <- pmf_bn_group_imperfect(ty, Tx+1, b, barN, N, delta, lambda)
          P_ty_Tx <- pmf_bn_group_imperfect(ty, Tx, b, barN, N, delta, lambda)
        } else {
          P_ty_Tx_plus1 <- pmf_bn_item_imperfect(ty, Tx+1, b, barN, N, delta, lambda)
          P_ty_Tx <- pmf_bn_item_imperfect(ty, Tx, b, barN, N, delta, lambda)
        }

        if (P_ty_Tx_plus1 > 0 && P_ty_Tx > 0) {  # Avoid division by zero
          FI_Tx <- FI_Tx + 4 * (sqrt(P_ty_Tx_plus1) - sqrt(P_ty_Tx))^2
        }
      }

      # Method: PMF-based approximation (PMF) - Using ty=0:(b-1)
      if (method == "PMF-SM") {
        # PMF values for Tx and Tx + 1
        if (type == "group") {
          P_ty_Tx_plus1 <- pmf_bn_group_imperfect(ty, Tx+1, b, barN, N, delta, lambda)
          P_ty_Tx <- pmf_bn_group_imperfect(ty, Tx, b, barN, N, delta, lambda)
        } else {
          P_ty_Tx_plus1 <- pmf_bn_item_imperfect(ty, Tx+1, b, barN, N, delta, lambda)
          P_ty_Tx <- pmf_bn_item_imperfect(ty, Tx, b, barN, N, delta, lambda)
        }

        if (P_ty_Tx_plus1 > 0 && P_ty_Tx > 0 && ty < b) {  # Avoid division by zero
          FI_Tx <- FI_Tx + 4 * (sqrt(P_ty_Tx_plus1) - sqrt(P_ty_Tx))^2
        } else {
          FI_Tx <- FI_Tx + 0
        }
      }

      # Method: Numerical Derivative (ND)
      if (method == "ND") {
        # Forward difference for numerical derivative
        if (type == "group") {
          P_plus <- pmf_bn_group_imperfect(ty, Tx+1, b, barN, N, delta, lambda)
          P_current <- pmf_bn_group_imperfect(ty, Tx, b, barN, N, delta, lambda)
        } else {
          P_plus <- pmf_bn_item_imperfect(ty, Tx+1, b, barN, N, delta, lambda)
          P_current <- pmf_bn_item_imperfect(ty, Tx, b, barN, N, delta, lambda)
        }

        if (P_plus > 0 && P_current > 0) {  # Avoid division by zero
          FI_Tx <- FI_Tx + ((P_plus - P_current)^2 / P_current)
        }
      }
    }

  }

  if (method == "AD") {
    if (type == "group") {
      FI_Tx <- FI_Tx + FIpmfBN.Tx.imperfect.group(Tx, b, barN, N,delta,lambda)} else {
        FI_Tx <- FI_Tx + FIpmfBN.Tx.imperfect.item(Tx, b, barN, N,delta,lambda)
      }
  }
  return(FI_Tx)
}



#' Expected Value and Variance for Grouped Binomial (Group-level Imperfect Case)
#'
#' This function calculates the expected value and variance for a grouped binomial distribution in the imperfect case,
#' where false positive and false negative rates are considered at group-level. The function computes the expectation and variance based
#' on the number of target individuals in the population, total population size, and group-level detection rates.
#'
#' @param N Numeric. The total population size.
#' @param Tx Numeric. The number of target individuals in the population.
#' @param barN Numeric. The average number of individuals per group.
#' @param b Numeric. The total number of groups.
#' @param delta Numeric. The group-level false negative rate (probability that a positive case is classified as negative).
#' @param lambda Numeric. The group-level false positive rate (probability that a negative case is classified as positive).
#'
#' @return A list with two components:
#' \itemize{
#'   \item \code{expectation}: The expected number of detected target individuals across all groups.
#'   \item \code{variance}: The variance of the number of detected target individuals across all groups.
#' }
#' The expectation and variance are based on the given parameters, considering the detection probabilities influenced by
#' the false positive and false negative rates.
#'
#' @details
#' This function calculates the expectation and variance for the number of detected target individuals in a binomially distributed
#' scenario where imperfect detection is considered. The detection probabilities are adjusted using the false positive rate
#' (\code{lambda}) and the false negative rate (\code{delta}). The function calculates the expected number of detections
#' and the variance using the given parameters.
#'
#' @examples
#' # Example calculation for expectation and variance with group-level imperfect detection
#' ExpVarBN.imperfect.group(N = 1000, Tx = 10, barN = 50, b = 5, delta = 0.1, lambda = 0.05)
#'
#' @export
ExpVarBN.imperfect.group <- function(N, Tx, barN, b, delta, lambda) {

  prob_non_detect <- lambda * choose_custom(N-barN,Tx) / choose_custom(N,Tx) + (1 - delta) * (1 - choose_custom(N-barN,Tx) / choose_custom(N,Tx))
  prob_detect <- 1 - prob_non_detect

  # Expectation
  expectation <- b * prob_detect

  # Variance
  variance <- b * prob_detect * prob_non_detect

  return(list(expectation = expectation, variance = variance))
}



#' Expected Value and Variance for Grouped Binomial (Item-level Imperfect Case)
#'
#' This function calculates the expected value and variance for a grouped binomial distribution in the imperfect case,
#' where false positive and false negative rates are considered at item-level. The function computes the expectation and variance based
#' on the number of target individuals in the population, total population size, and group-level detection rates.
#'
#' @param N Numeric. The total population size.
#' @param Tx Numeric. The number of target individuals in the population.
#' @param barN Numeric. The average number of individuals per group.
#' @param b Numeric. The total number of groups.
#' @param delta Numeric. The group-level false negative rate (probability that a positive case is classified as negative).
#' @param lambda Numeric. The group-level false positive rate (probability that a negative case is classified as positive).
#'
#' @return A list with two components:
#' \itemize{
#'   \item \code{expectation}: The expected number of detected target individuals across all groups.
#'   \item \code{variance}: The variance of the number of detected target individuals across all groups.
#' }
#' The expectation and variance are based on the given parameters, considering the detection probabilities influenced by
#' the false positive and false negative rates.
#'
#' @details
#' This function calculates the expectation and variance for the number of detected target individuals in a binomially distributed
#' scenario where imperfect detection is considered. The detection probabilities are adjusted using the false positive rate
#' (\code{lambda}) and the false negative rate (\code{delta}). The function calculates the expected number of detections
#' and the variance using the given parameters.
#'
#' @examples
#' # Example calculation for expectation and variance with group-level imperfect detection
#' ExpVarBN.imperfect.group(N = 1000, Tx = 10, barN = 50, b = 5, delta = 0.1, lambda = 0.05)
#'
#' @export
ExpVarBN.imperfect.item <- function(N, Tx, barN, b, delta, lambda) {

  prob_non_detect <- ( lambda*(1-Tx / N) + (1-delta) * (Tx / N) )^barN
  prob_detect <- 1 - prob_non_detect

  # Expectation
  expectation <- b * prob_detect

  # Variance
  variance <- b * prob_detect * prob_non_detect

  return(list(expectation = expectation, variance = variance))
}


#' Expected Value and Variance for Grouped Binomial (Imperfect Case at either group or item level)
#'
#' This function calculates the expected value and variance for a binomial distribution in the imperfect detection case,
#' with the option to compute the results at either the group level or item level. The function routes the calculation
#' to the appropriate sub-function depending on the specified type ("group" or "item").
#'
#' @param N Numeric. The total population size.
#' @param Tx Numeric. The number of target individuals in the population.
#' @param barN Numeric. The average number of individuals per group.
#' @param b Numeric. The total number of groups.
#' @param delta Numeric. The group or item-level false negative rate (probability that a positive case is classified as negative).
#' @param lambda Numeric. The group or item-level false positive rate (probability that a negative case is classified as positive).
#' @param type Character. The type of calculation: either "group" for group-level or "item" for item-level calculation. Default is "group".
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{expectation}: The expected number of detected target individuals based on the chosen type.
#'   \item \code{variance}: The variance of the number of detected target individuals based on the chosen type.
#' }
#' The returned list contains the expectation and variance for the selected model type (group-level or item-level).
#'
#' @details
#' This unified function selects the appropriate method for calculating expectation and variance based on the \code{type} parameter.
#' If \code{type} is "group", the function calls \code{ExpVarBN.imperfect.group}, which calculates the expectation and variance at the group level.
#' If \code{type} is "item", the function calls \code{ExpVarBN.imperfect.item}, which calculates the expectation and variance at the item level.
#' The function then returns the results for the selected type of model.
#'
#' @examples
#' # Example: Expectation and variance calculation for group-level model
#' ExpVarBN.imperfect(N = 1000, Tx = 10, barN = 50, b = 5, delta = 0.1, lambda = 0.05, type = "group")
#'
#' # Example: Expectation and variance calculation for item-level model
#' ExpVarBN.imperfect(N = 1000, Tx = 10, barN = 50, b = 5, delta = 0.1, lambda = 0.05, type = "item")
#'
#' @export
ExpVarBN.imperfect <- function(N, Tx, barN, b, delta, lambda, type = c("group", "item")) {
  # Validate input type
  type <- match.arg(type)

  if (type == "group") {
    # Group-level calculation
    ExpVarHG <- ExpVarBN.imperfect.group(N, Tx, barN, b, delta, lambda)
  } else {
    # Item-level calculation
    ExpVarHG <- ExpVarBN.imperfect.item(N, Tx, barN, b, delta, lambda)
  }

  # Return calculated results
  return(list(
    expectation = ExpVarHG$expectation,
    variance = ExpVarHG$variance
  ))
}

#' Simulate Grouped Binomial Sampling with Sensitivity and Specificity Adjustments
#'
#' This function simulates a grouped binomial sampling process, accounting for imperfect sensitivity and specificity in detection.
#'
#' @param N Integer. The total population size.
#' @param barN Integer. The group size.
#' @param Tx Integer. The number of infectives in the population.
#' @param b Integer. The number of group samples.
#' @param delta Numeric (0 to 1). Sensitivity of detection (delta = 1 removes sensitivity issues).
#' @param lambda Numeric (0 to 1). Specificity of detection (lambda = 1 removes specificity issues).
#' @param NoSim Integer. The number of simulation iterations.
#' @param type Character. Either `"item"` or `"group"`, specifying the detection model.
#'
#' @return A list containing simulation results:
#'   - If `type = "item"`:
#'     - `tydat`: Matrix of total detections.
#'     - `tydatdelta`: Matrix of detections considering sensitivity.
#'     - `tydatlambda`: Matrix of detections considering specificity.
#'     - `tydatdeltalambda`: Matrix of detections considering both sensitivity and specificity.
#'   - If `type = "group"`:
#'     - `tydat`: Matrix of total detections.
#'     - `tydatDelta`: Matrix of detections considering sensitivity.
#'     - `tydatLambda`: Matrix of detections considering specificity.
#'     - `tydatDeltaLambda`: Matrix of detections considering both sensitivity and specificity.
#'
#' @examples
#' simGroupedBN(N = 1000, barN = 4, Tx = 20, b = 4, delta = 0.7, lambda = 0.8, NoSim = 100, type = "item")
#'
#' @export
simGroupedBN <- function(N,            # Population size
                         barN,           # Group size
                         Tx,            # Number of infectives in the population
                         b,              # Number of group samples
                         delta,        # Sensitivity (0 to 1, =1 removes sensitivity issues)
                         lambda,        # Specificity (0 to 1, =1 removes specificity issues)
                         NoSim,     # Number of simulations
                         type=c("item","group")) {

  tydat <- matrix(0, nrow = NoSim, ncol = b + 2)
  if(type=="item"){
    tydatdelta <- matrix(0, nrow = NoSim, ncol = b + 2)
    tydatlambda <- matrix(0, nrow = NoSim, ncol = b + 2)
    tydatdeltalambda <- matrix(0, nrow = NoSim, ncol = b + 2)
  } else {
    tydatDelta <- matrix(0, nrow = NoSim, ncol = b + 2)
    tydatLambda <- matrix(0, nrow = NoSim, ncol = b + 2)
    tydatDeltaLambda <- matrix(0, nrow = NoSim, ncol = b + 2)
  }

  set.seed(123456)


  for (j in 1:NoSim) {
    if(type=="item"){
      ty <- tydelta <- tylambda <- tydeltalambda <- 0 } else {
        ty <- tyDelta <- tyLambda <- tyDeltaLambda <- 0
      }

    # sX <- 0

    for (i in 1:b) {

      # Hypergeometric sampling for the number of infectives in the current group
      # getX <- rhyper(1, Tx - sX, N - Tx - sX, barN)
      # sX <- sX + getX

      # Binomial sampling for the number of infectives in the current group
      getX <- rbinom(1, size = barN, prob = Tx / N)

      # Sensitivity-based detection
      Y <- rbinom(getX, 1, delta)
      # Specificity-based detection
      Z <- rbinom(barN - getX, 1, 1 - lambda)


      # Total detection
      getY <- sum(Y)
      getZ <- sum(Z)
      getYZ <- getY + getZ

      if (getX > 0) ty <- ty + 1

      if(type=="item"){
        # Update counts for various detection cases
        if (getY > 0) tydelta <- tydelta + 1
        if (getZ > 0) tylambda <- tylambda + 1
        if (getYZ > 0) tydeltalambda <- tydeltalambda + 1
      } else {
        # Delta and Lambda detections
        Dv <- rbinom(1, 1, delta)
        Gv <- rbinom(1, 1, lambda)

        if (getX > 0 && Dv > 0) tyDelta <- tyDelta + 1
        if (getX == 0 && Gv == 0) tyLambda <- tyLambda + 1
        if ((getX > 0 && Dv > 0) || (getX == 0 && Gv == 0)) tyDeltaLambda <- tyDeltaLambda + 1

      }

      # Store results
      tydat[j, b + 2] <- ty

      if(type=="item"){

        tydatdelta[j, b + 2] <- tydelta
        tydatlambda[j, b + 2] <- tylambda
        tydatdeltalambda[j, b + 2] <- tydeltalambda } else {
          tydatDelta[j, b + 2] <- tyDelta
          tydatLambda[j, b + 2] <- tyLambda
          tydatDeltaLambda[j, b + 2] <- tyDeltaLambda
        }

    }

    for (i in 1:(b + 1)) {

      if (ty == (i - 1)) tydat[j, i] <- 1

      if(type=="item"){
        if (tydelta == (i - 1)) tydatdelta[j, i] <- 1
        if (tylambda == (i - 1)) tydatlambda[j, i] <- 1
        if (tydeltalambda == (i - 1)) tydatdeltalambda[j, i] <- 1}else {

          if (tyDelta == (i - 1)) tydatDelta[j, i] <- 1
          if (tyLambda == (i - 1)) tydatLambda[j, i] <- 1
          if (tyDeltaLambda == (i - 1)) tydatDeltaLambda[j, i] <- 1 }

    }
  }

  if(type=="item"){
    list(tydat=tydat,tydatdelta=tydatdelta,tydatlambda=tydatlambda,tydatdeltalambda=tydatdeltalambda)} else {
      list(tydat=tydat,tydatDelta=tydatDelta,tydatLambda=tydatLambda,tydatDeltaLambda=tydatDeltaLambda)}

}
