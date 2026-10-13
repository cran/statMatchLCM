# ============================================================
# Imports
# ============================================================

#' @importFrom stats as.formula coef fitted rmultinom setNames
NULL

# ============================================================
# Factor / numeric utilities
# ============================================================

#' Convert factor/character variables to numeric with mapping
#'
#' Converts all factor or character columns in a data.frame
#' to numeric codes and stores the mapping tables.
#'
#' @param df data.frame with factor or character columns
#' @return A list with:
#' \describe{
#'   \item{data}{data.frame with numeric-coded variables}
#'   \item{levels}{list of factor levels}
#'   \item{tables}{list of mapping tables (factor → numeric)}
#' }
#' @export
#' @examples
#' data(datA)
#' fact_to_num(datA)
fact_to_num <- function(df) {
  maps <- list()
  tables <- list()
  out <- df
  for (nm in names(df)) {
    x <- df[[nm]]
    if (is.factor(x) || is.character(x)) {
      f <- factor(x)
      out[[nm]] <- as.numeric(f)
      maps[[nm]] <- levels(f)
      tables[[nm]] <- data.frame(
        level_fact = levels(f),
        level_num  = seq_along(levels(f))
      )
    }
  }
  list(data = out, levels = maps, tables = tables)
}

#' Convert numeric codes back to factor
#'
#' Restores a factor variable from numeric codes
#' using a mapping table created by \code{fact_to_num()}.
#'
#' @param x numeric vector
#' @param table data.frame with columns level_fact and level_num
#' @return A factor vector:
#' \describe{
#'   \item{levels}{Defined by \code{table$level_num}.}
#'   \item{labels}{Defined by \code{table$level_fact}.}
#' }
#' @export
#' @examples
#' # Convert factors to numeric codes and restore the original labels
#' data(datA)
#'
#' converted <- fact_to_num(datA)
#'
#' color_restored <- num_to_fact(
#'   x = converted$data$color,
#'   table = converted$tables$color
#' )
#'
#' color_restored
#' levels(color_restored)
num_to_fact <- function(x, table) {
  factor(x, levels = table$level_num, labels = table$level_fact)
}

# ============================================================
# Data preparation
# ============================================================

#' Create stacked A/B dataset for MDC
#'
#' Creates a joint data set with missing Y1/Z1,
#' required for mass data combination.
#'
#' @param datA data.frame A
#' @param datB data.frame B
#' @param Y the target variable in datA
#' @param Z the target variable in datB
#' @return data.frame of class() "datAB_SM" with harmonized structure
#' @export
#' @examples
#' data(datA)
#' data(datB)
#' datAB_to_SM(datA, datB, color, targetB)
datAB_to_SM <- function(datA, datB, Y, Z) {

  Y_expr <- substitute(Y)
  Z_expr <- substitute(Z)

  Y <- if (is.character(Y_expr)) {
    as.character(Y_expr)
  } else if (is.name(Y_expr)) {
    as.character(Y_expr)
  } else {
    stop("`Y` must be supplied as a column name, e.g. `color` or \"color\".")
  }

  Z <- if (is.character(Z_expr)) {
    as.character(Z_expr)
  } else if (is.name(Z_expr)) {
    as.character(Z_expr)
  } else {
    stop("`Z` must be supplied as a column name, e.g. `targetB` or \"targetB\".")
  }

  if (!is.data.frame(datA) || !is.data.frame(datB)) {
    stop("`datA` and `datB` must be data frames.")
  }

  if (!Y %in% names(datA)) {
    stop("`Y` must be the name of a variable in `datA`.")
  }

  if (!Z %in% names(datB)) {
    stop("`Z` must be the name of a variable in `datB`.")
  }

  datA_full <- datA
  datA_full[[Z]] <- NA
  datA_full$source <- "datA"

  datB_full <- datB
  datB_full[[Y]] <- NA
  datB_full$source <- "datB"

  datAB <- rbind(datA_full, datB_full)
  datAB <- data.frame(lapply(datAB, factor))

  attr(datAB, "Y") <- Y
  attr(datAB, "Z") <- Z

  class(datAB) <- c("datAB_SM", "data.frame")
  datAB
}

#' Create stacked A/B/C dataset for MDA
#'
#' Creates a joint data set with missing Y1/Z1,
#' required for mass data combination.
#'
#' @param datA data.frame A
#' @param datB data.frame B
#' @param datC data.frame C
#' @param Y the target variable in datA
#' @param Z the target variable in datB
#' @return data.frame of class() "datABC_SM" with harmonized structure
#' @export
#' @examples
#' data(datA)
#' data(datB)
#' data(datC)
#' datABC_to_SM(datA, datB, datC, color, targetB)
datABC_to_SM <- function(datA, datB, datC, Y, Z) {

  Y_expr <- substitute(Y)
  Z_expr <- substitute(Z)

  Y <- if (is.character(Y_expr)) {
    as.character(Y_expr)
  } else if (is.name(Y_expr)) {
    as.character(Y_expr)
  } else {
    stop("`Y` must be supplied as a column name, e.g. `color` or \"color\".")
  }

  Z <- if (is.character(Z_expr)) {
    as.character(Z_expr)
  } else if (is.name(Z_expr)) {
    as.character(Z_expr)
  } else {
    stop("`Z` must be supplied as a column name, e.g. `targetB` or \"targetB\".")
  }

  if (!is.data.frame(datA) ||
      !is.data.frame(datB) ||
      !is.data.frame(datC)) {
    stop("`datA`, `datB`, and `datC` must be data frames.")
  }

  if (Y == Z) {
    stop("`Y` and `Z` must refer to different variables.")
  }

  if (!Y %in% names(datA)) {
    stop("`Y` must be the name of a variable in `datA`.")
  }

  if (!Z %in% names(datB)) {
    stop("`Z` must be the name of a variable in `datB`.")
  }

  if (!all(c(Y, Z) %in% names(datC))) {
    stop("`datC` must contain both target variables.")
  }

  datA_full <- datA
  datA_full[[Z]] <- NA
  datA_full$source <- "datA"

  datB_full <- datB
  datB_full[[Y]] <- NA
  datB_full$source <- "datB"

  datC_full <- datC
  datC_full$source <- "datC"

  datABC <- rbind(datA_full, datB_full, datC_full)
  datABC <- data.frame(lapply(datABC, factor))

  attr(datABC, "Y") <- Y
  attr(datABC, "Z") <- Z

  class(datABC) <- c("datABC_SM", "data.frame")
  datABC
}

# ============================================================
# SMC – Step 1 (common for all SMC methods)
# ============================================================

#' SMC step 1: selection of best imputed dataset
#'
#' Identifies and selects the imputed dataset that
#' minimizes the Hellinger distance between the reference and synthetic distributions,
#' thereby achieving the highest level of distributional similarity and
#' improving the statistical consistency of the imputation.
#'
#'
#' @param datAB data.frame of the class datAB_SM return by the function \code{datAB_to_SM}
#' @param output_ll list with imputed datasets (impdata)
#' @return A list with imputed datasets:
#' \describe{
#'   \item{datAB_imp}{The full imputed dataset obtained after combining A and B.}
#'   \item{datA_imp}{Subset of \code{datAB_imp} corresponding to dataset A.}
#'   \item{datB_imp}{Subset of \code{datAB_imp} corresponding to dataset B.}
#' }
#' @export
#' @references
#' \itemize{
#'   \item Lewaa, I., Hafez, M. S., Ismail, M. A. (2023).
#'   \emph{Mixed Statistical Matching Approaches Using a Latent Class Model: Simulation Studies}.
#'   Journal of Statistics Applications  Probability, 12(1), 247--265.
#'
#'   \item D'Orazio, M., Di Zio, M., and Scanu, M. (2006).
#'   \emph{Statistical Matching: Theory and Practice}.
#'   John Wiley and Sons.
#'
#'   \item D'Orazio, M., Di~Zio, M., and Scanu, M. (2019).
#'   Auxiliary variable selection in a statistical matching problem.
#'
#'   \item Zhang, L.-C., and Chambers, R.
#'   Analysis of integrated data.
#'   CRC/Chapman and Hall, pp. 101--120.
#'
#'   \item Conti, P. L., Marella, D., and Scanu, M. (2016).
#'   Statistical matching analysis for complex survey data with applications.
#'   \emph{Journal of the American Statistical Association}, 111(516), 1715--1725.
#' }
#' @examples
#' if (requireNamespace("NPBayesImputeCat", quietly = TRUE)) {
#'   data(datA)
#'   data(datB)
#'
#'   datAB <- datAB_to_SM(datA, datB, "color", "targetB")
#'
#'   # call DPMPM (reduced settings for speed)
#'   output_list <- NPBayesImputeCat::DPMPM_nozeros_imp(
#'     X = datAB,
#'     nrun = 50,
#'     burn = 10,
#'     thin = 10,
#'     K = 20,
#'     aalpha = 0.25,
#'     balpha = 0.25,
#'     m = 2,
#'     seed = 1234,
#'     silent = TRUE
#'   )
#'
#'   step_first <- smc_step1(datAB, output_list)
#'   str(step_first$datA_imp)
#'   str(step_first$datB_imp)
#' }
smc_step1 <- function(datAB, output_ll){
  if (!inherits(datAB, "datAB_SM")) {
    stop("`datAB` must be an object of class `datAB_SM`.", call. = FALSE)
  }

  Y <- attr(datAB, "Y", exact = TRUE)
  Z <- attr(datAB, "Z", exact = TRUE)

  if (is.null(Y) || is.null(Z)) {
    stop("`datAB` must contain the `Y` and `Z` attributes.", call. = FALSE)
  }

  add_target_attributes <- function(x, Y, Z) {
    attr(x, "Y") <- Y
    attr(x, "Z") <- Z
    x
  }

  datAB_imp_ll <- lapply(
    output_ll$impdata,
    add_target_attributes,
    Y = Y,
    Z = Z
  )

  list1 <- lapply(datAB_imp_ll, function(x) {
    datA_imp <- x[x$source == "datA", , drop = FALSE]
    datB <- x[x$source == "datB", , drop = FALSE]
    prob_observed <- table(datB[[Z]])
    prob_synthetic <- table(factor(datA_imp[[Z]], levels = levels(datB[[Z]])))
    Hell_meas <- comp.prop(p1 = prob_synthetic, p2 = prob_observed, n1 = nrow(datA_imp), n2 = NULL, ref = TRUE)$meas
    list(Hell_meas = Hell_meas, datA_imp = datA_imp)
  })

  df1 <- data.frame(do.call(rbind, lapply(list1, `[[`, "Hell_meas")))
  minHell <- which.min(df1$Hell)
  datAB_imp <- datAB_imp_ll[[minHell]]
  rownames(datAB_imp) <- NULL

  out <- list(datAB_imp = datAB_imp, datA_imp = datAB_imp[datAB_imp$source == "datA", , drop = FALSE], datB_imp = datAB_imp[datAB_imp$source == "datB", , drop = FALSE])
  attr(out, "Y") <- Y
  attr(out, "Z") <- Z
  class(out) <- c("step1_SMC", "list")
  return(out)
}

# ============================================================
# SMC – Step 2 methods
# ============================================================

#' SMC1 - Nearest-neighbour hot deck on original variables
#'
#' Performs nearest-neighbour hot deck imputation on the selected imputed datasets by matching recipient units in dataset
#' A with donor units in dataset B based on common variables.
#' The procedure transfers the target variable from the nearest donor to construct a statistically coherent fused dataset.
#'
#' @param step1 output from \code{smc_step1()}
#' @importFrom StatMatch NND.hotdeck
#' @return A list with:
#' \describe{
#'   \item{datA_fused_2}{The final fused dataset A after step 2 of the procedure.}
#'   \item{out_nndd}{A list containing nearest-neighbour distance matching results.}
#'   \item{match_vars}{Matching variables in the nearest-neighbour distance matching procedure.}
#' }
#' @export
#' @references
#' \itemize{
#'   \item Lewaa, I., Hafez, M. S., and Ismail, M. A. (2023).
#'   \emph{Mixed Statistical Matching Approaches Using a Latent Class Model: Simulation Studies}.
#'   Journal of Statistics Applications and Probability, 12(1), 247--265.
#'
#'   \item D'Orazio, M., Di Zio, M., and Scanu, M. (2006).
#'   \emph{Statistical Matching: Theory and Practice}.
#'   John Wiley and Sons.
#'
#'   \item D'Orazio, M., Di~Zio, M., and Scanu, M. (2019).
#'   Auxiliary variable selection in a statistical matching problem.
#'
#'   \item Zhang, L.-C., and Chambers, R.
#'   Analysis of integrated data.
#'   CRC/Chapman and Hall, pp. 101--120.
#'
#'   \item Conti, P. L., Marella, D., and Scanu, M. (2016).
#'   Statistical matching analysis for complex survey data with applications.
#'   \emph{Journal of the American Statistical Association}, 111(516), 1715--1725.
#' }
#' @examples
#' if (requireNamespace("NPBayesImputeCat", quietly = TRUE)) {
#'   data(datA)
#'   data(datB)
#'
#'   datAB <- datAB_to_SM(datA, datB, "color", "targetB")
#'
#'   # call DPMPM (reduced settings for speed)
#'   output_list <- NPBayesImputeCat::DPMPM_nozeros_imp(
#'     X = datAB,
#'     nrun = 50,
#'     burn = 10,
#'     thin = 10,
#'     K = 20,
#'     aalpha = 0.25,
#'     balpha = 0.25,
#'     m = 2,
#'     seed = 1234,
#'     silent = TRUE
#'   )
#'
#'   step_first <- smc_step1(datAB, output_list)
#'   set.seed(1234)
#'   step_second <- smc1_step2(step_first)
#'
#'   result <- step_second$datA_fused_2
#'   result
#' }
smc1_step2 <- function(step1){
  if (!inherits(step1, "step1_SMC")) {
    stop("`step1` must be an object of class `step1_SMC`.", call. = FALSE)
  }
  datA_imp <- step1$datA_imp
  datB_imp <- step1$datB_imp
  Y <- attr(datB_imp, "Y", exact = TRUE)
  Z <- attr(datB_imp, "Z", exact = TRUE)
  datA_imp$source <- NULL; datB_imp$source <- NULL

  match.vars = setdiff(names(datA_imp), Z)
  out_nndd <- NND.hotdeck(
    data.rec = datA_imp,
    data.don = datB_imp,
    dist.fun="Gower",
    match.vars = match.vars
  )

  idxA <- out_nndd$mtc.ids[, 2]
  Z_map <- setNames(datB_imp[[Z]], rownames(datB_imp))
  datA_imp[[Z]] <- factor(Z_map[idxA], levels = levels(datB_imp[[Z]]))

  out <- list(datA_fused_2 = datA_imp, out_nndd = out_nndd, match_vars = match.vars)
  attr(out, "Y") <- Y
  attr(out, "Z") <- Z

  class(out) <- c("step2_SMC", "list")
  out
}

#' SMC2 – Hot deck on fitted multinomial probabilities
#'
#' Implements a model-based hot deck data fusion approach by estimating multinomial models on both datasets
#' and matching observations based on their fitted probability distributions
#'
#' @param step1 output from \code{smc_step1()}
#' @importFrom StatMatch NND.hotdeck
#' @return A list with:
#' \describe{
#'   \item{datA_fused_2}{A data frame containing the final fused (imputed) version of dataset A.}
#'   \item{out_nndd}{A list with results from the nearest-neighbour distance matching procedure.}
#'   \item{match_vars}{Matching variables in the nearest-neighbour distance matching procedure.}
#'   \item{model_form}{A formula object used to fit the models.}
#'   \item{fitA}{A data frame of fitted values obtained from the model for dataset A.}
#'   \item{fitB}{A data frame of fitted values obtained from the model for dataset B.}
#' }
#' @export
#' @references
#' \itemize{
#'   \item Lewaa, I., Hafez, M. S., and Ismail, M. A. (2023).
#'   \emph{Mixed Statistical Matching Approaches Using a Latent Class Model: Simulation Studies}.
#'   Journal of Statistics Applications and Probability, 12(1), 247--265.
#'
#'   \item D'Orazio, M., Di Zio, M., and Scanu, M. (2006).
#'   \emph{Statistical Matching: Theory and Practice}.
#'   John Wiley and Sons.
#'
#'   \item D'Orazio, M., Di~Zio, M., and Scanu, M. (2019).
#'   Auxiliary variable selection in a statistical matching problem.
#'
#'   \item Zhang, L.-C., and Chambers, R.
#'   Analysis of integrated data.
#'   CRC/Chapman and Hall, pp. 101--120.
#'
#'   \item Conti, P. L., Marella, D., and Scanu, M. (2016).
#'   Statistical matching analysis for complex survey data with applications.
#'   \emph{Journal of the American Statistical Association}, 111(516), 1715--1725.
#' }
#' @examples
#' if (requireNamespace("NPBayesImputeCat", quietly = TRUE)) {
#'   data(datA)
#'   data(datB)
#'
#'   datAB <- datAB_to_SM(datA, datB, "color", "targetB")
#'
#'   # call DPMPM (reduced settings for speed)
#'   output_list <- NPBayesImputeCat::DPMPM_nozeros_imp(
#'     X = datAB,
#'     nrun = 50,
#'     burn = 10,
#'     thin = 10,
#'     K = 20,
#'     aalpha = 0.25,
#'     balpha = 0.25,
#'     m = 2,
#'     seed = 1234,
#'     silent = TRUE
#'   )
#'
#'   step_first <- smc_step1(datAB, output_list)
#'   set.seed(1234)
#'   step_second <- smc2_step2(step_first)
#'
#'   result <- step_second$datA_fused_2
#'   result
#' }
smc2_step2 <- function(step1){
  if (!inherits(step1, "step1_SMC")) {
    stop("`step1` must be an object of class `step1_SMC`.", call. = FALSE)
  }

  datA_imp <- step1$datA_imp
  datB_imp <- step1$datB_imp
  Y <- attr(datA_imp, "Y", exact = TRUE)
  Z <- attr(datB_imp, "Z", exact = TRUE)

  datA_imp$source <- NULL; datB_imp$source <- NULL

  predictors <- setdiff(names(datA_imp), Z)
  formula_nnet <- as.formula(paste(Z, "~", paste(predictors, collapse = " + ")))
  modA <- multinom(formula_nnet, data = datA_imp, trace = FALSE)
  modB <- multinom(formula_nnet, data = datB_imp, trace = FALSE)

  fitted_colnames <- function(mod) {
    levs <- mod$lev
    fit <- as.data.frame(fitted(mod))
    if (length(levs) == 2L) {
      p2 <- fit[[1L]]
      fit <- data.frame(1 - p2, p2, check.names = FALSE, row.names = row.names(fit))
    }
    names(fit) <- levs
    fit
  }
  fitA <- fitted_colnames(modA)
  fitB <- fitted_colnames(modB)

  match.vars = colnames(fitA)
  out_nndd <- NND.hotdeck(
    data.rec = cbind(datA_imp, fitA),
    data.don = cbind(datB_imp, fitB),
    dist.fun = "Gower",
    match.vars = match.vars
  )

  idxA <- out_nndd$mtc.ids[, 2]
  Z_map <- setNames(datB_imp[[Z]], rownames(datB_imp))
  datA_imp[[Z]] <- factor(Z_map[idxA], levels = levels(datB_imp[[Z]]))

  out <- list(
    datA_fused_2 = datA_imp,
    out_nndd = out_nndd,
    match_vars = match.vars,
    model_form = formula_nnet,
    fitA = fitA,
    fitB = fitB
   )
  attr(out, "Y") <- Y
  attr(out, "Z") <- Z

  class(out) <- c("step2_SMC", "list")
  out
}

#' SMC3 – Multinomial simulation approach
#'
#' Generates the target variable using multinomial simulation
#' based on fitted probability distributions to preserve uncertainty and variability
#'
#' @param step1 output from \code{smc_step1()}
#' @return A list with:
#' \describe{
#'   \item{datA_fused_2}{A data frame containing the final fused (imputed) version of dataset A.}
#'   \item{fitA}{A data frame of fitted values obtained from the model for dataset A.}
#'   \item{model_form}{A formula object used to fit the models.}
#'   \item{modA}{A multinomial regression model on datA_imp.}
#' }
#' @export
#' @references
#' \itemize{
#'   \item Lewaa, I., Hafez, M. S., and Ismail, M. A. (2023).
#'   \emph{Mixed Statistical Matching Approaches Using a Latent Class Model: Simulation Studies}.
#'   Journal of Statistics Applications and Probability, 12(1), 247--265.
#'
#'   \item D'Orazio, M., Di Zio, M., and Scanu, M. (2006).
#'   \emph{Statistical Matching: Theory and Practice}.
#'   John Wiley and Sons.
#'
#'   \item D'Orazio, M., Di~Zio, M., and Scanu, M. (2019).
#'   Auxiliary variable selection in a statistical matching problem.
#'
#'   \item Zhang, L.-C., and Chambers, R.
#'   Analysis of integrated data.
#'   CRC/Chapman and Hall, pp. 101--120.
#'
#'   \item Conti, P. L., Marella, D., and Scanu, M. (2016).
#'   Statistical matching analysis for complex survey data with applications.
#'   \emph{Journal of the American Statistical Association}, 111(516), 1715--1725.
#' }
#' @examples
#' if (requireNamespace("NPBayesImputeCat", quietly = TRUE)) {
#'   data(datA)
#'   data(datB)
#'
#'   datAB <- datAB_to_SM(datA, datB, "color", "targetB")
#'
#'   # call DPMPM (reduced settings for speed)
#'   output_list <- NPBayesImputeCat::DPMPM_nozeros_imp(
#'     X = datAB,
#'     nrun = 50,
#'     burn = 10,
#'     thin = 10,
#'     K = 20,
#'     aalpha = 0.25,
#'     balpha = 0.25,
#'     m = 2,
#'     seed = 1234,
#'     silent = TRUE
#'   )
#'
#'   step_first <- smc_step1(datAB, output_list)
#'   set.seed(1234)
#'   step_second <- smc3_step2(step_first)
#'
#'   result <- step_second$datA_fused_2
#'   result
#' }
smc3_step2 <- function(step1){
  if (!inherits(step1, "step1_SMC")) {
    stop("`step1` must be an object of class `step1_SMC`.", call. = FALSE)
  }

  datA_imp <- step1$datA_imp
  datA_imp$source <- NULL
  Y <- attr(datA_imp, "Y", exact = TRUE)
  Z <- attr(datA_imp, "Z", exact = TRUE)

  predictors <- setdiff(names(datA_imp), Z)
  formula_nnet <- as.formula(paste(Z, "~", paste(predictors, collapse = " + ")))
  modA <- multinom(formula_nnet, data = datA_imp, trace = FALSE)

  fitted_colnames <- function(mod) {
    levs <- mod$lev
    fit <- as.data.frame(fitted(mod))
    if (length(levs) == 2L) {
      p2 <- fit[[1L]]
      fit <- data.frame(1 - p2, p2, check.names = FALSE, row.names = row.names(fit))
    }
    names(fit) <- levs
    fit
  }
  fitA <- fitted_colnames(modA)
  Z_new <- vapply(seq_len(nrow(fitA)), function(i) {
    indeks <- which(rmultinom(1, 1, prob = fitA[i, ]) == 1)
    names(fitA)[indeks]
  }, character(1)) |> as.factor()

  datA_imp[[Z]] <- Z_new

  out <- list(datA_fused_2 = datA_imp, fitA = fitA, model_form = formula_nnet, modA = modA)
  attr(out, "Y") <- Y
  attr(out, "Z") <- Z

  class(out) <- c("step2_SMC", "list")
  out
}


# ============================================================
# SMA – Step 1 (common for all SMA methods)
# ============================================================

#' SMA step 1: selection of best imputed dataset
#'
#' Selects the imputed dataset minimizing the Hellinger distance
#' between the reference distribution of Z from dataset B and the
#' synthetic distribution of Z from dataset A, in a three-sample
#' statistical matching framework (A, B, C).
#'
#' @param datABC data.frame of class datABC_SM returned by \code{datABC_to_SM}
#' @param output_ll list with imputed datasets (impdata)
#' @return A list with imputed datasets:
#' \describe{
#'   \item{datABC_imp}{The full imputed dataset combining A, B, and C.}
#'   \item{datA_imp}{Subset of the imputed data corresponding to dataset A.}
#'   \item{datB_imp}{Subset corresponding to dataset B.}
#'  }
#'
#' @importFrom StatMatch comp.prop
#' @export
#' @references
#' \itemize{
#'   \item Lewaa, I., Hafez, M. S., and Ismail, M. A. (2023).
#'   \emph{Mixed Statistical Matching Approaches Using a Latent Class Model: Simulation Studies}.
#'   Journal of Statistics Applications and Probability, 12(1), 247--265.
#'
#'   \item D'Orazio, M., Di Zio, M., and Scanu, M. (2006).
#'   \emph{Statistical Matching: Theory and Practice}.
#'   John Wiley and Sons.
#'
#'   \item D'Orazio, M., Di~Zio, M., and Scanu, M. (2019).
#'   Auxiliary variable selection in a statistical matching problem.
#'
#'   \item Zhang, L.-C., and Chambers, R.
#'   Analysis of integrated data.
#'   CRC/Chapman and Hall, pp. 101--120.
#'
#'   \item Conti, P. L., Marella, D., and Scanu, M. (2016).
#'   Statistical matching analysis for complex survey data with applications.
#'   \emph{Journal of the American Statistical Association}, 111(516), 1715--1725.
#' }
#' @examples
#' if (requireNamespace("NPBayesImputeCat", quietly = TRUE)) {
#'   data(datA)
#'   data(datB)
#'   data(datC)
#'
#'   datABC <- datABC_to_SM(datA, datB, datC, "color", "targetB")
#'
#'   # call DPMPM
#'   output_list_AII <- NPBayesImputeCat::DPMPM_nozeros_imp(
#'     X = datABC, nrun = 500, burn = 50, thin = 50,
#'     K = 80, aalpha = 0.25, balpha = 0.25,
#'     m = 2, seed = 1234, silent = FALSE
#'   )
#'
#'   step_first_aii <- sma_step1(datABC, output_list_AII)
#'   names(step_first_aii)
#' }
sma_step1 <- function(datABC, output_ll) {

  if (!inherits(datABC, "datABC_SM")) {
    stop("`datABC` must be an object of class `datABC_SM`.",
         call. = FALSE)
  }

  Y <- attr(datABC, "Y", exact = TRUE)
  Z <- attr(datABC, "Z", exact = TRUE)

  if (!all(c(Y, Z, "source") %in% names(datABC))) {
    stop("Variables specified by `Y`, `Z`, and `source` must occur in `datABC`.",
         call. = FALSE)
  }

  if (is.null(Y) || is.null(Z)) {
    stop("`datABC` must contain the `Y` and `Z` attributes.",
         call. = FALSE)
  }

  add_target_attributes <- function(x, Y, Z) {
    attr(x, "Y") <- Y
    attr(x, "Z") <- Z
    x
  }

  datABC_imp_ll <- lapply(
    output_ll$impdata,
    add_target_attributes,
    Y = Y,
    Z = Z
  )

  list1 <- lapply(datABC_imp_ll, function(x) {
    datA_imp <- x[x$source == "datA", , drop = FALSE]
    datB <- x[x$source == "datB", , drop = FALSE]

    prob_observed <- table(datB[[Z]])
    prob_synthetic <- table(
      factor(datA_imp[[Z]], levels = levels(datB[[Z]]))
    )

    Hell_meas <- comp.prop(
      p1 = prob_synthetic,
      p2 = prob_observed,
      n1 = nrow(datA_imp),
      n2 = NULL,
      ref = TRUE
    )$meas

    list(
      Hell_meas = Hell_meas,
      datA_imp = datA_imp
      )
  })

  df1 <- data.frame(
    do.call(rbind, lapply(list1, `[[`, "Hell_meas"))
  )

  minHell <- which.min(df1$Hell)

  datABC_imp <- datABC_imp_ll[[minHell]]
  rownames(datABC_imp) <- NULL

  out <- list(
    datABC_imp = datABC_imp,
    datA_imp = datABC_imp[datABC_imp$source == "datA", , drop = FALSE],
    datB_imp = datABC_imp[datABC_imp$source == "datB", , drop = FALSE]
    )

  attr(out, "Y") <- Y
  attr(out, "Z") <- Z

  class(out) <- c("step1_SMA", "list")
  out
}

#' SMA1 – Nearest-neighbour hot deck on shared variables
#'
#' Second step of the SMA procedure using nearest-neighbour hot deck matching on shared variables (Y1, Z1)
#' to fuse datasets while preserving their joint distribution.
#'
#' @param step1 list returned by the SMA step 1 procedure
#' @return A list with:
#' \describe{
#'   \item{datA_fused_2}{The final fused dataset A.}
#'   \item{out_nndd}{NNDD matching results.}
#'   \item{match.vars}{Matching variables in the nearest-neighbour distance matching procedure.}
#' }
#' @importFrom StatMatch NND.hotdeck
#' @export
#' @references
#' \itemize{
#'   \item Lewaa, I., Hafez, M. S., Ismail, M. A. (2023).
#'   \emph{Mixed Statistical Matching Approaches Using a Latent Class Model: Simulation Studies}.
#'   Journal of Statistics Applications  Probability, 12(1), 247--265.
#'
#'   \item D'Orazio, M., Di Zio, M., and Scanu, M. (2006).
#'   \emph{Statistical Matching: Theory and Practice}.
#'   John Wiley and Sons.
#'
#'   \item D'Orazio, M., Di~Zio, M., and Scanu, M. (2019).
#'   Auxiliary variable selection in a statistical matching problem.
#'
#'   \item Zhang, L.-C., and Chambers, R.
#'   Analysis of integrated data.
#'   CRC/Chapman and Hall, pp. 101--120.
#'
#'   \item Conti, P. L., Marella, D., and Scanu, M. (2016).
#'   Statistical matching analysis for complex survey data with applications.
#'   \emph{Journal of the American Statistical Association}, 111(516), 1715--1725.
#' }
#' @examples
#' if (requireNamespace("NPBayesImputeCat", quietly = TRUE)) {
#'   data(datA)
#'   data(datB)
#'
#'   data(datC)
#'
#'   datABC <- datABC_to_SM(datA, datB, datC, "color", "targetB")
#'
#'   # call DPMPM (reduced settings for speed)
#'   output_list_AII <- NPBayesImputeCat::DPMPM_nozeros_imp(
#'     X = datABC,
#'     nrun = 50,
#'     burn = 10,
#'     thin = 10,
#'     K = 20,
#'     aalpha = 0.25,
#'     balpha = 0.25,
#'     m = 2,
#'     seed = 1234,
#'     silent = TRUE
#'   )
#'
#'   step_first_aii <- sma_step1(datABC, output_list_AII)
#'   set.seed(1234)
#'   step_second_aii1 <- sma1_step2(step_first_aii)
#'
#'   result_aii1 <- step_second_aii1$datA_fused_2
#'   head(result_aii1)
#' }
sma1_step2 <- function(step1){
  if (!inherits(step1, "step1_SMA")) {
    stop("`step1` must be an object of class `step1_SMA`.", call. = FALSE)
  }

  datA_imp <- step1$datA_imp
  Y <- attr(datA_imp, "Y", exact = TRUE)
  Z <- attr(datA_imp, "Z", exact = TRUE)
  datC <- step1$datABC_imp[step1$datABC_imp$source == "datC", , drop = FALSE]

  datA_imp$source <- NULL
  match.vars = c(Y, Z)
  out_nndd <- NND.hotdeck(
    data.rec  = datA_imp,
    data.don  = datC,
    dist.fun = "Gower",
    match.vars = match.vars
  )

  idxA <- out_nndd$mtc.ids[, 2]
  Z_map <- setNames(datC[[Z]], rownames(datC))
  datA_imp[[Z]] <- factor(Z_map[idxA], levels = levels(datA_imp[[Z]]))

  out <- list(datA_fused_2 = datA_imp, out_nndd = out_nndd, match.vars = match.vars)
  attr(out, "Y") <- Y
  attr(out, "Z") <- Z

  class(out) <- c("step2_SMA", "list")
  out
}

#' SMA2 – Hot deck on fitted multinomial probabilities
#'
#' Second step of SMA using multinomial models and nearest-neighbour hot deck on fitted probabilities
#' to improve data fusion accuracy.
#'
#' @param step1 list returned by the SMA step 1 procedure
#' @return A list with:
#' \describe{
#'   \item{datA_fused_2}{A data frame containing the final fused (imputed) version of dataset A.}
#'   \item{out_nnd}{A list with results from the nearest-neighbour distance matching procedure.}
#'   \item{match_vars}{Matching variables in the nearest-neighbour distance matching procedure.}
#'   \item{modA}{A multinomial logistic model obtained for dataset datA_imp.}
#'   \item{modC}{A multinomial logistic model obtained for dataset datC.}
#' }
#' @importFrom StatMatch NND.hotdeck
#' @importFrom nnet multinom
#' @export
#' @references
#' \itemize{
#'   \item Lewaa, I., Hafez, M. S., Ismail, M. A. (2023).
#'   \emph{Mixed Statistical Matching Approaches Using a Latent Class Model: Simulation Studies}.
#'   Journal of Statistics Applications  Probability, 12(1), 247--265.
#'
#'   \item D'Orazio, M., Di Zio, M., and Scanu, M. (2006).
#'   \emph{Statistical Matching: Theory and Practice}.
#'   John Wiley and Sons.
#'
#'   \item D'Orazio, M., Di~Zio, M., and Scanu, M. (2019).
#'   Auxiliary variable selection in a statistical matching problem.
#'
#'   \item Zhang, L.-C., and Chambers, R.
#'   Analysis of integrated data.
#'   CRC/Chapman and Hall, pp. 101--120.
#'
#'   \item Conti, P. L., Marella, D., and Scanu, M. (2016).
#'   Statistical matching analysis for complex survey data with applications.
#'   \emph{Journal of the American Statistical Association}, 111(516), 1715--1725.
#' }
#' @examples
#' if (requireNamespace("NPBayesImputeCat", quietly = TRUE)) {
#'   data(datA)
#'   data(datB)
#'   data(datC)
#'
#'   datABC <- datABC_to_SM(datA, datB, datC, "color", "targetB")
#'
#'     output_list_AII <- NPBayesImputeCat::DPMPM_nozeros_imp(
#'     X = datABC,
#'     nrun = 50,
#'     burn = 10,
#'     thin = 10,
#'     K = 20,
#'     aalpha = 0.25,
#'     balpha = 0.25,
#'     m = 2,
#'     seed = 1234,
#'     silent = TRUE
#'   )
#'
#'   step_first_aii <- sma_step1(datABC, output_list_AII)
#'   set.seed(1234)
#'   step_second_aii2 <- sma2_step2(step_first_aii)
#'
#'   result_aii2 <- step_second_aii2$datA_fused_2
#'   head(result_aii2)
#' }
sma2_step2 <- function(step1){
  if (!inherits(step1, "step1_SMA")) {
    stop("`step1` must be an object of class `step1_SMA`.", call. = FALSE)
  }

  datA_imp <- step1$datA_imp
  Y <- attr(datA_imp, "Y", exact = TRUE)
  Z <- attr(datA_imp, "Z", exact = TRUE)
  datC <- step1$datABC_imp[step1$datABC_imp$source == "datC", , drop = FALSE]
  datA_imp$source <- NULL; datC$source <- NULL

  predictors <- setdiff(names(datA_imp), Z)
  formula_nnet <- as.formula(paste(Z, "~", paste(predictors, collapse = " + ")))

  modA <- multinom(formula_nnet, data = datA_imp, trace = FALSE)
  modC <- multinom(formula_nnet, data = datC, trace = FALSE)

  fitted_colnames <- function(mod) {
    levs <- mod$lev
    fit <- as.data.frame(fitted(mod))

    if (length(levs) == 2L) {
      p2 <- fit[[1L]]

      fit <- data.frame(
        1 - p2,
        p2,
        check.names = FALSE,
        row.names = row.names(fit)
      )
    }

    names(fit) <- levs
    fit
  }

  fitA <- fitted_colnames(modA)
  fitC <- fitted_colnames(modC)

  match.vars = colnames(fitA)
  out_nnd <- NND.hotdeck(
    data.rec  = cbind(datA_imp, fitA),
    data.don  = cbind(datC, fitC),
    dist.fun = "Gower",
    match.vars = match.vars
  )

  idxA <- out_nnd$mtc.ids[, 2]
  Z_map <- setNames(datC[[Z]], rownames(datC))
  datA_imp[[Z]] <- factor(Z_map[idxA], levels = levels(datA_imp[[Z]]))

  out <- list(
    datA_fused_2 = datA_imp,
    out_nnd = out_nnd,
    match_vars = match.vars,
    modA = modA,
    modC = modC
  )
  attr(out, "Y") <- Y
  attr(out, "Z") <- Z

  class(out) <- c("step2_SMA", "list")
  out
}

#' SMA3 – Multinomial simulation approach
#'
#'Second step of SMA is to generate the target variable using multinomial probability distributions derived from a donor-based model
#'to achieve statistically coherent data fusion
#'
#' @param step1 list returned by the SMA step 1 procedure
#' @return A list with:
#' \describe{
#'   \item{datA_fused_2}{A data frame containing the final fused version of dataset A.}
#'   \item{modB}{A multinomial logistic model obtained for dataset datB_imp.}
#'   \item{zero_one}{A binary matrix (one-hot encoded) sampled from \code{prob_Z_all}, with categories in rows and observations in columns; each column contains a single 1 indicating the selected category of Z.}
#' }
#' @importFrom nnet multinom
#' @export
#' @references
#' \itemize{
#'   \item Lewaa, I., Hafez, M. S., Ismail, M. A. (2023).
#'   \emph{Mixed Statistical Matching Approaches Using a Latent Class Model: Simulation Studies}.
#'   Journal of Statistics Applications  Probability, 12(1), 247--265.
#'
#'   \item D'Orazio, M., Di Zio, M., and Scanu, M. (2006).
#'   \emph{Statistical Matching: Theory and Practice}.
#'   John Wiley and Sons.
#'
#'   \item D'Orazio, M., Di~Zio, M., and Scanu, M. (2019).
#'   Auxiliary variable selection in a statistical matching problem.
#'
#'   \item Zhang, L.-C., and Chambers, R.
#'   Analysis of integrated data.
#'   CRC/Chapman and Hall, pp. 101--120.
#'
#'   \item Conti, P. L., Marella, D., and Scanu, M. (2016).
#'   Statistical matching analysis for complex survey data with applications.
#'   \emph{Journal of the American Statistical Association}, 111(516), 1715--1725.
#' }
#' @examples
#' if (requireNamespace("NPBayesImputeCat", quietly = TRUE)) {
#'   data(datA)
#'   data(datB)
#'   data(datC)
#'
#'   datABC <- datABC_to_SM(datA, datB, datC, "color", "targetB")
#'
#'   # call DPMPM (reduced settings for speed)
#'   output_list_AII <- NPBayesImputeCat::DPMPM_nozeros_imp(
#'     X = datABC,
#'     nrun = 50,
#'     burn = 10,
#'     thin = 10,
#'     K = 20,
#'     aalpha = 0.25,
#'     balpha = 0.25,
#'     m = 2,
#'     seed = 1234,
#'     silent = TRUE
#'   )
#'
#'   step_first_aii <- sma_step1(datABC, output_list_AII)
#'   set.seed(1234)
#'   step_second_aii3 <- sma3_step2(step_first_aii)
#'
#'   result_aii3 <- step_second_aii3$datA_fused_2
#'   head(result_aii3)
#' }
sma3_step2 <- function(step1){
  datA_imp <- step1$datA_imp
  datB_imp <- step1$datB_imp
  datA_imp$source <- NULL; datB_imp$source <- NULL
  Y <- attr(datB_imp, "Y", exact = TRUE)
  Z <- attr(datB_imp, "Z", exact = TRUE)

  predictors <- setdiff(names(datB_imp), Z)
  formula_nnet <- as.formula(paste(Z, "~", paste(predictors, collapse = " + ")))
  modB <- multinom(formula_nnet, data = datB_imp, trace = FALSE)

  predict_prob_matrix <- function(mod, newdata) {
    prob <- predict(mod, newdata = newdata, type = "probs")
    levs <- mod$lev

    if (length(levs) == 2L) {
      p2 <- as.numeric(prob)
      prob <- cbind(1 - p2, p2)
    } else {
      prob <- matrix(
        prob,
        nrow = nrow(newdata),
        ncol = length(levs)
      )
    }

    colnames(prob) <- levs
    rownames(prob) <- rownames(newdata)
    prob
  }

  prob_Z_all <- predict_prob_matrix(modB, datA_imp)

  zero_one <- apply(prob_Z_all, 1, function(prob) {
    names(prob) <- colnames(prob_Z_all)
    rmultinom(1, 1, prob)
  })

  rownames(zero_one) <- colnames(prob_Z_all)
  Z_new <- apply(zero_one, 2, function(x) rownames(zero_one)[x == 1])
  Z_new <- factor(Z_new, levels = levels(datB_imp[[Z]]))
  datA_imp[[Z]] <- Z_new

  out <- list(
    datA_fused_2 = datA_imp,
    modB = modB,
    zero_one = zero_one
  )
  attr(out, "Y") <- Y
  attr(out, "Z") <- Z

  class(out) <- c("step2_SMA", "list")
  out
}

# ============================================================
# Quality assessment
# ============================================================

#' Quality assessment of synthetic Z1
#'
#' Evaluates the quality of the synthetic target variable by computing the Hellinger distance
#'  between the reference and synthetic distributions. This measure quantifies the degree of similarity
#'  between the two distributions, providing an assessment of the accuracy and coherence of the data fusion process
#'
#' @param Z_obs the variable Z_obs from datB
#' @param Z_fused the variable Z_fused from step2$datA_fused_2
#' @importFrom StatMatch comp.prop
#' @return A list with:
#' \describe{
#'   \item{heli_latent}{The complete list returned by \code{StatMatch::comp.prop()}, containing \code{meas}, \code{chi.sq}, and \code{p.exp}. The Hellinger distance is available as \code{heli_latent$meas["Hell"]}.}
#'   \item{ref_distr}{A numeric vector or table representing the reference (original) distribution.}
#'   \item{synt_distr}{A numeric vector or table representing the synthetic (generated) distribution.}
#' }
#' @export
sm_quality <- function(Z_obs, Z_fused) {
  ref_distr  <- table(Z_obs)
  synt_distr <- table(Z_fused)

  heli_latent <- comp.prop(
    p1 = synt_distr,
    p2 = ref_distr,
    n1 = length(Z_obs),
    n2 = NULL,
    ref = TRUE
  )

  list(
    heli_latent = heli_latent,
    ref_distr = ref_distr,
    synt_distr = synt_distr
  )
}

#' Dataset datA
#'
#' A simple dataset with categorical variables.
#'
#' @format A data frame with 20 observations and 2 variables:
#' \describe{
#'   \item{gender}{Category "F", "M"}
#'   \item{color}{Color category "blue", "green", "red"}
#' }
#' @source Simulated data
"datA"

#' Dataset datB
#'
#' A simple dataset with categorical variables.
#'
#' @format A data frame with 25 observations and 2 variables:
#' \describe{
#'   \item{gender}{Category "F", "M"}
#'   \item{targetB}{Letter category "d", "e", "f"}
#' }
#' @source Simulated data
"datB"

#' Dataset datC
#'
#' A simple dataset with categorical variables.
#'
#' @format A data frame with 30 observations and 3 variables:
#' \describe{
#'   \item{gender}{Category "F", "M"}
#'   \item{color}{Color category "blue", "green", "red"}
#'   \item{targetB}{Letter category "d", "e", "f"}
#' }
#' @source Simulated data
"datC"

