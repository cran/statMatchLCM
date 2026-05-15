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
#' @return data.frame with harmonized structure
#' @export
#' @examples
#' data(datA)
#' data(datB)
#' datAB_to_SM(datA, datB)
datAB_to_SM <- function(datA, datB){
  datA_full <- data.frame(datA, Z1 = NA)
  datB_full <- data.frame(datB, Y1 = NA)
  datAB <- rbind(datA_full, datB_full)
  datAB <- data.frame(lapply(datAB, function(x) factor(x)))
  return(datAB)
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
#' @param datA data.frame A
#' @param datB data.frame B
#' @param output_ll list with imputed datasets (impdata)
#' @return A list with imputed datasets:
#' \describe{
#'   \item{datAB_imp1}{The full imputed dataset obtained after combining A and B.}
#'   \item{datA_imp1}{Subset of \code{datAB_imp1} corresponding to dataset A.}
#'   \item{datB_imp1}{Subset of \code{datAB_imp1} corresponding to dataset B.}
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
#'   datAB <- datAB_to_SM(datA, datB)
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
#'   step_first <- smc_step1(datA, datB, output_list)
#'   str(step_first$datA_imp1)
#'   str(step_first$datB_imp1)
#' }
smc_step1 <- function(datA, datB, output_ll){
  splitA <- nrow(datA)
  splitB <- nrow(datA) + 1
  nrowAB <- nrow(datA) + nrow(datB)
  datAB_impll <- output_ll$impdata

  list1 <- lapply(datAB_impll, function(x) {
    datA_imp <- x[1:splitA, ]
    datB_imp <- x[splitB:nrowAB, ]
    tt0  <- table(datA_imp[, ncol(x)])
    tts1 <- table(datB_imp[, ncol(x)])
    Hell_meas <- comp.prop(
      p1 = tts1, p2 = tt0, n1 = nrow(datA),
      n2 = NULL, ref = TRUE
    )$meas
    list(Hell_meas = Hell_meas,
         datA_imp = datA_imp,
         datB_imp = datB_imp)
  })

  df1 <- data.frame(do.call(rbind, lapply(list1, `[[`, "Hell_meas")))
  minHell <- which.min(df1$Hell)

  datAB_imp1 <- output_ll$impdata[[minHell]]
  rownames(datAB_imp1) <- NULL

  list(
    datAB_imp1 = datAB_imp1,
    datA_imp1  = datAB_imp1[1:splitA, ],
    datB_imp1  = datAB_imp1[splitB:nrowAB, ]
  )
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
#'   datAB <- datAB_to_SM(datA, datB)
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
#'   step_first <- smc_step1(datA, datB, output_list)
#'   step_second <- smc1_step2(step_first)
#'
#'   result <- step_second$datA_fused_2
#'   result
#' }
smc1_step2 <- function(step1){
  datA_imp1 <- step1$datA_imp1
  datB_imp1 <- step1$datB_imp1

  datA_imp1_num <- fact_to_num(datA_imp1)$data
  datB_imp1_num <- fact_to_num(datB_imp1)$data

  X_mtcc_datA <- names(datA_imp1)[-ncol(datA_imp1)]

  out_nndd <- NND.hotdeck(
    data.rec = datA_imp1_num,
    data.don = datB_imp1_num,
    match.vars = X_mtcc_datA
  )

  idxA <- out_nndd$mtc.ids[, 2]
  Z1_map <- setNames(datB_imp1$Z1, rownames(datB_imp1))
  datA_imp1$Z1 <- factor(Z1_map[idxA], levels = levels(datB_imp1$Z1))

  list(datA_fused_2 = datA_imp1, out_nndd = out_nndd)
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
#'   \item{model_form}{A formula object used to fit the models.}
#'   \item{FitteddatA_imp1}{A data frame of fitted values obtained from the model for dataset A.}
#'   \item{FitteddatB_imp1}{A data frame of fitted values obtained from the model for dataset B.}
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
#'   datAB <- datAB_to_SM(datA, datB)
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
#'   step_first <- smc_step1(datA, datB, output_list)
#'   step_second <- smc2_step2(step_first)
#'
#'   result <- step_second$datA_fused_2
#'   result
#' }
smc2_step2 <- function(step1){
  datA_imp1 <- step1$datA_imp1
  datB_imp1 <- step1$datB_imp1

  datA_imp1_num <- fact_to_num(datA_imp1)$data
  datB_imp1_num <- fact_to_num(datB_imp1)$data

  predictors <- setdiff(names(datA_imp1_num), "Z1")
  formula_nnet <- as.formula(paste("Z1 ~", paste(predictors, collapse = " + ")))

  modA <- multinom(formula_nnet, data = datA_imp1_num)
  modB <- multinom(formula_nnet, data = datB_imp1_num)

  fitA <- as.data.frame(fitted(modA))
  fitB <- as.data.frame(fitted(modB))

  out_nndd <- NND.hotdeck(
    data.rec = cbind(datA_imp1_num, fitA),
    data.don = cbind(datB_imp1_num, fitB),
    match.vars = colnames(fitA)
  )

  idxA <- out_nndd$mtc.ids[, 2]
  Z1_map <- setNames(datB_imp1$Z1, rownames(datB_imp1))
  datA_imp1$Z1 <- factor(Z1_map[idxA], levels = levels(datB_imp1$Z1))

  list(
    datA_fused_2 = datA_imp1,
    out_nndd = out_nndd,
    model_form = formula_nnet,
    FitteddatA_imp1 = fitA,
    FitteddatB_imp1 = fitB
  )
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
#'   \item{FitteddatA_imp1}{A data frame of fitted values obtained from the model for dataset A.}
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
#'   datAB <- datAB_to_SM(datA, datB)
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
#'   step_first <- smc_step1(datA, datB, output_list)
#'   step_second <- smc3_step2(step_first)
#'
#'   result <- step_second$datA_fused_2
#'   result
#' }
smc3_step2 <- function(step1){
  datA_imp1 <- step1$datA_imp1

  datA_to_num <- fact_to_num(datA_imp1)
  datA_imp1_num <- datA_to_num$data

  predictors <- setdiff(names(datA_imp1_num), "Z1")
  formula_nnet <- as.formula(paste("Z1 ~", paste(predictors, collapse = " + ")))

  modA <- multinom(formula_nnet, data = datA_imp1_num)
  fitA <- as.data.frame(fitted(modA))

  p <- ncol(fitA)
  Z1_new <- sapply(seq_len(nrow(fitA)), function(i)
    which(rmultinom(1, 1, fitA[i, 1:p]) == 1)
  )

  datA_imp1$Z1 <- num_to_fact(Z1_new, datA_to_num$tables$Z1)

  list(datA_fused_2 = datA_imp1, FitteddatA_imp1 = fitA)
}


# ============================================================
# SMA – Step 1 (common for all SMA methods)
# ============================================================

#' SMA step 1: selection of best imputed dataset
#'
#' Selects the imputed dataset minimizing the Hellinger distance
#' between the reference distribution from dataset A and the
#' synthetic distribution from dataset B, in a three-sample
#' statistical matching framework (A, B, C).
#'
#' @param datA data.frame A
#' @param datB data.frame B
#' @param datC data.frame C
#' @param output_ll list with imputed datasets (impdata)
#' @return A list with imputed datasets:
#' \describe{
#'   \item{datABC_imp1}{The full imputed dataset combining A, B, and C.}
#'   \item{datA_imp1}{Subset of the imputed data corresponding to dataset A.}
#'   \item{datB_imp1}{Subset corresponding to dataset B.}
#'   \item{datC_imp1}{Subset corresponding to dataset C.}
#' }
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
#'
#'   datAB <- datAB_to_SM(datA, datB)
#'   datC <- data.frame(
#'     X = datB$X[1:4],
#'     Y1 = datA$Y1[1:4],
#'     Z1 = datB$Z1[1:4]
#'   )
#'
#'   # adding auxiliary information
#'   datABC <- rbind(datAB, datC)
#'
#'   # call DPMPM
#'   output_list_AII <- NPBayesImputeCat::DPMPM_nozeros_imp(
#'     X = datABC, nrun = 500, burn = 50, thin = 50,
#'     K = 80, aalpha = 0.25, balpha = 0.25,
#'     m = 2, seed = 1234, silent = FALSE
#'   )
#'
#'   step_first_aii <- sma_step1(datA, datB, datC, output_list_AII)
#'   names(step_first_aii)
#' }
sma_step1 <- function(datA, datB, datC, output_ll){
  splitA <- nrow(datA)
  splitB <- nrow(datA) + 1
  nrowAB <- nrow(datA) + nrow(datB)
  splitC <- nrowAB + 1
  nrowABC <- nrowAB + nrow(datC)

  datABC_impll <- output_ll$impdata

  list1 <- lapply(datABC_impll, function(x) {

    datA_imp <- x[1:splitA, ]
    datB_imp <- x[splitB:nrowAB, ]

    tt0  <- table(datA_imp[, ncol(x)])   # reference distribution
    tts1 <- table(datB_imp[, ncol(x)])   # synthetic distribution

    Hell_meas <- comp.prop(
      p1 = tts1,
      p2 = tt0,
      n1 = nrow(datA),
      n2 = NULL,
      ref = TRUE
    )$meas

    list(Hell_meas = Hell_meas)
  })

  df1 <- data.frame(do.call(rbind, lapply(list1, `[[`, "Hell_meas")))
  minHell <- which.min(df1$Hell)

  datABC_imp1 <- datABC_impll[[minHell]]
  rownames(datABC_imp1) <- NULL

  list(
    datABC_imp1 = datABC_imp1,
    datA_imp1   = datABC_imp1[1:splitA, ],
    datB_imp1   = datABC_imp1[splitB:nrowAB, ],
    datC_imp1   = datABC_imp1[splitC:nrowABC, ]
  )
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
#'   datAB <- datAB_to_SM(datA, datB)
#'   datC <- data.frame(
#'     X = datB$X[1:4],
#'     Y1 = datA$Y1[1:4],
#'     Z1 = datB$Z1[1:4]
#'   )
#'
#'   # adding auxiliary information
#'   datABC <- rbind(datAB, datC)
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
#'   step_first_aii <- sma_step1(datA, datB, datC, output_list_AII)
#'   step_second_aii1 <- sma1_step2(step_first_aii)
#'
#'   result_aii1 <- step_second_aii1$datA_fused_2
#'   head(result_aii1)
#' }
sma1_step2 <- function(step1){
  datA_imp1 <- step1$datA_imp1
  datC_imp1 <- step1$datC_imp1

  datA_imp1_num <- fact_to_num(datA_imp1)$data
  datC_imp1_num <- fact_to_num(datC_imp1)$data

  X_mtcc_datA <- c("Y1", "Z1")

  out_nndd <- NND.hotdeck(
    data.rec  = datA_imp1_num,
    data.don  = datC_imp1_num,
    match.vars = X_mtcc_datA
  )

  idxA <- out_nndd$mtc.ids[, 2]
  Z1_map <- setNames(datC_imp1$Z1, rownames(datC_imp1))
  datA_imp1$Z1 <- factor(Z1_map[idxA], levels = levels(datC_imp1$Z1))

  list(datA_fused_2 = datA_imp1, out_nndd = out_nndd)
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
#'   \item{model_form}{A formula object specifying the model used.}
#'   \item{FitteddatA_imp1}{A data frame of fitted values obtained for dataset A.}
#'   \item{FitteddatC_imp1}{A data frame of fitted values obtained for dataset C.}
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
#'
#'   datAB <- datAB_to_SM(datA, datB)
#'   datC <- data.frame(
#'     X = datB$X[1:4],
#'     Y1 = datA$Y1[1:4],
#'     Z1 = datB$Z1[1:4]
#'   )
#'
#'   # adding auxiliary information
#'   datABC <- rbind(datAB, datC)
#'
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
#'   step_first_aii <- sma_step1(datA, datB, datC, output_list_AII)
#'   step_second_aii2 <- sma2_step2(step_first_aii)
#'
#'   result_aii2 <- step_second_aii2$datA_fused_2
#'   head(result_aii2)
#' }
sma2_step2 <- function(step1){
  datA_imp1 <- step1$datA_imp1
  datC_imp1 <- step1$datC_imp1

  datA_imp1_num <- fact_to_num(datA_imp1)$data
  datC_imp1_num <- fact_to_num(datC_imp1)$data

  predictors <- setdiff(names(datA_imp1_num), "Z1")
  formula_nnet <- as.formula(paste("Z1 ~", paste(predictors, collapse = " + ")))

  modA <- multinom(formula_nnet, data = datA_imp1_num)
  modC <- multinom(formula_nnet, data = datC_imp1_num)

  fitA <- as.data.frame(fitted(modA))
  fitC <- as.data.frame(fitted(modC))

  out_nnd <- NND.hotdeck(
    data.rec  = cbind(datA_imp1_num, fitA),
    data.don  = cbind(datC_imp1_num, fitC),
    match.vars = colnames(fitA),
    dist.fun = "Gower"
  )

  idxA <- out_nnd$mtc.ids[, 2]
  Z1_map <- setNames(datC_imp1$Z1, rownames(datC_imp1))
  datA_imp1$Z1 <- factor(Z1_map[idxA], levels = levels(datC_imp1$Z1))

  list(
    datA_fused_2 = datA_imp1,
    out_nnd = out_nnd,
    model_form = formula_nnet,
    FitteddatA_imp1 = fitA,
    FitteddatC_imp1 = fitC
  )
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
#'   \item{beta_BZ}{A numeric vector of estimated model coefficients.}
#'   \item{AA_dummy}{A matrix or data frame of dummy variables constructed for dataset A.}
#'   \item{BB_dummy}{A matrix or data frame of dummy variables constructed for dataset B.}
#'   \item{prob_Z_all}{A numeric matrix of predicted probabilities for each category of variable Z (rows correspond to observations).}
#'   \item{zero_one}{A binary matrix (one-hot encoded) sampled from \code{prob_Z_all}, where each row contains a single 1 indicating the selected category of Z.}
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
#'
#'   datAB <- datAB_to_SM(datA, datB)
#'   datC <- data.frame(
#'     X = datB$X[1:4],
#'     Y1 = datA$Y1[1:4],
#'     Z1 = datB$Z1[1:4]
#'   )
#'
#'   # adding auxiliary information
#'   datABC <- rbind(datAB, datC)
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
#'   step_first_aii <- sma_step1(datA, datB, datC, output_list_AII)
#'   step_second_aii3 <- sma1_step2(step_first_aii)
#'
#'   result_aii3 <- step_second_aii3$datA_fused_2
#'   head(result_aii3)
#' }
sma3_step2 <- function(step1){
  datA_imp1 <- step1$datA_imp1
  datB_imp1 <- step1$datB_imp1

  AB_shared <- intersect(names(datA_imp1), names(datB_imp1))
  p <- ncol(datA_imp1[AB_shared])

  dummy_ll <- function(dat) {
    dummy <- lapply(seq_len(p), function(j) {
      xname <- colnames(dat)[j]
      levs  <- levels(dat[[j]])
      used_levs <- if (j == p) levs else levs[-length(levs)]
      sapply(used_levs, function(lv) as.numeric(dat[[j]] == lv)) |>
        `colnames<-`(paste(xname, used_levs, sep = "_"))
    })
    do.call(cbind, dummy)
  }

  AA_dummy <- dummy_ll(datA_imp1)
  BB_dummy <- dummy_ll(datB_imp1)

  modB <- multinom(datB_imp1$Z1 ~ BB_dummy - 1)
  beta_BZ <- coef(modB)
  ref_Z1 <- levels(datB_imp1$Z1)[1]

  exp_AA_beta <- exp(AA_dummy %*% t(beta_BZ))
  prob_ref <- 1 / (1 + rowSums(exp_AA_beta))
  prob_Z <- exp_AA_beta * prob_ref

  prob_Z_all <- data.frame(prob_ref, prob_Z)
  colnames(prob_Z_all)[1] <- ref_Z1

  zero_one <- apply(prob_Z_all, 1, function(prob) {
    names(prob) <- colnames(prob_Z_all)
    rmultinom(1, 1, prob)
  })

  rownames(zero_one) <- colnames(prob_Z_all)
  Z1_new <- apply(zero_one, 2, function(x) rownames(zero_one)[x == 1])

  datA_fused_2 <- data.frame(datA_imp1[-ncol(datA_imp1)], Z1 = Z1_new)

  list(
    datA_fused_2 = datA_fused_2,
    beta_BZ = beta_BZ,
    AA_dummy = AA_dummy,
    BB_dummy = BB_dummy,
    prob_Z_all = prob_Z_all,
    zero_one = zero_one
  )
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
#' @param step1 output from \code{smc_step1()}
#' @param step2 output from step 2 method
#' @importFrom StatMatch comp.prop
#' @return A list with:
#' \describe{
#'   \item{heli_latent}{A numeric value representing the Hellinger distance between the reference and synthetic distributions.}
#'   \item{ref_distr}{A numeric vector or table representing the reference (original) distribution.}
#'   \item{synt_distr}{A numeric vector or table representing the synthetic (generated) distribution.}
#' }
#' @export
sm_quality <- function(step1, step2) {
  ref_distr  <- table(step1$datB_imp1$Z1)
  synt_distr <- table(step2$datA_fused_2$Z1)

  heli_latent <- comp.prop(
    p1 = synt_distr,
    p2 = ref_distr,
    n1 = nrow(step1$datA_imp1),
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
#'   \item{X}{Category (e.g., "F", "M")}
#'   \item{Y1}{Color category (e.g., "blue")}
#' }
#' @source Simulated data
"datA"

#' Dataset datB
#'
#' A simple dataset with categorical variables.
#'
#' @format A data frame with 15 observations and 2 variables:
#' \describe{
#'   \item{X}{Category (e.g., "F", "M")}
#'   \item{Z1}{Color category (e.g., "red", "green", "blue")}
#' }
#' @source Simulated data
"datB"
