#' @useDynLib lsReg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods is
NULL

#' Allocate memory for large-scale regression
#'
#' Prepares and caches data structures from a fitted base GLM for use in
#' repeated large-scale hypothesis tests via \code{\link{runtest}}.
#'
#' @param basemdl Base model of the form \code{y ~ xl}, fitted with \code{glm}.
#'   Must be of family \code{gaussian}, \code{binomial}, or \code{poisson}.
#' @param colstoadd Number of columns in \code{xr}. The full model tested will
#'   be \code{y ~ xl + xr}.
#' @param testtype Character string specifying the test type. One of
#'   \code{"lrt"} (default), \code{"score"}, \code{"robustscore"},
#'   \code{"wald"}, or \code{"robustwald"}.
#' @return An object of class \code{lsregmem} containing pre-allocated matrices
#'   and cached quantities from the base model, for use with \code{\link{runtest}}.
#' @export
#'
#' @examples
#' linmdlfile <- system.file("extdata", "linmodel.rds", package = "lsReg")
#' linmodel <- readRDS(linmdlfile)
#' link <- lsReg(linmodel, 1)
lsReg <- function(basemdl, colstoadd, testtype) {
  if (missing(basemdl) == TRUE)
    stop("No base model specified")
  if (is(basemdl, "glm") == FALSE)
    stop("Base model is not of class glm")
  mdltype <- match(basemdl$family$family, c("gaussian", "binomial", "poisson"))
  if (is.na(mdltype) == TRUE)
    stop("Family of model must be gaussian, binomial or poisson")
  
  if (missing(colstoadd) == TRUE)
    stop("colstoadd not specified")
  if (is.numeric(colstoadd) == FALSE)
    stop("colstoadd must be an integer value")
  if (length(colstoadd) != 1)
    stop("colstoadd must be a single integer value")
  if (floor(colstoadd) != colstoadd) {
    stop("colstoadd must be an integer value")
  }
  colstoadd <- as.integer(colstoadd)
  if (colstoadd < 1)
    stop("colstoadd must be a postive integer")
  if (missing(testtype) == TRUE)
    testtype <- "lrt"
  if (is.character(testtype) == FALSE)
    stop("testtype must be a character value")
  if (length(testtype) != 1)
    stop("testtype must be a character vector of length 1")
  tt <- match(testtype, c("lrt", "score", "robustscore", "wald", "robustwald"))
  if (is.na(tt) == TRUE)
    stop("Unknown testtype")
  
  lsregmem <- allocate.lsregmem(mdl = basemdl,
                                mdltype = mdltype,
                                q = colstoadd,
                                testtype = testtype)
  
  return (lsregmem)
}