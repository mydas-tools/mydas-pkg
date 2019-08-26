# data.R - DESC
# mydas/R/data.R
#' mydas datasets
#'
#' Example datasets for the classes defined in FLCore.
#'
#' \itemize{
#'    \item{\code{lh}, \code{\link{FLPar}}}{A dataset with life history parameteers,
#'     extracted from fish base.}
#' \item{\code{prior}, \code{\link{FLPar}}}{A dataset with the priors for the life
#'     history parameters, estimated from \code{prior} method.}
#' \item{om, \code{\link{FLStock}}}{An Operation Model based on turbot conditioned from
#'  \code{lh}}
#' \item{eq, \code{\link{FLBRP}}}{Turbot equilibrium object with reference points}
#'
#' Datasets can be loaded by issuing the \code{data} command, as in
#' \code{data(lh)}.
#'
#' @name datasets
#' @aliases lh prior om eq
#' @seealso \linkS4class{FLStock}, \linkS4class{FLBRP}, \linkS4class{FLPar}
#' @keywords datasets
#' @examples
#' \dontrun{
#' data(om)
#' summary(om)
#' }
NULL
