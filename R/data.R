#' C-Reactive Protein (CRP) Clinical Sample Meausurements
#'
#' @format ## `test_data`
#' A \code{data.table} with 750 rows and 5 columns:
#' \describe{
#'   \item{comparison}{The IVD-MD comparison identifiers.}
#'   \item{SampleID}{The clinical sample identifiers.}
#'   \item{ReplicateID}{The replicate measurement identifiers.}
#'   \item{MP_A}{Measurement results from the first IVD-MD in comparison.}
#'   \item{MP_B}{Measurement results from the second IVD-MD in comparison.}
#' }
#' @source Noklus
"crp_cs_data"

#' C-Reactive Protein (CRP) EQA Material Sample Meausurements
#'
#' @format ## `test_data`
#' A data frame with 120 rows and 5 columns:
#' \describe{
#'   \item{comparison}{The IVD-MD comparison identifiers.}
#'   \item{SampleID}{The EQA material sample identifiers.}
#'   \item{ReplicateID}{The replicate measurement identifiers.}
#'   \item{MP_A}{Measurement results from the first IVD-MD in comparison.}
#'   \item{MP_B}{Measurement results from the second IVD-MD in comparison.}
#' }
#' @source Noklus
"crp_eqam_data"
