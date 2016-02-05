#' Genomic filters for the somatic amplicon analysis
#'
#' When paired normal samples are unavailable, we use germline
#' lymphoblast cell lines to identify CNVs and outliers that are
#' common in the germline and unlikely to be somatic.  Filters
#' relevant for the somatic amplicon analysis include germline CNVs,
#' centromeric regions, and assembly gaps.  We additionally defined a
#' 'paired bin filter' containing bins that appear to be linked by
#' improperly spaced reads in the lymphoblast samples.  Aberrantly
#' spaced linked bins could indicate a new junction in the lymphoblast
#' cell line or a sequencing artifact.
#'
#' REFACTORING: Separate filters from the parameter list. Create a
#' single filter object for both the amplicon and deletion analyses.
#'
#' @return a named list
#'
#' @examples
#' filters <- listGenomeFilters()
#' 
#' @export
#' 
#' @param ucsc_build currently ignored as many of the filters are only
#'   developed for UCSC build hg19.
#' 
listGenomeFilters <- function(ucsc_build="hg19"){
  if(ucsc_build != "hg19") stop("Only available for build hg19")
  data(tx_hg19, envir=environment())
  tx <- get("tx_hg19")
  data(binAssemblyGaps_hg19, envir=environment())
  binAssemblyGaps <- get("binAssemblyGaps_hg19")
  data(gaps_hg19, envir=environment())
  gaps <- get("gaps_hg19")
  centromeres <- gaps[gaps$type=="centromere"]
  data(lymphoblast_filters_hg19, envir=environment())
  lymphoblast_filters <- get("lymphoblast_filters_hg19")

  data(normalblood_filters_hg19, envir=environment())
  cnv <- reduce(c(lymphoblast_filters[["amplicon"]],
                  lymphoblast_filters[["deletion"]],
                  normalblood_filters[["amplicon"]],
                  normalblood_filters[["deletion"]]))
  out <- reduce(c(lymphoblast_filters[["outlier"]],
                  normalblood_filters[["outlier"]]))
  list(centromeres=centromeres,
       assembly_gaps=binAssemblyGaps,
       germline_cnv=cnv,
       outliers=out,
       transcripts=tx)
}
