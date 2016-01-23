#' 1kb bins with metadata for GC and mappability
#'
#' \code{GRanges} objects for 1kb bins of UCSC genome builds hg19 and
#' hg18 with percent GC and average mappability stored in the
#' \code{GRanges} metadata.
#'
#' @docType data
#' @keywords datasets
#' @name bins
#' @usage data(bins1kb_hg19)
#' @usage data(bins1kb_hg18)
#' @aliases bins1kb_hg18 bins1kb_hg19
#'
#' @format a \code{GRanges} object containing 1kb intervals with
#'   metadata columns 'gc' and 'map'.  See details
#'
#' @details
#'
#' The metadata column 'gc' contains the percent GC for the 1kb
#'   intervals.  The percent GC was multiplied by 1000 and rounded to
#'   the nearest integer.
#'
#' The metadata column 'map' contains the average mappability for the
#'   1kb bins.  The average mappability was multiplied by 1000 and
#'   rounded to the nearest integer.
#'
#' @seealso See
#'   \url{ftp://hgdownload1.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig}
#'   for 100mer mappability tracks.
#'
NULL

#' Genome assembly gaps 
#'
#' A \code{GRanges} object representing assembly gaps.  Assembly gaps
#' were defined as 1kb bins with very low GC content (less than 10
#' percent), including 'N' sequences.  The function \code{reduce} was
#' applied to the 1kb bins with low GC, merging adjacent 1kb bins with
#' low GC.
#'
#' @docType data
#' @name binAssemblyGaps
#' @usage data(binAssemblyGaps_hg19)
#' @aliases binAssemblyGaps_hg19
#'
#' @format a \code{GRanges} object with metadata element 'gc'
#'   providing the corresponding GC content for the interval.
#'
NULL

#' Low mappability regions
#'
#' A \code{GRanges} object containing intervals of mappability less
#' than 0.75.
#'
#' @docType data
#' @name lowMappabilityBins
#' @usage data(lowMappabilityBins_hg19)
#' @aliases lowMappabilityBins_hg19
#'
#' @format a \code{GRanges} object
NULL


#' Lymphoblast coverage-based filters
#'
#' Genomic intervals for outliers, deletions, and amplifications
#' identified in lymphoblast cell lines.
#'
#' @details
#'
#' Coverage estimates were preprocessed using \code{sv_preprocess} in
#'   1kb bins.  The autosomal median absolute deviation (MAD) was used
#'   as a robust measure of coverage variation. Bins for which two or
#'   more of the 10 lymphoblast cell lines had a preprocessed coverage
#'   estimate 5 or more MADs from zero were categorized as outliers.
#'   Segmentation of the preprocessed coverage estimates was performed
#'   with circular binary segmentation with default parameters from
#'   \code{SegmentParam}.  Any segment with mean less than -1 or
#'   greater than 1 was considered a deletion or amplification,
#'   respectively.  The genomic intervals for outliers, deletions, and
#'   amplifications were combined in a single object
#'
#' @docType data
#' @name lymphoblast_filters_hg19
#' @usage data(lymphoblast_filters_hg19)
#' @aliases lymphoblast_filters_hg19
#' @format a \code{GRanges} object
NULL
