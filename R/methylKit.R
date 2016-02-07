#' @import methods GenomicRanges qvalue emdbook
#' 
#' @importMethodsFrom GenomeInfoDb seqnames seqlengths sortSeqlevels
#'
#' @importClassesFrom IRanges RangedData RangedDataList RangedSelection
#' @importMethodsFrom IRanges nearest as.data.frame
#' @importFrom IRanges IRanges RangedDataList RangedData RangesList findOverlaps RangedSelection ranges countOverlaps
#'
#' @importClassesFrom S4Vectors DataTable Annotated Vector List DataTableORNULL characterORNULL SimpleList DataFrame
#' @importFrom S4Vectors DataFrame Rle levels elementMetadata "elementMetadata<-"
#' 
#' @importFrom Rsamtools TabixFile scanTabix yieldSize "yieldSize<-"
#' 
#' @importFrom data.table data.table setnames setcolorder as.data.table tables setkey setkeyv key "key<-" haskey CJ SJ copy rbindlist
#' 
#' @importFrom mclust  densityMclust Mclust
#' 
#' @importFrom rtracklayer export.bed
#' 
## http://stackoverflow.com/questions/8637993/better-explanation-of-when-to-use-imports-depends
## fastseg depends on GenomicRanges, we need to depend on that too in Description
#' @importFrom fastseg fastseg
#' 
#' @importFrom parallel mclapply
#' 
#' @importFrom Rcpp sourceCpp
#' @useDynLib methylKit 
NULL