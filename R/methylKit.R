#' @import methods GenomicRanges qvalue emdbook
#' 
#' @importMethodsFrom GenomeInfoDb seqnames seqlengths sortSeqlevels seqlevels
#'
#' @importMethodsFrom IRanges nearest as.data.frame
#' @importFrom IRanges IRanges findOverlaps RangedSelection ranges countOverlaps
#'
#' @importFrom Rsamtools TabixFile scanTabix yieldSize "yieldSize<-"
#' 
#' @importFrom data.table data.table fread setnames setcolorder as.data.table tables 
#'             setkey setkeyv key "key<-" haskey CJ SJ copy rbindlist setorder
#' 
#' @importFrom mclust  densityMclust Mclust
#' 
#' @importFrom rtracklayer export.bed
#' 
#' @importFrom gtools mixedsort chr
#' 
#' @importFrom R.utils isGzipped
#' 
#' @importFrom limma squeezeVar
#' @importFrom grDevices colorRamp colorRampPalette densCols rainbow rgb topo.colors 
#' @importFrom graphics abline barplot boxplot hist legend lines mtext pairs par pie plot 
#'            points rect smoothScatter strwidth text 
#' @importFrom stats IQR as.dendrogram as.formula binomial 
#'            coefficients cor cor.test dendrapply dist 
#'            dnorm fitted formula glm.fit hclust 
#'            is.leaf kruskal.test lm median model.matrix 
#'            optimize p.adjust pchisq pf pnorm prcomp 
#'            quantile rbeta rnbinom wilcox.test dhyper phyper uniroot
#' @importFrom utils head read.table write.table 
#' 
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
