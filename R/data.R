#' Aging acinar cells
#'
#' The variable 'donor_age' (age of donor in years) was used as sequential labels. Taken from publication 
#' "Single-Cell Analysis of Human Pancreas Reveals Transcriptional Signatures of Aging and Somatic Mutation Patterns"
#' GSE81547
#'
#' @format A SingleCellExperiment object with [genes] by [cells].
#' @source \url{https://www.sciencedirect.com/science/article/pii/S009286741731053X}
"acinar_sce"

#' Developing human female germline cells
#'
#' The variable 'time' (age in weeks) was used as sequential labels. Taken from publication 
#' "Single-Cell RNA-Seq Analysis Maps Development of Human Germline Cells and Gonadal-Niche Interactions"
#' GSE86146
#'
#' @format A SingleCellExperiment object with [genes] by [cells].
#' @source \url{https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(17)30078-4}
"germ_sce"

#' Developing beta cells
#'
#' The variable 'age' (developmental stage) was used as sequential labels. Taken from publication 
#' "Deciphering Pancreatic Islet β Cell and α Cell Maturation Pathways and Characteristic Features at the Single-Cell Level"
#' GSE87375
#'
#' @format A SingleCellExperiment object with [genes] by [cells].
#' @source \url{https://www.cell.com/cell-metabolism/fulltext/S1550-4131(17)30208-5}
"beta_sce"

#' Human embryonic stem cells
#'
#' The variable 'esc_day' (embryonic day) was used as sequential labels. Taken from publication 
#' "Single-cell RNA-seq reveal lineage formation and X-chromosome dosage compensation in human preimplantation embryos"
#' E-MTAB-3929
#'
#' @format A SingleCellExperiment object with [genes] by [cells].
#' @source \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3929/}
"hesc_sce"

#' MEFs reprogrammed to neurons
#'
#' The variable 'time_point' (days since induction) was used as sequential labels. Taken from publication 
#' "Dissecting direct reprogramming from fibroblast to neuron using single-cell RNA-seq"
#' GSE67310
#'
#' @format A SingleCellExperiment object with [genes] by [cells].
#' @source \url{https://www.nature.com/articles/nature18323}
"mef_sce"

#' [human?] colon cells
#'
#' The sequential labels were derived from unsupervised clustering, using Seurat. Taken from publication 
#' "Unsupervised Trajectory Analysis of Single-Cell RNA-Seq and Imaging Data Reveals Alternative Tuft Cell Origins in the Gut"
#' GSE102698
#'
#' @format A SingleCellExperiment object with [genes] by [cells].
#' @source \url{https://www.cell.com/cell-systems/fulltext/S2405-4712(17)30449-0}
"colon_sce"