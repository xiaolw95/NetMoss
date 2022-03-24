#' @title Dividing Modules.
#' @description This is a function to divide integrated networks into different modules.
#'
#' @param case_union a matrix of integrated diseased network.
#' @param control_union a matrix of integrated healthy network.
#' @param deepSplit For method "hybrid", can be either logical or integer in the range 0 to 4. For method "tree", must be logical. In both cases, provides a rough control over sensitivity to cluster splitting. The higher the value (or if TRUE), the more and smaller clusters will be produced. For the "hybrid" method, a finer control can be achieved via maxCoreScatter and minGap below.
#' @param minModuleSize Minimum module size.
#'
#' @return divided modules.
#' @export
#' @importFrom WGCNA TOMsimilarity
#'
#' @examples
#' Data(TestData)
#' divModule(TestData[[1]], TestData[[2]])
#'
divModule <- function(case_union, control_union, deepSplit = 4, minModuleSize = 20) {

  print ("dividing modules...")

  adj_control <- as.matrix(0.5 + 0.5 * control_union)
  adj_case <- as.matrix(0.5 + 0.5 * case_union)

  pow_adjacency <- function(mat, pow) {
    mat_soft <- mat
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        mat_soft[i, j] <- mat[i, j] ^ pow
      }
    }
    return(mat_soft)
  }

  ##health network modularization
  for (i in 7) {
    soft_pow_control <- 4
    adj_soft_control <- pow_adjacency(adj_control, soft_pow_control)
    TOM_control <- TOMsimilarity(adj_soft_control)
    dissTOM_control <- 1 - TOM_control
    # Call the hierarchical clustering function
    geneTree_control = hclust(as.dist(dissTOM_control), method = "average")

    minModuleSize = minModuleSize

    # Module identification using dynamic tree cut:
    dynamicMods_control = cutreeDynamic(
      dendro = geneTree_control,
      distM = dissTOM_control,
      deepSplit = deepSplit,
      pamRespectsDendro = FALSE,
      minClusterSize = minModuleSize
    )

    print(table(dynamicMods_control))
    dynamicColors_control = labels2colors(dynamicMods_control)
  }
  ##disease network modularization
  for (i in 7) {
    soft_pow_case <- 4
    adj_soft_case <- pow_adjacency(adj_case, soft_pow_case)
    TOM_case <- TOMsimilarity(adj_soft_case)
    dissTOM_case <- 1 - TOM_case
    # Call the hierarchical clustering function
    geneTree_case = hclust(as.dist(dissTOM_case), method = "average")

    # We like large modules, so we set the minimum module size relatively high:
    minModuleSize = minModuleSize

    # Module identification using dynamic tree cut:
    dynamicMods_case = cutreeDynamic(
      dendro = geneTree_case,
      distM = dissTOM_case,
      deepSplit = deepSplit,
      pamRespectsDendro = FALSE,
      minClusterSize = minModuleSize
    )

    print(table(dynamicMods_case))
    dynamicColors_case = labels2colors(dynamicMods_case)
  }

  modDivision = list(dissTOM_case,
                     dissTOM_control,
                     dynamicMods_case,
                     dynamicMods_control)

  return(modDivision)

}
