#' @title Calculating NetMoss Score
#' @description This is a function to calculate NetMoss score for each bacterium.
#'
#' @param control_mat a matrix of integrated healthy network.
#' @param case_mat a matrix of integrated diseased network.
#' @param control_dist a distance matrix among healthy models.
#' @param case_dist a distance matrix among diseased models.
#' @param control_mod module division of each microbe in healthy network.
#' @param case_mod module division of each microbe diseased network
#' @param scaled logical.If TURE then the NetMoss score are scaled.
#'
#' @return NetMoss score and module division results of each bacterium.
#' @export
#'
#' @examples
#' data(TestData)
#' NetzGo(control_mat = TestData[[2]], case_mat = TestData[[1]],
#'     control_dist = TestData[[4]], case_dist = TestData[[3]],
#'     control_mod = TestData[[6]], case_mod = TestData[[5]], scaled = T)
#'
NetzGO <-
  function(control_mat,
           case_mat,
           control_dist,
           case_dist,
           control_mod,
           case_mod,
           scaled) {

    if (scaled) {
      print("distance matrix scaled")
      control_dist_scaled <-
        (control_dist - min(control_dist)) / (max(control_dist) - min(control_dist))
      case_dist_scaled <-
        (case_dist - min(case_dist)) / (max(case_dist) - min(case_dist))
      control_dist <- control_dist_scaled
      case_dist <- case_dist_scaled
    }

    dynamicMods_inter <-
      paste("control_mod_", control_mod, "-", "case_mod_", case_mod, sep = "")
    dynamicMods_inter_names <- names(table(dynamicMods_inter))

    taxon_names <- rownames(control_mat)
    nodes <-
      data.frame(taxon_names, control_mod, case_mod, dynamicMods_inter)
    nodes$taxon_names <- as.character(nodes$taxon_names)
    nodes$control.degree <- rowSums(control_mat) - 1
    nodes$case.degree <- rowSums(case_mat) - 1

    dynamicMods_control_names <- as.numeric(names(table(control_mod)))
    dynamicMods_case_names <- as.numeric(names(table(case_mod)))

    nodes$control.mod.degree <- NA
    nodes$case.mod.degree <- NA
    for (i in dynamicMods_control_names) {
      row_num <- which(control_mod == i)
      col_num <- which(control_mod == i)
      nodes$control.mod.degree[row_num] <-
        rowSums(control_mat[row_num, col_num]) - 1 - rowSums(control_mat[row_num, -col_num])
    }

    for (i in dynamicMods_case_names) {
      row_num <- which(case_mod == i)
      col_num <- which(case_mod == i)
      nodes$case.mod.degree[row_num] <-
        rowSums(case_mat[row_num, col_num]) - 1 - rowSums(case_mat[row_num, -col_num])
    }

    nodes$diff.mod.degree <-
      (nodes$case.mod.degree - min(nodes$case.mod.degree)) / (max(nodes$case.mod.degree) -
                                                                min(nodes$case.mod.degree)) - (nodes$control.mod.degree - min(nodes$control.mod.degree)) /
      (max(nodes$control.mod.degree) - min(nodes$control.mod.degree))
    nodes$diff.mod.degree.abs <-
      abs((nodes$case.mod.degree - min(nodes$case.mod.degree)) / (max(nodes$case.mod.degree) -
                                                                    min(nodes$case.mod.degree)) - (nodes$control.mod.degree - min(nodes$control.mod.degree)) /
            (
              max(nodes$control.mod.degree) - min(nodes$control.mod.degree)
            )
      )

    control.nodes.to.mods.dist <-
      as.data.frame(matrix(
        nrow = length(taxon_names),
        ncol = length(dynamicMods_inter_names),
        NA
      ))
    rownames(control.nodes.to.mods.dist) <- taxon_names
    colnames(control.nodes.to.mods.dist) <- dynamicMods_inter_names
    case.nodes.to.mods.dist <- control.nodes.to.mods.dist

    nodes$netzgo <- 0
    for (i in dynamicMods_inter_names) {
      row_num <- which(dynamicMods_inter == i)
      for (j in row_num) {
        for (l in dynamicMods_inter_names) {
          col_num <- which(dynamicMods_inter == l)
          control.nodes.to.mods.dist[j, l] <-
            mean(control_dist[j, col_num])
          case.nodes.to.mods.dist[j, l] <- mean(case_dist[j, col_num])
        }
      }
    }
    diff.nodes.to.modes.dist <-
      case.nodes.to.mods.dist - control.nodes.to.mods.dist

    for (i in dynamicMods_inter_names) {
      target_num <- which(dynamicMods_inter == i)
      target_control_mod <- strsplit(i, "-")[[1]][1]
      target_case_mod <- strsplit(i, "-")[[1]][2]

      case_mod_names <-
        dynamicMods_inter_names[grep(target_control_mod, dynamicMods_inter_names)]
      case_mod_names <- setdiff(case_mod_names, i)
      control_mod_names <-
        dynamicMods_inter_names[grep(target_case_mod, dynamicMods_inter_names)]
      control_mod_names <- setdiff(control_mod_names, i)

      for (j in case_mod_names) {
        nodes$netzgo[target_num] <-
          nodes$netzgo[target_num] + diff.nodes.to.modes.dist[target_num, dynamicMods_inter_names %in% j]
      }
      for (k in control_mod_names) {
        nodes$netzgo[target_num] <-
          nodes$netzgo[target_num] - diff.nodes.to.modes.dist[target_num, dynamicMods_inter_names %in% k]
      }
    }
    nodes$netzgo.scale <- NA
    for (i in names(table(dynamicMods_inter))) {
      target_num <- which(nodes$dynamicMods_inter == i)
      for (j in target_num) {
        nodes$netzgo.scale[j] <-
          (nodes$netzgo[j] - min(nodes$netzgo[target_num])) / (max(nodes$netzgo[target_num]) - min(nodes$netzgo[target_num]))
      }
    }
    if (length(which(is.na(nodes$netzgo.scale))) != 0)
    {
      nodes = nodes[-which(is.na(nodes$netzgo.scale)),]
    }


    nodes = nodes[, c("taxon_names",
                      "control_mod",
                      "case_mod",
                      "netzgo.scale")]
    colnames(nodes) = c("taxon_names",
                        "control_mod",
                        "case_mod",
                        "NetMoss_Score")

    return(nodes)

  }
