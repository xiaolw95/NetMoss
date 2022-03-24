#' @title Identifying Differential Bacteria Based on a Network Workflow.
#' @description This is a function to identify differential bacteria between case and control data set.
#'
#' @param case_dir string.The directory of diseased data set.
#' @param control_dir string.The directory of healthy data set.
#' @param net_case_dir string.The directory of network correlation of diseased data set.
#' @param net_control_dir string.The directory of network correlation of healthy data set.
#' @param deepSplit For method "hybrid", can be either logical or integer in the range 0 to 4. For method "tree", must be logical. In both cases, provides a rough control over sensitivity to cluster splitting. The higher the value (or if TRUE), the more and smaller clusters will be produced. For the "hybrid" method, a finer control can be achieved via maxCoreScatter and minGap below.
#' @param minModuleSize Minimum module size.
#'
#' @return NetMoss score and module division results of each bacterium.
#' @export
#'
#' @examples
#' case_dir = '/tests/case_dir'
#' control_dir = '/tests/control_dir'
#' net_case_dir = '/tests/net_case_dir'
#' net_control_dir = '/tests/net_control_dir'
#' NetMoss(case_dir = case_dir,
#'      control_dir = control_dir,
#'      net_case_dir = net_case_dir,
#'      net_control_dir = net_control_dir)
#'
NetMoss <-
  function(case_dir,
           control_dir,
           net_case_dir,
           net_control_dir,
           deepSplit = 4,
           minModuleSize = 20) {

    my.wd = getwd()

    netAll = getNetwork(
      case_dir = case_dir,
      control_dir = control_dir,
      net_case_dir = net_case_dir,
      net_control_dir = net_control_dir
    )

    case_union = netAll[[1]]
    control_union = netAll[[2]]

    modAll = divModule(case_union = case_union, control_union = control_union)
    dissTOM_case = modAll[[1]]
    dissTOM_control = modAll[[2]]
    dynamicMods_case = modAll[[3]]
    dynamicMods_control = modAll[[4]]

    nodes_result = NetzGO(
      control_mat = control_union,
      case_mat = case_union,
      control_dist = dissTOM_control,
      case_dist = dissTOM_case,
      control_mod = as.numeric(dynamicMods_control),
      case_mod = as.numeric(dynamicMods_case),
      scaled = T
    )

    setwd(my.wd)
    return(nodes_result)

  }
