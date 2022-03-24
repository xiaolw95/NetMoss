#' @title Integrating Networks.
#' @description This is a function to integrate data set from different studies based on network.
#'
#' @param case_dir string.The directory of diseased data set.
#' @param control_dir string.The directory of healthy data set.
#' @param net_case_dir string.The directory of network correlation of diseased data set.
#' @param net_control_dir string.The directory of network correlation of healthy data set.
#'
#' @return integrated diseased and healthy networks.
#' @export
#'
#' @examples
#' case_dir = 'tests/case_dir'
#' control_dir = 'tests/control_dir'
#' net_case_dir = 'tests/net_case_dir'
#' net_control_dir = 'tests/net_control_dir'
#' getNetwork(case_dir = case_dir,
#'      control_dir = control_dir,
#'      net_case_dir = net_case_dir,
#'      net_control_dir = net_control_dir)
#'
getNetwork <-
  function(case_dir,
           control_dir,
           net_case_dir,
           net_control_dir) {
    ######################################1. read data#################################################
    print ("importing datasets...")
    #case network
    setwd(case_dir)
    case_samples <- c()
    dir1 = list.files()
    n1 = length(dir1)
    for (i in 1:n1)
    {
      x.case = read.table(file = dir1[i],
                          header = T,
                          sep = "\t")
      case_samples <- c(case_samples, length(colnames(x.case)) - 1)
    }
    #control data
    setwd(control_dir)
    control_samples <- c()
    dir2 = list.files()
    n2 = length(dir2)
    for (i in 1:n2)
    {
      x.control = read.table(file = dir2[i],
                             header = T,
                             sep = "\t")
      control_samples <-
        c(control_samples, length(colnames(x.control)) - 1)
    }
    #case network
    case_data_list <- list()
    setwd(net_case_dir)
    dir3 = list.files()
    n3 = length(dir3)
    for (m in 1:n3)
    {
      x.case.net = read.table(
        file = dir3[m],
        header = T,
        sep = "\t",
        row.names = 1
      )
      case_data_list[[m]] = x.case.net
    }
    #control network
    control_data_list <- list()
    setwd(net_control_dir)
    dir4 = list.files()
    n4 = length(dir4)
    for (m in 1:n4)
    {
      x.control.net = read.table(
        file = dir4[m],
        header = T,
        sep = "\t",
        row.names = 1
      )
      control_data_list[[m]] = x.control.net
    }

    #union genus
    union_genus1 <- c()
    for (k in 1:length(case_data_list))
    {
      union_genus1 <- union(union_genus1, rownames(case_data_list[[k]]))
    }

    union_genus2 <- c()
    for (k in 1:length(control_data_list))
    {
      union_genus2 <- union(union_genus2, rownames(control_data_list[[k]]))
    }

    union_genus <- union(union_genus1,union_genus2)

    union_matrix <-
      matrix(nrow = length(union_genus),
             ncol = length(union_genus),
             0)

    rownames(union_matrix) <- union_genus
    colnames(union_matrix) <- union_genus
    diag(union_matrix) <- 1


    ######################################2. network construction#################################################
    print ("constructing networks...")
    case_union_data_list <- list()
    for (i in n1) {
      fi <- case_data_list[[i]]
      re_union <- union_matrix
      diff_genus <- setdiff(union_genus, rownames(fi))
      rownames(re_union) <- c(rownames(fi), diff_genus)
      colnames(re_union) <- c(rownames(fi), diff_genus)
      for (j in 1:nrow(fi)) {
        re_union[j, 1:nrow(fi)] <- as.numeric(fi[j, ])
      }
      re_union <- re_union[union_genus, union_genus]
      diag(re_union) <- 1
      case_union_data_list[[i]] <- re_union
    }

    control_union_data_list <- list()
    for (i in n2) {
      fi <- control_data_list[[i]]
      re_union <- union_matrix
      diff_genus <- setdiff(union_genus, rownames(fi))
      rownames(re_union) <- c(rownames(fi), diff_genus)
      colnames(re_union) <- c(rownames(fi), diff_genus)
      for (j in 1:nrow(fi)) {
        re_union[j, 1:nrow(fi)] <- as.numeric(fi[j, ])
      }
      re_union <- re_union[union_genus, union_genus]
      diag(re_union) <- 1
      control_union_data_list[[i]] <- re_union
    }


    ######################################3. combination#################################################
    print ("integrating networks...")
    pool_union <- union_matrix
    for (i in (1:(nrow(pool_union) - 1))) {
      for (j in ((i + 1):nrow(pool_union))) {
        son <- 0
        mom <- 0
        for (k in n1) {
          v <- (1 - case_union_data_list[[k]][i, j] ^ 2) / (case_samples[k] - 1)
          w <- 1 / v
          son <- son + w * case_union_data_list[[k]][i, j]
          mom <- mom + w
        }
        pool_union[i, j] <- son / mom
        pool_union[j, i] <- pool_union[i, j]
      }
    }
    case_union <- pool_union

    pool_union <- union_matrix
    for (i in (1:(nrow(pool_union) - 1))) {
      for (j in ((i + 1):nrow(pool_union))) {
        son <- 0
        mom <- 0
        for (k in n2) {
          v <-
            (1 - control_union_data_list[[k]][i, j] ^ 2) / (control_samples[k] - 1)
          w <- 1 / v
          son <- son + w * control_union_data_list[[k]][i, j]
          mom <- mom + w
        }
        pool_union[i, j] <- son / mom
        pool_union[j, i] <- pool_union[i, j]
      }
    }
    control_union <- pool_union

    integratedNet = list(case_union, control_union)

    return(integratedNet)

  }
