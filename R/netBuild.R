#' @title Building microbial networks
#' @description This is a function to build microbial networks using abundance matrix.
#'
#' @param case_dir string.The directory of diseased data set.
#' @param control_dir string.The directory of healthy data set.
#' @param net_case_dir string.The directory of network correlation of diseased data set.
#' @param net_control_dir string.The directory of network correlation of healthy data set.
#' @param method string. The method chosen to build microbial networks which can be "sparcc" or "pearson".
#'
#' @return microbial networks of case and control data set.
#' @export
#'
#' @examples
#' case_dir = '/tests/case_dir'
#' control_dir = '/tests/control_dir'
#' net_case_dir = '/tests/net_case_dir'
#' net_control_dir = '/tests/net_control_dir'
#' netBuild(case_dir = case_dir,
#'      control_dir = control_dir,
#'      net_case_dir = net_case_dir,
#'      net_control_dir = net_control_dir,
#'      method = "sparcc")
#'
#'
netBuild <- function(case_dir,
                     control_dir,
                     net_case_dir,
                     net_control_dir,
                     method = "sparcc"
) {
  my.wd1 = getwd()

  print ("building networks...")

  #case dir
  case_data_list <- list()
  case_samples <- c()
  setwd(case_dir)
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
    case_samples <-
      c(case_samples, length(colnames(x.case.net)) - 1)
  }
  #control dir
  control_data_list <- list()
  control_samples <- c()
  setwd(control_dir)
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
    control_samples <-
      c(control_samples, length(colnames(x.control.net)) - 1)
  }

  ###sparcc
  if (method == "sparcc")
  {
    ###case network
    for (kk in 1:n3)
    {
      eval(parse(text = paste(
        "d = case_data_list[[", kk, "]]", sep = ""
      )))
      aa = t(d)
      sparcc.amgut = sparcc(aa)
      adjacency_weight_sparcc = sparcc.amgut$CORR
      colnames(adjacency_weight_sparcc) <- colnames(aa)
      rownames(adjacency_weight_sparcc) <- colnames(aa)
      write.table(
        adjacency_weight_sparcc,
        paste0(net_case_dir, "/d_study_net", kk, ".txt"),
        col.names = NA,
        sep = '\t',
        quote = FALSE
      )

    }
    ###control network
    for (kk in 1:n4)
    {
      eval(parse(text = paste(
        "h = control_data_list[[", kk, "]]", sep = ""
      )))
      aa = t(h)
      sparcc.amgut = sparcc(aa)
      adjacency_weight_sparcc = sparcc.amgut$CORR
      colnames(adjacency_weight_sparcc) <- colnames(aa)
      rownames(adjacency_weight_sparcc) <- colnames(aa)
      write.table(
        adjacency_weight_sparcc,
        paste0(net_control_dir, "/h_study_net", kk, ".txt"),
        col.names = NA,
        sep = '\t',
        quote = FALSE
      )
    }
  }

  ###pearson
  if (method == "pearson")
  {
    for (kk in 1:n3)
    {
      eval(parse(text = paste(
        "d = case_data_list[[", kk, "]]", sep = ""
      )))
      cor = corr.test(t(d))$r
      pv = corr.test(t(d))$p
      cor[pv > 0.05] = 0
      write.table(
        cor,
        paste0(net_case_dir, "/d_study_net", kk, ".txt"),
        col.names = NA,
        sep = '\t',
        quote = FALSE
      )
    }
    ###control network
    for (kk in 1:n3)
    {
      eval(parse(text = paste(
        "h = control_data_list[[", kk, "]]", sep = ""
      )))
      cor = corr.test(t(d))$r
      pv = corr.test(t(d))$p
      cor[pv > 0.05] = 0
      write.table(
        cor,
        paste0(net_control_dir, "/h_study_net", kk, ".txt"),
        col.names = NA,
        sep = '\t',
        quote = FALSE
      )
    }
  }

  setwd(my.wd1)
}
