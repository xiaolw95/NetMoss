#' @title Plot the combined ROC
#' @description This is a function to plot ROC of combined study.
#'
#' @param case_dir string.The directory of diseased data set.
#' @param control_dir string.The directory of healthy data set.
#' @param marker the combined markers identified by NetMoss.
#' @param metadata the metadata of all input studies.
#' @param plot.roc logical.If TURE then the combined ROC will be plotted.
#' @param train.num numerical. training times for the module. 20 was set by default.
#'
#' @return a TPR and FPR table of the combined study.
#' @export
#'
#' @examples
#' myROC = netROC(case_dir = case_dir,
#'      control_dir = control_dir,
#'      marker = marker,
#'      metadata = metadata,
#'      plot.roc = T,
#'      train.num = 20)
#'
netROC <-
  function(case_dir,
           control_dir,
           marker,
           metadata,
           plot.roc = T,
           train.num = 20)
  {
    my.wd2 = getwd()

    print ("training datasets...")

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

    ##loop
    roc.all = list()
    for (kk in 1:n3)
    {
      eval(parse(text = paste(
        "h = control_data_list[[", kk, "]]", sep = ""
      )))
      eval(parse(text = paste("d = case_data_list[[", kk, "]]", sep = "")))

      ##trim samples
      hh = h
      hh$genus = rownames(hh)
      hh = melt(hh, id.vars = "genus")
      dd = d
      dd$genus = rownames(dd)
      dd = melt(dd, id.vars = "genus")
      all1 = rbind(hh, dd)
      all2 = dcast(all1, genus ~ variable, mean, fill = 0)
      sample.all = all2[, -1]
      rownames(sample.all) = all2$genus
      aa = data.frame(t(sample.all))
      aa$sum = rowSums(aa)
      aa = aa[which(aa$sum != 0), ]
      sample.all2 = data.frame(aa[, -ncol(aa)])
      m.marker = intersect(as.character(colnames(sample.all2)), as.character(rownames(marker)))
      m.meta = metadata[as.character(rownames(metadata)) %in% as.character(rownames(sample.all2)), ]

      eval(parse(text = paste(
        "auc_crc", kk, "_all = data.frame()", sep = ""
      )))
      for (mm in 1:train.num)
      {
        ##train and test data
        train.data = sample(rownames(sample.all2), nrow(sample.all2) * 0.7, replace = F)
        sample.train = sample.all2[as.character(rownames(sample.all2)) %in% as.character(train.data), ]
        sample.train$group = m.meta[as.character(rownames(sample.train)), "type"]
        test.data = setdiff(rownames(sample.all2), train.data)
        sample.test = sample.all2[as.character(rownames(sample.all2)) %in% as.character(test.data), ]
        sample.test$group = m.meta[as.character(rownames(sample.test)), "type"]

        ##training
        if (length(intersect(m.marker, colnames(sample.train))) == 1)
        {
          sample.train2 = data.frame(sample.train[, as.character(colnames(sample.train)) %in% as.character(m.marker)])
          colnames(sample.train2) = as.character(intersect(m.marker, colnames(sample.train)))
          rownames(sample.train2) = rownames(sample.train)
        }else
        {
          sample.train2 = sample.train[, as.character(colnames(sample.train)) %in% as.character(m.marker)]
        }
        sample.train2$group = sample.train$group
        rf.train = randomForest(
          group ~ .,
          sample.train2,
          ntree = 1000,
          nperm = 100,
          importance = T
        )
        importance_rf = data.frame(importance(rf.train))
        importance_rf = importance_rf[order(importance_rf$MeanDecreaseAccuracy, decreasing = T), ]

        #########10-fold cross validation
        train.5_10 = replicate(5,
                               rfcv(
                                 sample.train2[-ncol(sample.train2)],
                                 sample.train2$group,
                                 cv.fold = 10,
                                 step = 1.5
                               ),
                               simplify = F)
        train.5_10.2 = data.frame(sapply(train.5_10, '[[', 'error.cv'))
        train.5_10.2$names = rownames(train.5_10.2)
        train.5_10.2 = melt(train.5_10.2, id = 'names')
        train.5_10.2$names = as.numeric(as.character(train.5_10.2$names))
        train.5_10.2 = summaryBy(value ~ names, train.5_10.2, FUN = mean)

        ####re-trainning
        marker.num = min(train.5_10.2[which(train.5_10.2$value.mean == min(train.5_10.2$value.mean)), 1])
        marker.re = data.frame(rownames(importance_rf[1:marker.num, ]))
        colnames(marker.re) = "Name"
        if (length(intersect(marker.re$Name, colnames(sample.train))) == 1)
        {
          sample.train3 = data.frame(sample.train[, as.character(colnames(sample.train)) %in% as.character(marker.re$Name)])
          colnames(sample.train3) = as.character(intersect(marker.re$Name, colnames(sample.train)))
          rownames(sample.train3) = rownames(sample.train)
        }else
        {
          sample.train3 = sample.train[, as.character(colnames(sample.train)) %in% as.character(marker.re$Name)]
        }
        sample.train3$group = sample.train$group
        rf.train = randomForest(
          group ~ .,
          sample.train3,
          ntree = 1000,
          nperm = 100,
          importance = T
        )

        ##testing
        test_pred.select = predict(rf.train, sample.test)
        test_freq.select = table(test_pred.select, sample.test$group)

        ##ROC
        group_select = data.frame(test_pred.select)
        group_select$test_orig.select = sample.test[as.character(rownames(group_select)), 'group']
        e11 = group_select
        e11$test_pred.select = as.numeric(e11$test_pred.select)
        e11$test_orig.select = as.numeric(e11$test_orig.select)
        e11[which(e11$test_orig.select != 1), 2] = 0
        e11[which(e11$test_pred.select != 1), 1] = 0
        for (i in 1:nrow(e11))
          ifelse(e11[i, 1] == e11[i, 2], e11[i, 3] <- 1, e11[i, 3] <- 0)
        b = data.frame(predict(rf.train, sample.test, type = 'prob'))
        e11$prob = b[as.character(rownames(e11)), "disease"]

        ##plot
        p = ggplot(e11, aes(d = test_orig.select, m = prob)) +
          geom_roc(n.cuts = 0) +
          style_roc()
        auc <- calc_auc(p)

        #roc
        eval(parse(text = paste("roc_crc", kk, "_", mm, " = e11", sep = "")))
        eval(parse(text = paste(
          "auc_crc", kk, "_all[", mm, ",1] = ", mm, sep = ""
        )))
        eval(parse(text = paste(
          "auc_crc", kk, "_all[", mm, ",2] = auc$AUC", sep = ""
        )))

        #pick
        if (auc$AUC > 0.7)
        {
          eval(parse(
            text = paste("roc.all[[", kk, "]] = roc_crc", kk, "_", mm, sep = "")
          ))
          break
        }
        if (mm == train.num)
        {
          eval(parse(
            text = paste(
              "m.auc = which(auc_crc",
              kk,
              "_all$V2 == max(auc_crc",
              kk,
              "_all$V2))",
              sep = ""
            )
          ))
          eval(parse(
            text = paste("roc.all[[", kk, "]] = roc_crc", kk, "_", m.auc, sep = "")
          ))
        }
      }

    }

    ##
    for (nn in 1:n3)
    {
      eval(parse(text = paste("crc", nn, " = roc.all[[", nn, "]]", sep = "")))
      eval(parse(text = paste("crc", nn, "$type = 'CRC", nn, "'", sep = "")))
    }

    if (n3 > 1 )
    {
      mydata = crc1
      for (nn in 2:n3)
      {
        eval(parse(text = paste(
          "mydata = rbind(mydata,crc", nn, ")", sep = ""
        )))
      }

      ##weighted
      w.case = data.frame(sample = case_samples)
      w.case$w1 = 1 / w.case$sample * sum(w.case$sample)
      w.case$w2 = w.case$w1 / sum(w.case$w1)
      y <- mydata$test_orig.select
      w = rep(w.case[1, 3], nrow(crc1))
      for (nn in 2:n3)
      {
        eval(parse(
          text = paste("w = c(w,rep(w.case[", nn, ",3],nrow(crc", nn, ")))", sep = "")
        ))
      }
    }else
    {
      mydata = crc1

      ##weighted
      w.case = data.frame(sample = case_samples)
      w.case$w1 = 1 / w.case$sample * sum(w.case$sample)
      w.case$w2 = w.case$w1 / sum(w.case$w1)
      y <- mydata$test_orig.select
      w = rep(w.case[1, 3], nrow(crc1))
    }

    y.hat <- mydata$prob
    tp.fp <- WeightedROC(y.hat, y, w)

    ##plot combined ROC
    if (plot.roc)
    p2 =  ggplot() +
      geom_path(aes(FPR, TPR), data = tp.fp) +
      coord_equal() +
      annotate("text",
               x = .75,
               y = .15,
               ## 注释text的位置
               label = paste("AUC =", round(WeightedAUC(tp.fp), 2))) +
      labs(x = "False positive fraction", y = "True positive fraction", title = "Combined ROC") +
      theme_bw()

    setwd(my.wd2)
    ggsave("NetMoss_ROC.png",p2)

    tp.fp2 = tp.fp[,c(3,1,2)]
    return(tp.fp2)

  }
