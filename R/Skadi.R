skadi <- function(x, y,
                  max.distance = 1,
                  method = "spearman",
                  grubbs.threshold = 0.05,
                  diagnostic.plot = T,
                  euclid.outlier.check = F,
                  give.uncorrected.p.value = F,
                  xlab = "x",
                  ylab = "y"){


  library(bootstrap)
  library(outliers)

  if (max.distance > 1) {
    warning("Detection of more than one outlier is still in development! Be extra critical.")
  }

  cor.data = data.frame(x = x, y = y)
  correlation.stats = c()
  outlier.suspects  = list()
  grubbres.full = data.frame(index   = (rep(NA, max.distance)),
                             p.value = (rep(NA, max.distance)))
  counter = 0

  #####pseudojackknife
  for(filtersize in 1:max.distance){
    filtermatrix = combn(1:nrow(cor.data),filtersize)
    #print(filtermatrix)


    for(filtercol in 1:ncol(filtermatrix)){
      cor.stat =      cor.test(x      = cor.data[-filtermatrix[,filtercol], 1],
                          y      = cor.data[-filtermatrix[,filtercol], 2],
                          method = method)$p.value
      # print(cor.stat)
      #####Add the correlation strength per iteration

      correlation.stats = c(correlation.stats, cor.stat)
      #print(correlation.stats)
      outlier.suspects[[length(correlation.stats)]] = filtermatrix[,filtercol]
      counter = counter +1
      #print(counter)
    }
    grub.res.part = grubbs.test(correlation.stats[c(((counter-ncol(filtermatrix)):counter))[-1]])
    grubbres.full[filtersize,2] = grub.res.part$p.value

    crit.word = unlist(strsplit(grub.res.part$alternative, " "))[1]
    if(crit.word == "highest"){
      grubbres.full[filtersize,1] = which.max(correlation.stats[c(((counter-ncol(filtermatrix)):counter))[-1]])
    } else if(crit.word == "lowest"){
      grubbres.full[filtersize,1] = which.min(correlation.stats[c(((counter-ncol(filtermatrix)):counter))[-1]])
    }

  }

  #####Decide on the best answer
  if(length(which(grubbres.full$p.value < grubbs.threshold)) > 0){
    if(euclid.outlier.check){

    prop.line   = lm(y ~ x, data = cor.data[-outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]],])

    line.b = unlist(prop.line$coefficients[1])
    line.c = unlist(prop.line$coefficients[1])/-(unlist(prop.line$coefficients[2]))
    dist.points = c()
    for(orig.point in 1:nrow(cor.data)){
      dist.points[orig.point] = distance_from_2d_line(a = c(cor.data[orig.point,1], cor.data[orig.point,2]),
                                                      b = c(0,line.b),
                                                      c = c(line.c,0))
    }

    if(grubbs.test(dist.points)$p.value > grubbs.threshold | !outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]] %in%
      c(which.min(dist.points), which.max(dist.points))){
      print("Found a suspicious value, but grubbs to line was not significant")
      return.this = cor.test(x = cor.data$x,
                             y = cor.data$y,
                             method = method)
      if(max.distance == 1){
      return.this$outlier = NA
      }
      if(!give.uncorrected.p.value){
        return.this$p.value = min((return.this$p.value * (nrow(cor.data)/(nrow(cor.data)-1)))/(1-grubbs.threshold), 1)
      }
      return(return.this)
    }
    }

    #print("hue")
    #print(outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]])
    if(diagnostic.plot){
    plot(x = cor.data[-outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]],]$x,
         y = cor.data[-outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]],]$y,
         xlim = c(min(cor.data$x, na.rm = T), max(cor.data$x, na.rm = T)),
         ylim = c(min(cor.data$y, na.rm = T), max(cor.data$y, na.rm = T)),
         main = "outliers in red",
         xlab = xlab,
         ylab = ylab)
      points(x = cor.data[outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]],]$x,
             y = cor.data[outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]],]$y,
             col = 'red')
      abline(lm(y ~ x, data = cor.data[-outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]],]))


    }
    print("Outlier found")
    return.this = cor.test(x = cor.data[-outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]],]$x,
                   y = cor.data[-outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]],]$y,
                   method = method)
    if(max.distance == 1){
    return.this$outlier = outlier.suspects[[grubbres.full[which.min(grubbres.full$p.value),1]]]
    }
    if(!give.uncorrected.p.value){
      return.this$p.value = min((return.this$p.value * (nrow(cor.data)/(nrow(cor.data)-1)))/(1-grubbs.threshold), 1)
    }
    return(return.this)


  } else {
    print("Did not find outlier")
    return.this = cor.test(x = cor.data$x,
                           y = cor.data$y,
                           method = method)
    if(max.distance == 1){
    return.this$outlier = NA
    }
    if(!give.uncorrected.p.value){
      return.this$p.value = min((return.this$p.value * (nrow(cor.data)/(nrow(cor.data)-1)))/(1-grubbs.threshold), 1)
    }
    return(return.this)


  }

}



######Testing skadi
#resk    <- c()
#resp    <- c()
#respear <- c()
#for(num in 1:1000){
#  print(num)
#  x = rnorm(10)
#  y = rnorm(10)
#  ressk = skadi(x,
#                y,
#                max.distance = 1,
#                grubbs.threshold = 0.05,
#                method = "pearson",
#                euclid.outlier.check = T,
#                give.uncorrected.p.value = F)
#  ressp = cor.test(x, y, method = "spearman")
#  respe = cor.test(x, y, method = "pearson")
#  resk[num] = ressk$p.value
#  resp[num] = ressp$p.value
#  respear[num] = respe$p.value
#}
#
#######

#hist(resk,    breaks = 20, col = rgb(1,0,0,0.5), main = "red = skadi, blue = spearman, green = pearson")
#hist(resp,    breaks = 20, col = rgb(0,0,1,0.5), add = T)
#hist(respear, breaks = 20, col = rgb(0,1,0,0.5), add = T)

#sum(resk    < 0.05)
#sum(resp    < 0.05)
#sum(respear < 0.05)
