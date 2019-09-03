pairwise_DA_wrapper  <- function(reads, groups, comparisons,
                                 mc.samples = 1000, denom          = "all",
                                 verbose    = TRUE, useMC          = F,
                                 parametric = F   , ignore.posthoc = F,
                                 paired.test = FALSE,   use.splosh = F,
                                 p.threshold = c(0.1),
                                 e.threshold = 1,
                                 xmin = -2.5,
                                 xmax =  2.5,
                                 ymin =  0.01,
                                 ymax =  1){

  library(ALDEx2)
  #generate named output data.frame
  out_df      = data.frame(microbe = row.names(reads))
  test.type   = "wi.eBH"
  ptype       = "BH.adjusted.p.value"
  if(parametric){
    test.type = "we.eBH"
  }
  if(ignore.posthoc){
    test.type = gsub(pattern     = "BH",
                     replacement = "p",
                     x           = test.type)
    ptype     = gsub(pattern     = "BH.adjusted.",
                     replacement = "",
                     x           = ptype)
  }
  if(use.splosh){
    test.type = "cFDR"
    ptype     = 'SPLOSH'
  }
  lst.df = comparisons
  for(comparison in 1:nrow(lst.df)){
    comp_name   <- paste(lst.df[comparison,1], lst.df[comparison,2], sep = " vs ")
    readsA      <- reads[, groups == lst.df[comparison,1]]
    readsB      <- reads[, groups == lst.df[comparison,2]]
    labA        <- rep(paste("a_order_", lst.df[comparison,1]), ncol(readsA))
    labB        <- rep(paste("b_order_", lst.df[comparison,2]), ncol(readsB))
    grouplabels <- c(labA, labB)
    reads.clr   <- aldex.clr(reads = data.frame(readsA, readsB),
                             conds = grouplabels,
                             mc.samples = mc.samples, denom="all", verbose=TRUE, useMC= useMC)
    print("done with clr transform")
    reads.eff   <- aldex.effect(reads.clr, verbose = TRUE, include.sample.summary = FALSE, useMC = useMC)
    print("done with effect size")
    reads.tes   <- aldex.ttest(reads.clr, paired.test = paired.test)
    print("done with ttest")
    if(use.splosh){
      if(parametric) {
      reads.tes$cFDR   <- unlist(splosh(reads.tes$we.ep, dplaces = 10)[1]) # reads.tes$we.ep # unlist(qvalue(reads.tes$we.ep)[3])#
      }
      if(!parametric){
        reads.tes$cFDR <- unlist(splosh(reads.tes$wi.ep, dplaces = 10)[1]) # reads.tes$wi.ep # unlist(qvalue(reads.tes$wi.ep)[3])#
      }
    }
    low.p <-  which(reads.tes[,test.type] <  p.threshold[1])
    high.e <- which(abs(reads.eff$effect) >= e.threshold)
    print(paste("low", ptype, sep = "."))
    print(low.p)
    print("high.e")
    print(high.e)

    plot(x = reads.eff$effect, y = reads.tes[,test.type], log="y",
         pch=19, col=rgb(0,0,0,0.3),
         main=paste("E vs p ", comp_name),
         xlab="effect size",
         ylab=paste("E(", ptype, ")") ,
         cex=0.5,
         ylim = c(ymin, ymax),
         xlim = c(xmin, xmax))
    abline(h =  p.threshold[1], col = "#41b6c4", lty = 2)
    if(length(p.threshold) > 1){
        extra.p.cols = c("#41b6c4",
                       "#1d91c0",
                       "#225ea8",
                       "#253494")
      for(extra.p in 2:min(length(p.threshold), 4)){
        abline(h =  p.threshold[extra.p], col = extra.p.cols[extra.p], lty = 2)
      }
    }
    abline(v =  e.threshold, col = "red",  lty = 2)
    abline(v = -e.threshold, col = "red",  lty = 2)
    if(length(low.p)  > 0){
      points(reads.eff$effect[low.p], reads.tes[,test.type][low.p],
             pch=19,
             col=rgb(0,0,1,0.5),
             cex=0.5)
    }
    if(length(high.e) > 0){
      points(reads.eff$effect[high.e], reads.tes[,test.type][high.e],
             col=rgb(1,0,0,0.5),
             cex=0.8)
    }

    add_df           <- data.frame(rownames(reads.tes),reads.tes[,test.type], reads.eff$effect)
    names(add_df)    <- c("microbe", paste(comp_name, ptype) , paste(comp_name))
    rownames(add_df) <- rownames(reads.tes)


    out_df   <-  merge(out_df, add_df, by.x = "microbe", by.y = "microbe", all.x= T)
  }
  return(out_df)
}
