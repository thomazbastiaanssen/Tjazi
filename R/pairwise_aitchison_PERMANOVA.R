pairwise_aitchison_PERMANOVA <- function (clr.samples, groups, relevant.comparisons, adjust.p = T, adj.method = "holm" ) 
{
  library(vegan)
  library(ALDEx2)
  out_df <- data.frame()
  for (number in 1:nrow(relevant.comparisons)) {
    relevant.groups = groups[groups == relevant.comparisons[number, 
                                                            1] | groups == relevant.comparisons[number, 2]]
    comparison_title <- paste(relevant.comparisons[number, 
                                                   1], " VS ", relevant.comparisons[number, 2])
    species.selected <- clr.samples[, groups == relevant.comparisons[number, 
                                                                     1] | groups == relevant.comparisons[number, 2]]
    dist.clr <- dist(t(species.selected))
    ado <- adonis(dist.clr ~ relevant.groups, permutations = 1000)
    out_df[1, number] = ado$aov.tab$`Pr(>F)`[1]
    out_df[2, number] = ado$aov.tab$R2[1]
    out_df[3, number] = ado$aov.tab$F.Model[1]
    colnames(out_df)[number] <- comparison_title
  }
  row.names(out_df) <- c("p.value", "R2", "F.value")
  
  if(adjust.p){
    out_df = rbind(p.adjust(out_df[1,], method = adj.method), 
                   out_df)
    row.names(out_df) <- c("p.adjusted", "p.value", "R2", "F.value")
    
  }
  return(out_df)
}

