pairwise_aitchison_PERMANOVA <- function(clr.samples, groups, relevant.comparisons){
  library(vegan)
  library(ALDEx2)
  out_df <- data.frame()
  for(number in 1:nrow(relevant.comparisons)){
    relevant.groups = groups[groups == relevant.comparisons[number,1] |groups == relevant.comparisons[number,2]]
    comparison_title <- paste(relevant.comparisons[number,1], " VS ", relevant.comparisons[number,2])
    species.selected <- clr.samples[,groups == relevant.comparisons[number,1] |groups == relevant.comparisons[number,2]]
    dist.clr <- dist(t(species.selected))
    ado <- adonis(dist.clr~relevant.groups, permutations = 999)
    out_df[1,number] = ado$aov.tab$`Pr(>F)`[1]
    colnames(out_df)[number] <- comparison_title
  }
  return(out_df)
}
