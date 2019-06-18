gg_volcano_wrapper <- function(DA_df,
                               p.vals           = c(0.01),
                               e.vals           = c(-1, 1),
                               pal.name         = "YlGnBu",
                               xlab             = "Effect Size",
                               ylab             = "log(p-value)",
                               title            = NULL,
                               e.name.threshold = 1,
                               p.name.threshold = 0.05,
                               repel            = T,
                               labelsize        = 2,
                               ylim             = c(0.01, 1),
                               xlim             = c(-3.5, 3.5), 
                               pointsize        = 4
                               ){
  library(RColorBrewer)
  library(ggrepel)
  plot_df = DA_df
  plot_df[,1][with(plot_df, !(abs(plot_df[3]) > e.name.threshold & plot_df[2] < p.name.threshold))] = ""

  p.cols    = brewer.pal(n    = 9,
                         name = pal.name)
  cols.used = p.cols[c(length(p.cols):(length(p.cols)-length(p.vals)))-1]

  volcano = ggplot(data = plot_df,
         aes(x     = plot_df[,3],
             y     = plot_df[,2],
             fill  = cut((plot_df)[,2], c(0, p.vals, 1)),
             label = plot_df[,1])) +
    geom_point(shape = 21, size = pointsize) +
    geom_hline(yintercept = p.vals,
               col = cols.used[length(p.vals):1],
               linetype = "dashed") +
    geom_vline(xintercept = e.vals,
               col        = "red",
               linetype   = "dashed") +
    scale_y_log10(limits = ylim)+
    scale_x_continuous(limits = xlim)+
    scale_fill_manual(values = cols.used,
                      name   = "Legend") +
    theme_bw() +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(aspect.ratio = 1)
  if(repel){
    volcano = volcano + geom_text_repel(size = labelsize)
  } else {
    volcano = volcano + geom_text(size = labelsize)
  }
  volcano
}

