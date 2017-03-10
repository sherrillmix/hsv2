library(dnar)
source('readAnnotation.R')

plot(1,1,type='n',xlim=c(-1000,max(genes$end)+1000),ylim=range(genes$row))
polygon(arrow(ifelse(genes$strand,genes$start,genes$end),ifelse(genes$strand,genes$end,genes$start),genes$row,arrowLength=2000))
