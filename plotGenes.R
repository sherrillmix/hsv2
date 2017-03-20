library(dnar)
source('readAnnotation.R')

plot(1,1,type='n',xlim=c(-1000,max(genes$end)+1000),ylim=range(genes$row)+c(-.5,.5))
spacer<-500
strWidths<-strwidth(genes$name)
rows<-stackRegions(genes$start-spacer-strWidths-100,genes$end+spacer)
plot(1,1,type='n',xlim=c(-1000,max(genes$end)+1000),ylim=range(rows)+c(-.5,.5),ylab='',xlab='Genome position',yaxt='n',bty='n')
polygon(arrow(ifelse(genes$strand,genes$start,genes$end),ifelse(genes$strand,genes$end,genes$start),rows,arrowLength=2000))
text(genes$start-100,rows,genes$name,adj=1)
