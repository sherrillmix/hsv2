library(dnar)
source('readAnnotation.R')
if(!file.exists('aaChanges.csv'))source('compareGenomes.R')
aaChanges<-read.csv('aaChanges.csv',stringsAsFactors=FALSE)
#note some genes end early in all three references e.g. UL15. Not sure how to plot
diffs<-lapply(split(aaChanges,1:nrow(aaChanges)),function(xx){
	aas<-seqSplit(xx[,c('msAA','bernAA')],fill='.')
	diffs<-apply(aas,2,function(x)x[1]!=x[2])
	return(which(diffs))
})
names(diffs)<-aaChanges$name
genes[1,]


spacer<-500
pdf('out/genes.pdf',width=11,height=6)
  #ylim magic number
  plot(1,1,type='n',xlim=c(-1000,max(genes$end)+1000),ylim=c(1,10)+c(-.5,.5),ylab='',xlab='Genome position',yaxt='n',bty='n')
  strWidths<-strwidth(genes$name)
  rows<-stackRegions(genes$start-spacer-strWidths-100,genes$end+spacer)
  diffYs<-sapply(aaChanges$name,function(xx)rows[genes$name==xx])
  bases<-sapply(aaChanges$name,function(xx)genes[genes$name==xx,'start'])
  polygon(arrow(ifelse(genes$strand,genes$start,genes$end),ifelse(genes$strand,genes$end,genes$start),rows,arrowLength=2000))
  text(genes$start-100,rows,genes$name,adj=1)
  mapply(function(xx,yy,base)rect(base+xx*3,yy-.5,base+xx*3+2,yy+.5,col='red',border=NA),diffs,diffYs,bases)
dev.off()
