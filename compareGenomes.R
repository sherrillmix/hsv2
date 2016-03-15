source("~/scripts/R/dna.R")
#library(levenR)
library(dnaplotr)
#depends on directories:
#  work (stores temp files)
#  out (stores output files)
#  ref (stores consensuses, annotations and reference sequences)
#  dat (raw reads, dont need yet)

ref<-read.fa('ref/NC_001798 Reference Assembly.fasta')[1,]
bern<-read.fa('ref/BernStein_Consensus.fasta')
ms<-read.fa('ref/MS_Consensus.fasta')

#align<-levenAlign(c(ref$seq,bern$seq,ms$seq),ref$seq,nThreads=3)
if(!file.exists('work/align.blat'))blatReadsVsRefs(c('ref'=ref$seq,'bern'=bern$seq,'ms'=ms$seq),c('ref'=ref$seq),'work/align.blat',faToTwoBit='~/installs/blat/faToTwoBit','gfServer'='~/installs/blat/gfServer','gfClient'='~/installs/blat/gfClient',options='-maxNtSize=1000000')
align<-readBlat('work/align.blat')
align<-align[align$score>100000,]
refAligns<-blockToAlign(c('bern'=bern$seq,'ms'=ms$seq,'ref'=ref$seq)[align$qName],ref$seq,align$qStarts,align$tStarts,align$blockSizes)
combo<-combineAligns(refAligns$tSeq,refAligns$qSeq)
names(combo)<-align$qName
seqMat<-do.call(rbind,strsplit(combo,''))
rownames(seqMat)<-names(combo)

mismatch<-outer(1:3,1:3,function(x,y){mapply(function(xx,yy)sum(seqMat[xx,]!=seqMat[yy,]),x,y)})
table(seqMat['ms',]==seqMat['ref',],seqMat['ms',]==seqMat['bern',])

geneFa<-read.fa('ref/ncbiGenes.txt')
annots<-read.csv('ref/Annotations_table_Hg52.csv',stringsAsFactors=FALSE)
annots<-annots[annots$Sequence.Name!='HSVtest',]
annots$start<-as.numeric(sub(',','',annots$Minimum))
annots$end<-as.numeric(sub(',','',annots$Maximum))

genes<-annots[annots$Type=='CDS',]
genes$refATG<-NA
refSeq<-paste(seqMat['ref',],collapse='')
genes$gapStart<-noGap2Gap(refSeq,genes[,'start'])
genes$gapStartPlus2<-noGap2Gap(refSeq,genes[,'start']+2)
genes$gapEnd<-noGap2Gap(refSeq,genes[,'end'])
genes$gapEndMinus2<-noGap2Gap(refSeq,genes[,'end']-2)

#find the codons
genes$startCodon<-substring(refSeq,ifelse(genes$Direction=='forward',genes$gapStart,genes$gapEndMinus2),ifelse(genes$Direction=='forward',genes$gapStartPlus2,genes$gapEnd))
genes$stopCodon<-substring(refSeq,ifelse(genes$Direction=='forward',genes$gapEndMinus2,genes$gapStart),ifelse(genes$Direction=='forward',genes$gapEnd,genes$gapStartPlus2))
genes$msAA<-substring(paste(seqMat['ms',],collapse=''),genes$gapStart,genes$gapEnd)
genes$msAA<-sub('X.*$','X',replaceAfterStop(dna2aa(degap(ifelse(genes$Direction=='forward',genes$msAA,revComp(genes$msAA))))))
genes$bernAA<-substring(paste(seqMat['bern',],collapse=''),genes$gapStart,genes$gapEnd)
genes$bernAA<-sub('X.*$','X',replaceAfterStop(dna2aa(degap(ifelse(genes$Direction=='forward',genes$bernAA,revComp(genes$bernAA))))))



#TODO missing one ATG in reverse. is this correct?
table(genes$startCodon,genes$Direction)

genes$bernMsDiff<-apply(genes[,c('gapStart','gapEnd')],1,function(x)sum(
	seqMat['bern',x[1]:x[2]]!=seqMat['ms',x[1]:x[2]] &
	seqMat['bern',x[1]:x[2]] %in% c('A','C','G','T','-') &
	seqMat['ms',x[1]:x[2]] %in% c('A','C','G','T','-')
))

table(genes$bernAA==genes$msAA|!grepl('z',genes$bernAA)|!grepl('z',genes$msAA),genes$bernMsDiff==0)
table(genes$bernAA==genes$msAA,genes$bernMsDiff==0,grepl('z',genes$bernAA)|grepl('z',genes$msAA))

aaChanges<-genes[genes$bernAA!=genes$msAA&!grepl('z',genes$bernAA)&!grepl('z',genes$msAA),]

pdf('out/aaChanges.pdf')
for(ii in 1:nrow(aaChanges)){
	message(ii)
	#thisDescription<-myTable[myTable$Gene==sub(' CDS','',genes[1,'Name']),'Description']
	plotAA(aaChanges[ii,c('msAA','bernAA')],main=paste(aaChanges[ii,'Name'],groups=c('ms','bern'))
	aas<-seqSplit(aaChanges[ii,c('msAA','bernAA')],fill='.')
	diffs<-apply(aas,2,function(x)x[1]!=x[2])
	axis(3,which(diffs),paste(aas[1,which(diffs)],which(diffs),aas[2,which(diffs)],sep=''),mgp=c(3,.35,0),tcl=-.3)
}
dev.off()





geneId<-61
interestingSeqs<-apply(seqMat[,genes[geneId,'gapStart']:genes[geneId,'gapEnd']],1,paste,collapse='')
png('test.png',width=4000,height=2000,res=200)
plotDNA(interestingSeqs,groups=rownames(seqMat))
title(main=genes[geneId,'Name'])
dev.off()

