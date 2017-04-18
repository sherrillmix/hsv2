library(dnar)

geneFa<-read.fa('ref/ncbiGenes.txt')
annots<-read.csv('ref/Annotations_table_Hg52.csv',stringsAsFactors=FALSE)
annots<-annots[annots$Sequence.Name!='HSVtest',]
annots$start<-as.numeric(sub(',','',annots$Minimum))
annots$end<-as.numeric(sub(',','',annots$Maximum))
ref<-read.fa('ref/NC_001798 Reference Assembly.fasta')[1,]

genes<-annots[annots$Type=='CDS',]
genes$strand<-genes$Direction=='forward'

genes$row<-stackRegions(genes$start-1000,genes$end+1000)
genes$name<-sub(' CDS$','',genes$Name)
seqs<-substring(ref$seq,genes$start,genes$end)
seqs[!genes$strand]<-revComp(seqs[!genes$strand])
genes$aa<-dna2aa(seqs)
