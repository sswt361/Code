library(limma)            
expFile="symbol.txt"      
geneFile="gene.txt"      
setwd("C:\\Users\\19806\\Desktop\\193TLS22\\07.TLSexp")    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

outTab=rbind(ID=colnames(geneExp), geneExp)
write.table(outTab, file="TLSexp.txt", sep="\t", quote=F, col.names=F)