library(limma)
library(reshape2)
library(ggpubr)

expFile="TLSexp.txt"     
setwd("C:\\Users\\19806\\Desktop\\193TLS22\\12.geneDiff")       

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       
treatNum=length(group[group==0])   
sampleType=c(rep(1,conNum), rep(2,treatNum))

exp=log2(data+1)
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
sigGene=c()
for(i in colnames(exp)[1:(ncol(exp)-1)]){
	if(sd(exp[,i])<0.001){next}
	wilcoxTest=wilcox.test(exp[,i] ~ exp[,"Type"])
	pvalue=wilcoxTest$p.value
	if(wilcoxTest$p.value<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene2=c(sigGene, "Type")
exp=exp[,sigGene2]
diffGeneExp=t(exp[,sigGene])
diffOut=cbind(id=row.names(diffGeneExp), diffGeneExp)
write.table(diffOut, file="diffGeneExp.txt", sep="\t", row.names=F, quote=F)

data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", fill = "Type", 
	     ylab="Gene expression",
	     xlab="",
	     legend.title="Type",
	     palette = c("#95CC5EFF", "#FD7446FF"),
	     width=0.7)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

pdf(file="boxplot.pdf", width=12, height=5.5)
print(p1)
dev.off()
