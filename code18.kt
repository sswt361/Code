library(limma)
library(ggpubr)

groupFile="TLS.group.txt"    
cliFile="clinical.txt"       
setwd("C:\\Users\\19806\\Desktop\\193TLS22\\18.clinicalDiff")    

data=read.table(groupFile, header=T, sep="\t", check.names=F, row.names=1)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

sameSample=intersect(row.names(cli), row.names(data))
rt=cbind(cli[sameSample,,drop=F], data[sameSample,"TLS_score",drop=F])

for(clinical in colnames(rt)[1:(ncol(rt)-1)]){
	data=rt[,c("TLS_score", clinical)]
	colnames(data)=c("TLS_score", "clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	group=levels(factor(data$clinical))
	data$clinical=factor(data$clinical, levels=group)
	comp=combn(group,2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
	boxplot=ggboxplot(data, x="clinical", y="TLS_score", fill="clinical",
		          xlab="",
		          ylab="TLS score",
		          palette="aaas",
		          legend.title=clinical)+ 
	    stat_compare_means(comparisons = my_comparisons)
	
	pdf(file=paste0("TLS_", clinical, ".pdf"), width=5, height=4.3)
	print(boxplot)
	dev.off()
}