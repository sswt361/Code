library(limma)
library(survival)
library(survminer)

expFile="ssgseaOut.txt"    
cliFile="time.txt"          
setwd("C:\\Users\\19806\\Desktop\\193TLS22\\15.survival")  

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tumorData=data[group==0,,drop=F]
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365

sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(cli, data)

res.cut=surv_cutpoint(rt, time = "futime", event = "fustat", variables =c("TLS_score"))
cutoff=as.numeric(res.cut$cutpoint[1])
Group=ifelse(rt[,"TLS_score"]>cutoff, "TLS-High", "TLS-Low")
Group=factor(Group, levels=c("TLS-Low", "TLS-High"))
rt=cbind(as.data.frame(rt), Group)

pdf(file="cutoff.pdf", width=6, height=5, onefile=F)
plot(res.cut, "TLS_score", palette="npg", lwd=3)
dev.off()

outTab=cbind(ID=row.names(rt), rt)
write.table(outTab, file="TLS.group.txt", sep="\t", quote=F, row.names=F)

diff=survdiff(Surv(futime, fustat) ~ Group, data=rt)
pValue=1-pchisq(diff$chisq, df=1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=", sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Group, data = rt)
#print(surv_median(fit))
	
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           #pval=pValue,
		           pval.size=6,
		           legend.title="Group",
		           legend.labs=c("TLS-Low","TLS-High"),
		           font.legend=13,
		           xlab="Time(years)",
		           ylab="Overall survival",
		           break.time.by=2,
		           palette=c("#95CC5EFF", "#FD7446FF"),
		           risk.table=F,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)

res_cox<-coxph(Surv(futime, fustat) ~ Group, data=rt)
surPlot$plot = surPlot$plot + ggplot2::annotate("text", x=2, y=0.12, cex=4.5,
                        label = paste0("log-rank: ", pValue)) +
			ggplot2::annotate("text", x=3.5, y=0.05, cex=4.5,
                        label = paste0("HR: ",round(summary(res_cox)$conf.int[1],3), " (95%CI: ", round(summary(res_cox)$conf.int[3],3), "-", round(summary(res_cox)$conf.int[4],3),")" ))
            #ggplot2::annotate("text", x=6, y=0.8, cex=5,
            #    label = paste0("cutoff: ", round(cutoff,3)))

pdf(file="survival.pdf", width=5, height=4.25, onefile=FALSE)
print(surPlot)
dev.off()
