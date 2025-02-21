expFile="geoMatrix.txt"     
gmtFile="TLS.gmt"         
setwd("C:\\Users\\19806\\Desktop\\193TLS\\16.GEOssGSEA")    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(mat)

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
	rt[rt<0]=0
	rt=log2(rt+1)}
mat=normalizeBetweenArrays(rt)

geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
	
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
	return((x-min(x))/(max(x)-min(x)))}

ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)