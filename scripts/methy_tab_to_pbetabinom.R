
library(VGAM)
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

wd="./"
minCT=1
logp=FALSE
for (e in commandArgs(TRUE)) {
        ta = strsplit(e,"=",fixed=TRUE)
        if(! is.na(ta[[1]][2])) {
                if(ta[[1]][1] == "wd"){
                        wd<-ta[[1]][2]
                }
                if(ta[[1]][1] == "file1"){
                        file1<-ta[[1]][2]
                }
                if(ta[[1]][1] == "file2"){
                        file2<-ta[[1]][2]
                }
                if(ta[[1]][1] == "output"){
                        output<-ta[[1]][2]
                }
                if(ta[[1]][1] == "minCT"){
                        minCT<-as.numeric(ta[[1]][2])
                }
		if(ta[[1]][1] == "logp"){
                        logp<-as.logical(ta[[1]][2])
                }
         }
}


setwd(wd)


d1<-read.table(file1,sep="\t",header=F)
rownames(d1)<-d1[,1]
d2<-read.table(file2,sep="\t",header=F)
rownames(d2)<-d2[,1]
common.names<-intersect(rownames(d1),rownames(d2))


d.common<-cbind(d1[common.names,c(3,4)],d2[common.names,c(3,4)])

f<-cbind(d.common[d.common[,3]>=minCT,2],d.common[d.common[,3]>=minCT,4])
f[,1]=f[,1]/f[,2]

mean1<-mean(f[,1])
var1<-var(f[,1])
shapes1<-estBetaParams(mean1,var1)

p.s1<-pzoibetabinom.ab(f[,1], f[,2], shape1 = shapes1$alpha, shape2 = shapes1$beta, lower.tail=T, log.p=logp)
p.s1=-p.s1
z.s1<-(p.s1-mean(p.s1))/sd(p.s1)
filter.names<-common.names[d.common[,3]>=minCT]

write.table(cbind(filter.names,p.s1,z.s1),output,row.names=F, col.names=F,quote=F,sep="\t")

