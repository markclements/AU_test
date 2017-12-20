#library(phangorn)
#library(scaleboot) 

AU.test<-function(trees = trees,data = data,file = file,Q=Q,bf=bf){

require(phangorn)
require(scaleboot)
trees<-trees
#trees<-read.nexus(file)
#data<-phyDat((read.nexus.data(file)),type="DNA",return.index=TRUE)
data<-data
fit<-pml(trees[[1]],data,k=4),#model="",k=4)
fit<-update(fit,Q=Q,bf=bf)
fiti<-as.vector(1:length(trees))
for (i in 1:length(trees)){
	fiti[i]<-update(fit,tree=trees[[i]])$logLik
	fiti
	}
best<-which.max(fiti)
cat("best tree =",best,"\n")

fit<-update(fit,tree=trees[[best]])

cat("optimizing paramaters","\n")
best.fit<-optim.pml(fit,optEdge=FALSE,optGamma=TRUE,optBf=FALSE)

siteLik<-matrix(0,nrow=length(attr(data,"index")),ncol=length(trees))
cat("getting site-wise likelihoods","\n")

for (i in 1:length(trees)){
	x <- update(best.fit,tree=trees[[i]])
	siteLik[,i]<-rep(x$siteLik,attr(data,"weight"))
	siteLik}
cat("performing AU test", "\n")
AU<-relltest(siteLik)
AU
}