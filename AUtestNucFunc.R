#library(phangorn)
#library(scaleboot) 

AU.test<-function(trees = trees,data = data, bf=bf){

require(phangorn)
require(scaleboot)
trees<-trees
siteLik<-matrix(0,nrow=length(attr(data,"index")),ncol=length(trees))
#bf<-baseFreq(data)
for (i in 1:length(trees)){
	fit<-pml(trees[[i]],data,k=4,bf=bf)
	cat("optimizing paramaters","\n")
	best.fit<-optim.pml(fit,optEdge=FALSE,optGamma=TRUE,optBf=FALSE,optQ=TRUE)
	cat("getting site-wise likelihoods","\n")
	siteLik[,i]<-rep(best.fit$siteLik,attr(data,"weight"))
	}

cat("performing AU test", "\n")
AU<-relltest(siteLik)
AU
}

baseFreq <- function(dat){
    if (class(dat) != "phyDat") 
        stop("data must be of class phyDat")
    levels <- attr(dat,"levels")
    weight <- attr(dat,"weight")
    n <- length(dat)    
    res <- numeric(length(levels))  
    dat <- new2old.phyDat(dat)   
    for(i in 1:n)res <- res+colSums(dat[[i]]*weight)    
    res <- res/sum(res)
    names(res) <- levels    
    res    
}