#library(phangorn)
#library(scaleboot) 

AU.test<-function(trees, data){ # returns the result of the relltest function of scaleboot
require(ape) # load necessary packages
require(phangorn)
require(scaleboot)

if (class(trees) != "multiPhylo") stop("trees must be of class multiPhylo, see ape documenation")
if (class(data) != "phyDat") stop("data must be in phyDat format, see phangon documentation")

bf<-as.vector(baseFreq(data)) #calculate emperical base freq

fit<-lapply(trees,pml,data=data,k=4,bf=bf) # calcualte likelihood of all trees

best<-(unlist(lapply(fit,FUN=function(x) x$logLik))) # find best likelihood and tree

cat("best tree =",max(best)," ","tree #",which.max(best), "\n") # report it to user

fit<-fit[[which.max(best)]] # extract pml model with best tree and use it to optimize tree model parameters

cat("optimizing paramaters","\n")
best.fit<-optim.pml(fit,optEdge=FALSE,optGamma=TRUE,optBf=FALSE,optQ=TRUE) # find best model; Q matrix and gamma

siteLik<-matrix(0,nrow=length(attr(data,"index")),ncol=length(trees)) # create a matrix to store site-wise likeliohoods for each tree
cat("getting site-wise likelihoods","\n")

for (i in 1:length(trees)){ # loop through all trees and calculate site-wise likelihoods using the best fit model from above
	x <- update(best.fit,tree=trees[[i]])
	siteLik[,i]<-rep(x$siteLik,attr(data,"weight")) # store site-wise likelihoods in columns of matrix
	siteLik}
cat("performing AU test", "\n")
AU<-relltest(siteLik) # performs multiscale bootstrap on the column elements of the matrix = site-wise likelihoods 
AU ## returns the result of the relltest function of scaleboot, must use summary() to get the p=values. 
}

baseFreq <- function(dat){ ## a function to calculate the emperical base freq of a alignment
    if (class(dat) != "phyDat") 
        stop("data must be of class phyDat")
    levels <- attr(dat,"levels")
    weight <- attr(dat,"weight")	
    n <- length(dat)	
    res <- numeric(length(levels)) 	
    for(i in 1:n)res <- res+table(rep(dat[[i]],attr(dat,"weight")))[1:4]	
    res <- res/sum(res)
    names(res) <- levels	
    res	
}
