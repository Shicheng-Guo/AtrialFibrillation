
ManhattanmyDMP<-function(myDMP){
  library("qqman")
  library("Haplin")
  SNP=rownames(myDMP)
  CHR=myDMP$CHR
  if(length(grep("X",CHR))>0){
    CHR<-sapply(CHR,function(x) gsub(pattern = "X",replacement = "23",x))
    CHR<-sapply(CHR,function(x) gsub(pattern = "Y",replacement = "24",x))
  }
  CHR<-as.numeric(CHR)
  BP=myDMP$MAPINFO
  P=myDMP$P.Value
  manhattaninput=data.frame(SNP,CHR,BP,P)
  max<-max(2-log(manhattaninput$P,10))
  genomewideline=0.05/nrow(manhattaninput)
  pdf("manhattan.pdf")
  manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,10),lwd=2, suggestiveline=F,genomewideline=FALSE)
  dev.off()
  pdf("qqplot.pdf")
  pQQ(P, nlabs =length(P), conf = 0.95)
  dev.off()
}

