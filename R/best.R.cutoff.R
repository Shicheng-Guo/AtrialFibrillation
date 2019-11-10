xx<-t(myNorm[na.omit(match(subset(map,map[1]=="chr22")[,4],rownames(myNorm))),])
R<-c()
P<-c()
for(i in 1:ncol(xx)){
  for(j in (i+1):ncol(xx)){
    x<-xx[,i]
    y<-xx[,j]
    fit<-cor.test(x,y)
      r=abs(fit$estimate)
      p=fit$p.value
      R<-c(R,r)
      P<-c(P,p)
    }
  }
pdf("P-R.pdf")
plot(-log(P,10),R)
dev.off()
