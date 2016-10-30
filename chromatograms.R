x<-read.table(file="E:/omicron-web-design/files-analysis/chemometrics/chromatograms.txt")
#dim(x)
#write.csv(x,file="C:/Documents and Settings/snlo/My Documents/other/DSTS-chromatography/chromatograms.csv")

################################################################
#Modelling
################################################################
library(pls)

glyphosate<-as.numeric(x[1,2:27])
retension.time<-x[2:24001,1]
signal<-t(as.matrix(x[2:24001,2:27]))
x1<-data.frame(glyphosate=I(glyphosate),signal=I(signal))
train<-x1[1:21,]
test<-x1[22:26,]

model1<-plsr(glyphosate~signal,ncomp=19,data=train,validation="LOO")
plot(RMSEP(model1))
p1<-predict(model1,ncomp=5,newdata=test)
p2<-predict(model1,ncomp=6,newdata=test)
p3<-predict(model1,ncomp=7,newdata=test)
p4<-predict(model1,ncomp=8,newdata=test)
p5<-predict(model1,ncomp=9,newdata=test)
p6<-predict(model1,ncomp=10,newdata=test)
PLS.unscaled<-cbind(p1,p2,p3,p4,p5,p6)
rownames(PLS.unscaled)<-c("16601","9942","4517","5227","5360")
colnames(PLS.unscaled)<-c("ncomp=5","ncomp=6","ncomp=7","ncomp=8","ncomp=9","ncomp=10")

model2<-plsr(glyphosate~signal,ncomp=19,scale=TRUE,data=train,validation="LOO")
windows()
plot(RMSEP(model2))
p1<-predict(model2,ncomp=5,newdata=test)
p2<-predict(model2,ncomp=6,newdata=test)
p3<-predict(model2,ncomp=7,newdata=test)
p4<-predict(model2,ncomp=8,newdata=test)
p5<-predict(model2,ncomp=9,newdata=test)
p6<-predict(model2,ncomp=10,newdata=test)
PLS.scaled<-cbind(p1,p2,p3,p4,p5,p6)
rownames(PLS.scaled)<-c("16601","9942","4517","5227","5360")
colnames(PLS.scaled)<-c("ncomp=5","ncomp=6","ncomp=7","ncomp=8","ncomp=9","ncomp=10")

model3<-pcr(glyphosate~signal,ncomp=19,data=train,validation="LOO")
windows()
plot(RMSEP(model3))
p1<-predict(model3,ncomp=5,newdata=test)
p2<-predict(model3,ncomp=6,newdata=test)
p3<-predict(model3,ncomp=7,newdata=test)
p4<-predict(model3,ncomp=8,newdata=test)
p5<-predict(model3,ncomp=9,newdata=test)
p6<-predict(model3,ncomp=10,newdata=test)
PCR.unscaled<-cbind(p1,p2,p3,p4,p5,p6)
rownames(PCR.unscaled)<-c("16601","9942","4517","5227","5360")
colnames(PCR.unscaled)<-c("ncomp=5","ncomp=6","ncomp=7","ncomp=8","ncomp=9","ncomp=10")

model4<-pcr(glyphosate~signal,ncomp=19,scale=TRUE,data=train,validation="LOO")
windows()
plot(RMSEP(model4))
p1<-predict(model4,ncomp=5,newdata=test)
p2<-predict(model4,ncomp=6,newdata=test)
p3<-predict(model4,ncomp=7,newdata=test)
p4<-predict(model4,ncomp=8,newdata=test)
p5<-predict(model4,ncomp=9,newdata=test)
p6<-predict(model4,ncomp=10,newdata=test)
PCR.scaled<-cbind(p1,p2,p3,p4,p5,p6)
rownames(PCR.scaled)<-c("16601","9942","4517","5227","5360")
colnames(PCR.scaled)<-c("ncomp=5","ncomp=6","ncomp=7","ncomp=8","ncomp=9","ncomp=10")

PLS.unscaled
PLS.scaled
PCR.unscaled
PCR.scaled

#PLS.scaled with ncomp=6 chosen and submitted

################################################################
#Plotting chromatograms for each plant
################################################################

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred_gly0.png",width=1400,height=380)
par(mfrow=c(1,3),cex=1.0) 
plot(x[2:24001,1],x[2:24001,2],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=0")
plot(x[2:24001,1],x[2:24001,3],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=0")
plot(x[2:24001,1],x[2:24001,4],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=0")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred_gly1.png",width=1400,height=380)
par(mfrow=c(1,3),cex=1.0) 
plot(x[2:24001,1],x[2:24001,5],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=1")
plot(x[2:24001,1],x[2:24001,6],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=1")
plot(x[2:24001,1],x[2:24001,7],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=1")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred_gly5.png",width=1400,height=380)
par(mfrow=c(1,3),cex=1.0) 
plot(x[2:24001,1],x[2:24001,8],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=5")
plot(x[2:24001,1],x[2:24001,9],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=5")
plot(x[2:24001,1],x[2:24001,10],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=5")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred_gly10.png",width=1400,height=380)
par(mfrow=c(1,3),cex=1.0) 
plot(x[2:24001,1],x[2:24001,11],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=10")
plot(x[2:24001,1],x[2:24001,12],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=10")
plot(x[2:24001,1],x[2:24001,13],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=10")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred_gly20.png",width=1400,height=380)
par(mfrow=c(1,3),cex=1.0) 
plot(x[2:24001,1],x[2:24001,14],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=20")
plot(x[2:24001,1],x[2:24001,15],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=20")
plot(x[2:24001,1],x[2:24001,16],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=20")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred_gly30.png",width=1400,height=380)
par(mfrow=c(1,3),cex=1.0) 
plot(x[2:24001,1],x[2:24001,17],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=30")
plot(x[2:24001,1],x[2:24001,18],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=30")
plot(x[2:24001,1],x[2:24001,19],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=30")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred_gly50.png",width=1400,height=380)
par(mfrow=c(1,3),cex=1.0) 
plot(x[2:24001,1],x[2:24001,20],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=50")
plot(x[2:24001,1],x[2:24001,21],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=50")
plot(x[2:24001,1],x[2:24001,22],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate=50")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred1.png",width=650,height=500)
par(cex=1.3) 
plot(x[2:24001,1],x[2:24001,23],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate to be predicted")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred2.png",width=650,height=500)
par(cex=1.3)
plot(x[2:24001,1],x[2:24001,24],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate to be predicted")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred3.png",width=650,height=500)
par(cex=1.3)
plot(x[2:24001,1],x[2:24001,25],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate to be predicted")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred4.png",width=650,height=500)
par(cex=1.3)
plot(x[2:24001,1],x[2:24001,26],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate to be predicted")
dev.off()

png(file="E:/omicron-web-design/files-analysis/chemometrics/pred5.png",width=650,height=500)
par(cex=1.3)
plot(x[2:24001,1],x[2:24001,27],type="l",xlab="Retention time",ylab="Signal",main="Glyphosate to be predicted")
dev.off()

