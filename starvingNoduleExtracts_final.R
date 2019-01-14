#####################################
#starvingNoduleExtracts_final.R
###This script is an edited version of the rough-draft script starvingNoduleExtract_PHBusevsModel.R.
#Load data from starving nodule extract:
load("starvingNoduleExtractPHBuse.RDATA")
# load model functions:
source("functionsForEstimatingSurvivalOnPHB.R")
##Load viability data from dilution plating:
load("viabilityCountsEdited.RData")

###To what extent did contamination affect rhizobia colony counts?
### Look at variation among replicates (same strain and sample date). Set mean of replicates as zero, and look at distance from mean 
ld$strainsample <- with(ld, paste(strain,sample,sep="."))
strainsample <- unique(ld$strainsample)
for(i in 1:length(strainsample)){
  ld[ld$strainsample==strainsample[i],"replicatemeanCFU"] <- 
    mean(ld[ld$strainsample==strainsample[i],"meancellsperml"])
}
ld$replicatemeanminusCFU <- with(ld, replicatemeanCFU-meancellsperml)

defaultmar <- c(5,4,4,2) + 0.1
#tiff("FigE1_effectOfContamination.tiff",width=450,height=800,res=100)
#par(mfcol=c(2,1),mar=c(4.1,4.1,2.1,2.1))
layout(matrix(1:2,ncol=1),widths=1,heights=c(2,2),respect=FALSE)
par(mar=c(4.1,2.1,0,0),oma=c(4,4,2,2))
plot(replicatemeanminusCFU ~ contaminantCFUperML, data=ld,xlim=c(0,2e7),lwd=2,cex=1.5,cex.axis=1.3,
     xlab="contaminant CFU/mL",cex.lab=1.4,ylab="")
abline(h=0,lwd=2,lty=2)
par(mar=c(0,2.1,3.1,0),xpd=NA)
plot(replicatemeanminusCFU ~ quality, data=ld,xlim=c(1,5),
     ylab="",xlab="quality index",cex.axis=1.3,cex.lab=1.4,lwd=2,cex=1.5)
mtext("distance from replicate mean rhizobia CFU/mL",side=2,line=1,outer=TRUE,las=0,cex=1.4)
abline(h=0,lwd=2,lty=2,xpd=FALSE)
par(mfcol=c(1,1),mar=defaultmar,xpd=FALSE)
#dev.off()

##Visualize patterns of PHB use and reproduction for all strains (put PHB and CFU/mL on the same
# graph, keep replicates separate)

flowdatacombinedld <- flowdatacombined[flowdatacombined$treatment=="LD",]
flowdatabystrain <- flowdatabystrain[order(flowdatabystrain$starvationdays),]
cfubystrain_ld <- cfubystrain_ld[order(cfubystrain_ld$starvationdays),]
focalstrain <- "A3"
maxdays <- 455
defaultmar <- c(5,4,4,2)+0.1
focalstrain <- unique(flowdatabystrain$strain)
for(i in 1:length(focalstrain)){
#tiff(file=paste(focalstrain[i],"phbAndCFU.tiff",sep=""),width=500,height=800)
  par(mfcol=c(2,1),mar=c(2,4,2,2)+0.1)
  plot(phbEst_high2 ~ starvationdays, data=flowdatacombinedld[flowdatacombinedld$strain==focalstrain[i],],
       xlim=c(0,maxdays),ylab="PHB per cell (pg)",xlab="starvation days",cex=2,lwd=2,
       cex.lab=1.5,cex.axis=1.3)
  title(main=focalstrain[i],cex=2)
  points(meanPHB_high ~ starvationdays, data=flowdatabystrain[flowdatabystrain$strain==focalstrain[i] &
                                                                flowdatabystrain$treatment=="LD",],type="l",lwd=2)
  par(mar=c(5,4,0,2)+0.1)  
  plot(log10(meancellsperml) ~ starvationdays, data=flowdatacombinedld[flowdatacombinedld$strain==focalstrain[i],],
     xlim=c(0,maxdays),ylab="log10(CFU/mL)",xlab="starvation days",cex=2,lwd=2,cex.lab=1.5,cex.axis=1.3)
points(log10(meanCFU) ~ starvationdays, data=cfubystrain_ld[cfubystrain_ld$strain==focalstrain[i],],lwd=2,type="l")
par(mfcol=c(1,1),mar=defaultmar)
#dev.off()
}

maxdays <- 32
defaultmar <- c(5,4,4,2)+0.1
focalstrain <- unique(flowdatabystrain$strain)
for(i in 1:length(focalstrain)){
#  tiff(file=paste(focalstrain[i],"phbAndCFU30day.tiff",sep=""),width=500,height=800)
  par(mfcol=c(2,1),mar=c(2,4,2,2)+0.1)
  plot(phbEst_high2 ~ starvationdays, data=flowdatacombinedld[flowdatacombinedld$strain==focalstrain[i],],
       xlim=c(0,maxdays),ylab="PHB per cell (pg)",xlab="starvation days",cex=2,lwd=2,
       cex.lab=1.5,cex.axis=1.3)
  title(main=focalstrain[i],cex=2)
  points(meanPHB_high ~ starvationdays, data=flowdatabystrain[flowdatabystrain$strain==focalstrain[i] &
                                                                flowdatabystrain$treatment=="LD",],type="l",lwd=2)
  par(mar=c(5,4,0,2)+0.1)  
  plot(log10(meancellsperml) ~ starvationdays, data=flowdatacombinedld[flowdatacombinedld$strain==focalstrain[i],],
       xlim=c(0,maxdays),ylab="log10(CFU/mL)",xlab="starvation days",cex=2,lwd=2,cex.lab=1.5,cex.axis=1.3)
  points(log10(meanCFU) ~ starvationdays, data=cfubystrain_ld[cfubystrain_ld$strain==focalstrain[i],],lwd=2,type="l")
  par(mfcol=c(1,1),mar=defaultmar)
#  dev.off()
}



######1) How does PHB-use and reproduction differ among strains differing in initial PHB per cell? Focus on the first month only.
####################################################################################################################################

##Make 2 graphs: one showing PHB use vs. initial PHB for the first month, the other showing reproduction vs. phb use
firstmonth <- flowdatabystrain[flowdatabystrain$sample==5 & flowdatabystrain$treatment=="LD",]
#linear regression
phbuse_1moregression <- lm(phbused ~ initialPHB, data=firstmonth)
summary(phbuse_1moregression)
#correlation coefficient
with(firstmonth, cor(initialPHB, phbused))
# r = 0.963399

#tiff(file="Fig5BphbuseOver1month.tiff",width=450,height=450,res=100)
par(mar=c(5,5,3,2)+0.1)
plot(phbused ~ initialPHB, data=firstmonth,pch=19,cex=2,ylab="1 month PHB use (pg/cell)",
     xlab="initial PHB per cell (pg)",cex.lab=1.4,cex.axis=1.3)
abline(phbuse_1moregression,lwd=2)
lm_coef <- round(coef(phbuse_1moregression),3)
r2 <- round(summary(phbuse_1moregression)$r.squared,2)
eq <- substitute(paste("y = ",m,"x - ",b),list(m=lm_coef[2],b=abs(lm_coef[1])))
eqr2 <- substitute(paste("R"^{2}," = ", rsquared),list(rsquared=r2))
text(0.11, 0.3, labels=eq,pos=4,cex=1.5)
text(0.11, 0.25, labels=eqr2,pos=4,cex=1.5)
par(mar=c(5,4,4,2) + 0.1)
#dev.off()

#
#Do strains that start with more PHB also keep more PHB after slowing down?
##maybe figure out how long it takes them to get to zero..
timetozero <- flowdatabystrainld[flowdatabystrainld$sample==6,]


plot(meanPHB_high ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==5,],pch=19, cex=2,
     xlab = "initial PHB/cell (pg)",ylab="remaining PHB/cell (pg)", main = "29 starvation days")
plot(meanPHB_high ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==6,],main="127 starvation days")
plot(meanPHB_high ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==8,],main="450 starvation days")

flowdatabystrainld$percentInitialPHB <- with(flowdatabystrainld, meanPHB_high/initialPHB) * 100
flowdatabystrainld$phbused2 <- flowdatabystrainld$phbused
flowdatabystrainld[flowdatabystrainld$phbused<0, "phbused2"] <- 0
flowdatabystrainld$percentPHBused <- with(flowdatabystrainld, phbused2/initialPHB) * 100

flowdatabystrainld[flowdatabystrainld$percentInitialPHB<0,"percentInitialPHB"] <- 0
plot(percentInitialPHB ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==5,])
plot(percentInitialPHB ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==6,])
plot(percentInitialPHB ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==8,])

plot(percentPHBused ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==5,])
plot(percentPHBused ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==6,])
plot(percentPHBused ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==8,])

###Do strains that use more PHB in the first month reproduce more?
# I think it's valid to ask that, among strains that multiplied, did those that use more PHB multiply more?
phbusevsreproduction <- 
  lm(CFUchange ~ phbused, data=firstmonth[firstmonth$CFUchange>=0,])
summary(phbusevsreproduction)
confint(phbusevsreproduction)
# calculate correlation coefficient
with(firstmonth[firstmonth$CFUchange>=0,],cor(phbused,CFUchange))
## r = 0.49958
# what are the regression parameters without the outlier?
phbusevsreproductionNoOutlier <- 
  lm(CFUchange ~ phbused, data=firstmonth[firstmonth$CFUchange>=0 & firstmonth$strain !="F5",])
summary(phbusevsreproductionNoOutlier)

##the association is not statistically significant (R2=0.31, F=2.22, df=5, p=0.19)
#tiff(file="Fig3BphbUsevsReproduction.tiff",width=500,height=500,res=100)
growers <- firstmonth[firstmonth$CFUchange>=0,]
shrinkers <- firstmonth[firstmonth$CFUchange<0,]
xspan <- range(growers$phbused,shrinkers$phbused) + c(-0.01,0.01)
nf <- layout(matrix(c(1,2),2,1),heights=c(1,0.5))
par(mar=c(0,5.1,3.1,2.1))
plot(CFUchange ~ phbused, data=growers,pch=19,cex=2,ylab= "",
     xlab="",cex.lab=1.4,ylim=c(25e5,50e5),xaxt="n",yaxt="n",xlim=xspan)
axis(2,at=seq(3e6,5e6,by=0.5e6), labels=seq(3,5,by=0.5),cex.axis=1.3)
abline(phbusevsreproduction,lwd=2)
lm_coef <- round(coef(phbusevsreproduction)/10^6,2)
r2 <- round(summary(phbusevsreproduction)$r.squared,2)
eq <- substitute(paste("y = ","(",m, "*",10^6,")",
                       " x"," + ", b,"*", 10^6 ), list(m=lm_coef[2], b=lm_coef[1]))
eqr2 <- substitute(paste("R"^{2}," = ", rsquared),list(rsquared=r2))
text(0.05,4.5e6,labels=eq, cex=1.5,pos=4)
text(0.05,4.1e6,labels=eqr2, cex=1.5,pos=4)
with(growers, text(phbused+0.02, CFUchange, labels=strain))
par(mar=c(5.1,5.1,1,2.1))
plot(range(firstmonth$phbused), c(-4.5e5,-3e5),type="n",xaxt="s",ylab="",xlab="PHB use (pg/cell)",
     cex.lab=1.4,cex.axis=1.3,yaxt="n",xlim=xspan)
points(CFUchange ~ phbused, data=shrinkers, cex=2,pch=19)
axis(2,at=c(-3.5e5,-4.5e5),cex.axis=1.4,labels=c(-3.5, -4.5))
mtext(side=2,text= expression("net population growth "~(10^6~ "cells/mL"))
      ,xpd=NA ,line=2.5,at=-4,cex=1.4)
layout(1)
par(mar=defaultmar,yaxt="s",xaxt="s",xpd=FALSE)
#dev.off()


#####relative population size vs. initial PHB:
# look at the relationship between relative population size and initial PHB over time:
popgrowth <- lm(relativepop ~ initialPHB*starvationdays, data=flowdatacombinedld)
summary(popgrowth)
##37 samples are missing CFU data. I should double check that the data are actually missing:
missingcfu <- flowdatacombinedld[is.na(flowdatacombinedld$meancellsperml),"uniqueID"]
ld$uniqueID <- with(ld, paste(strain,rep,sample,"LD",sep="."))
ld[ld$uniqueID %in% missingcfu,] #Yes, they're missing (probably due to contamination). One has a count, but
# no information for row or volume. I could try to find it, but it's not worth it. 

#look at relative population size over time, use different symbols for different strains
##I think it's clearer to present it as a series of graphs:
library(ggplot2)
ggplot(data=flowdatabystrainld,aes(starvationdays,relativepop)) + geom_point(aes(color=initialPHB)) +
  geom_line(aes(group=strain,color=initialPHB))
plot(relativepop ~ starvationdays, data=flowdatabystrainld)

## 1 month

#tiff(file="relativepopvsinitialPHB_1mo.tiff",width=450,height=450,res=100)
plot(relativepop ~ initialPHB, data=firstmonth,pch=19,cex=2,ylab= "1 month relative population",xlab="initial PHB per cell (pg)",
     cex.lab=1.4,cex.axis=1.3)
abline(h=1,lty=2,lwd=2)
#dev.off()

## 5 months
#tiff(file="relativepopvsinitialPHB_5mo.tiff",width=450,height=450,res=100)
plot(relativepop ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==6,],pch=19,cex=2,ylab= "4 month relative population",xlab="initial PHB per cell (pg)",
     cex.lab=1.4,cex.axis=1.3)
abline(h=1,lty=2,lwd=2)
#dev.off()
## 15month
#tiff(file="relativepopvsinitialPHB_15mo.tiff",width=450,height=450,res=100)
plot(relativepop ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==8,],pch=19,cex=2,ylab= "15 month relative population",xlab="initial PHB per cell (pg)",
     cex.lab=1.4,cex.axis=1.3)
abline(h=1,lty=2,lwd=2)

#dev.off()

# combine 5 and 15 month into a single graph:
fourmo <- flowdatabystrainld[flowdatabystrainld$sample==6 & flowdatabystrainld$strain!="F2",c("initialPHB","relativepop")]
fifteenmo <- flowdatabystrainld[flowdatabystrainld$sample==8 &flowdatabystrainld$strain!="F2",c("initialPHB","relativepop")]

#tiff("Fig3ArelativepopvsinitialPHB_4and15mo.tiff",width=500,height=500,res=100)
par(mar=c(5.1,5.1,3.1,2.1))
plot(relativepop ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==6,],pch=19,cex=2,
     ylab= expression("relative population "~ (N[t] / N[0])),xlab="initial PHB per cell (pg)",
     cex.lab=1.4,cex.axis=1.3, ylim=c(0,27))
points(relativepop ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==8,],pch=1,cex=2,lwd=2)
abline(h=1,lty=2,lwd=2)
xx = fourmo$initialPHB
y_0 = fourmo$relativepop
y_1 = fifteenmo$relativepop
arrows(xx,y_0, xx, (y_1+0.2),length=0.1,lwd=2,angle=20,col=grey(0.5))
legend("topleft",pt.cex=2,pch=c(19,1),pt.lwd=2,legend=c("127 days","450 days"),title="time in starvation culture")
par(mar=defaultmar)
#dev.off()

# part of resilience is the population decline.
# look at decline between months 4 and 15 vs. growth between months 1 and 4:
resilience <- flowdatabystrainld[flowdatabystrainld$starvationdays %in% c(0,29,127,450),c("strain","starvationdays","meanCFU")]
resilience.w <- reshape(resilience,timevar = "starvationdays",idvar="strain",direction="wide")
resilience.w$d0to29 <- resilience.w$meanCFU.29-resilience.w$meanCFU.0
resilience.w$d29to127 <- resilience.w$meanCFU.127-resilience.w$meanCFU.29
resilience.w$d127to450 <- resilience.w$meanCFU.450-resilience.w$meanCFU.127
plot(resilience.w$d29to127,-1*resilience.w$d127to450,xlab="growth from month 1 to 4",ylab="decline from month 4 to 15")
plot(resilience.w$d0to29,resilience.w$d29to127,xlab="growth from month 0 to 1",ylab="growth from month 1 to 4")

##Does starting with more PHB mean keeping more PHB after entering dormancy?
# PHB at one month:
plot(meanPHB_high ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==5,],ylim=c(-0.03,0.2))
plot(meanPHB_high ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==6,],pch=19)
plot(meanPHB_high ~ initialPHB, data=flowdatabystrainld[flowdatabystrainld$sample==8,],pch=2)
flowdatabystrainld <- flowdatabystrainld[order(flowdatabystrainld$sample,flowdatabystrainld$strain),]
plot(flowdatabystrainld[flowdatabystrainld$sample==6,"meanPHB_high"]
  ,flowdatabystrainld[flowdatabystrainld$sample==8,"relativepop"],ylab="15 month relative pop",
  xlab="high group PHB at 4 months")

##############################################################################################################
###2) Look at change in PHB use over time (Daily PHB use rate)
##Calculate dailyPHB use as change in PHB divided by days elapsed:
samples <- unique(flowdatabystrainld$sample)
##Skip sample # 3. For most strains, this sample (3 week) showed PHB levels lower than the subsequent sample.
# there might have been something wonky with my standards. I also need to exclude sample 7 because it only included 2 strains
# (A3 and 110). I'll calculate PHB use for those separately.
samples <- samples[!samples%in% c("3","7")]
for(i in 1:(length(samples)-1)){
  phbfinal <- flowdatabystrainld[flowdatabystrainld$sample==samples[i+1],"meanPHB_high"]
  phbinitial <- flowdatabystrainld[flowdatabystrainld$sample==samples[i],"meanPHB_high"]
  dayfinal <- flowdatabystrainld[flowdatabystrainld$sample==samples[i+1],"starvationdays"]
  dayinitial <- flowdatabystrainld[flowdatabystrainld$sample==samples[i],"starvationdays"]
  flowdatabystrainld[flowdatabystrainld$sample==samples[i],"dailyPHBuse"] <- 
    (phbinitial - phbfinal)/(dayfinal-dayinitial+1)
}

## Include sample 7 for USDA110 and A3:
samples <- unique(flowdatabystrainld$sample)
samples <- samples[samples!="3"]
for(i in 1:(length(samples)-1)){
  phbfinal <- flowdatabystrainld[flowdatabystrainld$strain %in% c("A3","110")& flowdatabystrainld$sample==samples[i+1],"meanPHB_high"]
  phbinitial <- flowdatabystrainld[flowdatabystrainld$strain %in% c("A3","110")&flowdatabystrainld$sample==samples[i],"meanPHB_high"]
  dayfinal <- flowdatabystrainld[flowdatabystrainld$strain %in% c("A3","110")&flowdatabystrainld$sample==samples[i+1],"starvationdays"]
  dayinitial <- flowdatabystrainld[flowdatabystrainld$strain %in% c("A3","110")&flowdatabystrainld$sample==samples[i],"starvationdays"]
  flowdatabystrainld[flowdatabystrainld$strain %in% c("A3","110")&flowdatabystrainld$sample==samples[i],"dailyPHBuse"] <- 
    (phbinitial-phbfinal)/(dayfinal-dayinitial+1)
}

##Do the same thing for high density treatments (this will be A3 and 110 only)
flowdatabystrainhd <- flowdatabystrain[flowdatabystrain$treatment=="HD" & flowdatabystrain$strain %in% c("A3","110"),]
samples <- unique(flowdatabystrainhd$sample)
for(i in 1:(length(samples)-1)){
  phbfinal <- flowdatabystrainhd[flowdatabystrainhd$sample==samples[i+1],"meanPHB_high"]
  phbinitial <- flowdatabystrainhd[flowdatabystrainhd$sample==samples[i],"meanPHB_high"]
  dayfinal <- flowdatabystrainhd[flowdatabystrainhd$sample==samples[i+1],"starvationdays"]
  dayinitial <- flowdatabystrainhd[flowdatabystrainhd$sample==samples[i],"starvationdays"]
  flowdatabystrainhd[flowdatabystrainhd$sample==samples[i],"dailyPHBuse"] <- 
    (phbinitial - phbfinal)/(dayfinal-dayinitial+1)
}

##make a graph showing PHB use-rate over time for all strains:
flowdatabystrainld <- flowdatabystrainld[order(flowdatabystrainld$strain,flowdatabystrainld$starvationdays),]
#tiff(file="FigE3_dailyPHBuse_allstrains.tiff",width=481,height=387, res=100)
p <- ggplot(flowdatabystrainld[!flowdatabystrainld$starvationdays %in% c(15,430,450),],aes(starvationdays,dailyPHBuse)) 
p <- p + geom_point(size=2) + geom_line(aes(group=strain),size=0.8)
p <- p + theme(axis.title=element_text(size=rel(1.4)))
p <- p + theme(axis.text = element_text(size=rel(1.3)))
p
#dev.off()

#tiff(file="PHBovertime_allstrains.tiff",width=500,height=400, res=100)
strainz <- unique(flowdatabystrainld$strain)
nf <- layout(matrix(c(1,2),1,2),width=c(2,1))
par(mar=c(5.1,4.1,4.1,0))
plot(meanPHB_high ~ starvationdays, data=flowdatabystrainld,type="n",xlim=c(0,126),cex.axis=1.3,cex.lab=1.4,
     ylab="PHB per cell (pg)",xlab="starvation days")
for(i in 1:length(strainz)){
  points(meanPHB_high ~ starvationdays, data=flowdatabystrainld[flowdatabystrainld$strain==strainz[i],],
         type="b",cex=2,lwd=2)
}
par(mar=c(5.1,0.5,4.1,2.1))
plot(meanPHB_high ~ starvationdays, data=flowdatabystrainld,type="n",xlim=c(150,451),yaxt="n",cex.axis=1.3,xlab="")
for(i in 1:length(strainz)){
  points(meanPHB_high ~ starvationdays, data=flowdatabystrainld[flowdatabystrainld$strain==strainz[i],],
         type="b",cex=2,lwd=2)
}
par(mar=defaultmar,yaxt="s")
layout(1)
#dev.off()


#tiff(file="CFUovertime_allstrains.tiff",width=500,height=400, res=100)
strainz <- unique(flowdatabystrainld$strain)
nf <- layout(matrix(c(1,2),1,2),width=c(2,1))
par(mar=c(5.1,4.1,4.1,0))
plot(log10(meanCFU) ~ starvationdays, data=flowdatabystrainld,type="n",xlim=c(0,126),cex.axis=1.3,cex.lab=1.4,
     ylab="log10(CFU/mL)",xlab="starvation days")
for(i in 1:length(strainz)){
  points(log10(meanCFU) ~ starvationdays, data=flowdatabystrainld[flowdatabystrainld$strain==strainz[i],],
         type="b",cex=2,lwd=2)
}
par(mar=c(5.1,0.5,4.1,2.1))
plot(log10(meanCFU) ~ starvationdays, data=flowdatabystrainld,type="n",xlim=c(150,451),yaxt="n",cex.axis=1.3,xlab="")
for(i in 1:length(strainz)){
  points(log10(meanCFU) ~ starvationdays, data=flowdatabystrainld[flowdatabystrainld$strain==strainz[i],],
         type="b",cex=2,lwd=2)
}
par(mar=defaultmar,yaxt="s")
layout(1)
#dev.off()


#tiff(file="relativepopovertime_allstrains.tiff",width=500,height=400, res=100)
strainz <- unique(flowdatabystrainld$strain)
nf <- layout(matrix(c(1,2),1,2),width=c(2,1))
par(mar=c(5.1,4.1,4.1,0))
plot(relativepop ~ starvationdays, data=flowdatabystrainld,type="n",xlim=c(0,126),cex.axis=1.3,cex.lab=1.4,
     ylab="relative population",xlab="starvation days")
for(i in 1:length(strainz)){
  points(relativepop ~ starvationdays, data=flowdatabystrainld[flowdatabystrainld$strain==strainz[i],],
         type="b",cex=2,lwd=2)
}
par(mar=c(5.1,0.5,4.1,2.1))
plot(relativepop ~ starvationdays, data=flowdatabystrainld,type="n",xlim=c(150,451),yaxt="n",cex.axis=1.3,xlab="")
for(i in 1:length(strainz)){
  points(relativepop ~ starvationdays, data=flowdatabystrainld[flowdatabystrainld$strain==strainz[i],],
         type="b",cex=2,lwd=2)
}
par(mar=defaultmar,yaxt="s")
layout(1)
#dev.off()


##make a graph showing PHB use rate over time for USDA110 and A3
a3110 <- flowdatabystrainld[flowdatabystrainld$strain %in% c("A3","110"),]
#calculate theoretical PHB use rates from model (assuming constant temperature of 26 deg C)
ceilingmaintenance <- calculatePHB_1phase(rate="high",temp=26,t=1) # ceiling maintenance rate (from Price 2004)
lowmaintenance <- calculatePHB_1phase(rate="med",midphb=0.1,temp=26,t=1) # 10% of ceiling
highdormancy <- calculatePHB_1phase(rate="med",midphb=0.01,temp=26,t=1) # 1% of ceiling
meddormancy <- calculatePHB_1phase(rate="med",midphb=0.001,temp=26,t=1) # 0.1% of ceiling
lowdormancy <- calculatePHB_1phase(rate="low",temp=26,t=1) # floor dormancy rate (from Price 2004)

#tiff(file="Fig5AphbUseA3110.tiff",width=450,height=400,res=100)
par(mar=c(5.1,5.1,2.1,2.1))
plot(log10(dailyPHBuse) ~ starvationdays, data=a3110[a3110$starvationdays<400,],
     type="n",ylab="PHB use rate (pg/day)",cex.lab=1.4,yaxt="n",cex.axis=1.3,ylim=c(-5,-1), xlab="starvation days")
##draw horizontal lines for theoretical PHB-use rates:
abline(h=log10(c(ceilingmaintenance,lowdormancy)),lwd=2,col="grey",lty=1)
abline(h=log10(c(highdormancy,meddormancy,lowmaintenance)),lwd=2,col="grey",lty=2)
points(log10(dailyPHBuse) ~ starvationdays, data=a3110[a3110$strain=="A3" & !a3110$starvationdays%in% c(15,430,450),],
       type="b",pch=2,cex=2,lwd=2)
points(log10(dailyPHBuse) ~ starvationdays, data=a3110[a3110$strain=="110" & !a3110$starvationdays %in% c(15,430,450),],
                                                       type="b",pch=1,cex=2,lwd=2)
par(yaxt="s")
labelsY = parse(text=paste(c("10"),c(-4,-3,-2,-1),sep="^"))
axis(2,at=c(-4,-3,-2,-1),labels=labelsY,las=0,cex.axis=1.3)
# text(c(20,20),c(log10(c(ceilingmaintenance,lowdormancy))+0.1),
#      labels=c("full somatic maintenance","molecular repair only"),cex=1.1,col="grey")
par(mar=defaultmar)
#dev.off()

##make a graph showing absolute change in PHB per cell over time for USDA 110 and A3 (with model estimates)
### the point of this graph is to demonstrate that the high-PHB field isolate (A3) appears to enter dormancy with more PHB/cell
# make a vector of predicted PHB use over 422 days (day 29 to day 450 of starvation)--assuming 26 degC
phbUse_highdormant <- calculatePHB_1phase(rate="med",t=422,midphb=0.01,temp=26)
phbUse_lowdormant <- calculatePHB_1phase(rate="low",t=422,temp=26)
phbUse_meddormant <- calculatePHB_1phase(rate="med",t=422,midphb=0.001,temp=26)
phbUse_highdormant <- phbUse_highdormant[c(29,127,450)-29+1] # subset to get PHB used at days 29, 127, and 450 (setting day 29 as starting point)
phbUse_lowdormant <- phbUse_lowdormant[c(29,127,450)-29+1]
phbUse_meddormant <- phbUse_meddormant[c(29,127,450)-29+1]

#tiff(file="starvingNoduleExtract_phbUseVsPredicted.tiff",width=450,height=450,res=100)
a3initial <- a3110[a3110$strain=="A3" & a3110$starvationdays==29,"meanPHB_high"] # use the 1 month sample as the starting point
# for dormant metabolism
initial110 <- a3110[a3110$strain=="110" & a3110$starvationdays==29,"meanPHB_high"] 
#make vectorsfor drawing a shaded region of predicted PHB under dormant metabolism (top=high dormancy, bottom= low dormancy)
predX <- c(29,127,450,450,127,29) 
phbusevector <- c(phbUse_highdormant,rev(phbUse_lowdormant))
par(mar=c(5.1,4.1,2.1,2.1))
plot(meanPHB_high ~ starvationdays, data=a3110[a3110$sample!=8,],type="n",ylim=c(0,0.5),ylab="estimated PHB/cell (pg)",xlab="starvation days",cex.lab=1.3,cex.axis=1.4)
#make x and y points to draw a shaded polygon showing model predictions
predYa3 <-a3initial-phbusevector
predY110 <- initial110 - phbusevector
polygon(predX,predYa3,col=gray(0.8,alpha=0.5),border=NA)
polygon(predX,predY110,col=grey(0.5,alpha=0.5),border=NA)
points(meanPHB_high ~ starvationdays, data=a3110[a3110$strain=="A3" & a3110$sample !=8,],type="b",lwd=2,pch=2,cex=2)
points(meanPHB_high ~ starvationdays, data=a3110[a3110$strain=="110" & a3110$sample != 8,],type="b",lwd=2,pch=1,cex=2)
points(c(29,127,450),c(a3initial-phbUse_meddormant),type="l",lty=2,col=grey(0.3),lwd=2)
points(c(29,127,450),c(initial110-phbUse_meddormant),type="l",lty=2,col=grey(0.3),lwd=2)
#dev.off()

##Make a graph showing how population size changes over the same time period:
#tiff(file="starvingNoduleExtract_logCFU.tiff",width=450,height=450,res=100)
par(mar=c(5.1,4.1,2.1,2.1))
plot(log10(meanCFU) ~ starvationdays, data=a3110[a3110$sample!=8,],type="n",
     xlab="starvation days",cex.lab=1.3,cex.axis=1.4,ylab="log10(CFU/mL)")
points(log10(meanCFU) ~ starvationdays, data=a3110[a3110$strain=="A3" & a3110$sample!=8,],type="b",pch=2,lwd=2,cex=2)
points(log10(meanCFU) ~ starvationdays, data=a3110[a3110$strain=="110" & a3110$sample!=8,],type="b",pch=1,lwd=2,cex=2)
par(mar=defaultmar)
#dev.off()

###3) compare high- and low-density starvation treatments (USDA110 and A3 only)
## graph of absolute level of PHB per cell over time

#tiff(file="starvingNoduleExtract_HighVsLowdensityPHB.tiff",width=450,height=450,res=100)
plot(meanPHB_high ~ starvationdays, data=a3110[a3110$starvationdays<50,],type="n",ylim=c(0,0.5),ylab="PHB/cell (pg)",xlab="starvation days",cex.lab=1.3,cex.axis=1.4)
points(c(25,25),c(0.48,0.4),cex=2,pch=c(19,1),lwd=2)
points(c(23,23),c(0.48,0.4),cex=2,pch=c(17,2),lwd=2)
segments(c(18,18),c(0.48,0.4),c(22,22),c(0.48,0.4),lty=c(2,1),lwd=2)
text(c(28,28),c(0.48,0.4),c("HD","LD"),cex=1.5)
points(meanPHB_high ~ starvationdays, data=a3110[a3110$strain=="A3"  & a3110$starvationdays<50,],type="b",lwd=2,pch=2,cex=2)
points(meanPHB_high ~ starvationdays, data=flowdatabystrainhd[flowdatabystrainhd$strain=="A3" & flowdatabystrainhd$starvationdays<50,],type="b",lwd=2,pch=17,cex=2,lty="dashed")
points(meanPHB_high ~ starvationdays, data=a3110[a3110$strain=="110" & a3110$starvationdays<50,],type="b",lwd=2,pch=1,cex=2)
points(meanPHB_high ~ starvationdays, data=flowdatabystrainhd[flowdatabystrainhd$strain=="110" & flowdatabystrainhd$starvationdays<50,],type="b",lwd=2,pch=19,cex=2,lty="dashed")
#dev.off()


## graph of viable population size over time
#tiff(file="starvingNoduleExtract_CFUhdvsld.tiff",width=450,height=450,res=100)
plot(log10(meanCFU) ~ starvationdays, data=flowdatabystrainhd,type="n",
     xlab="starvation days",cex.lab=1.3,cex.axis=1.4,ylab="population(log10 CFU/mL)",xlim=c(0,30),ylim=c(5,7.5))
points(log10(meanCFU) ~ starvationdays, data=a3110[a3110$strain=="A3"  & a3110$starvationdays<50,],type="b",pch=2,lwd=2,cex=2)
points(log10(meanCFU) ~ starvationdays, data=a3110[a3110$strain=="110"  & a3110$starvationdays<50,],type="b",pch=1,lwd=2,cex=2)
points(log10(meanCFU)~ starvationdays, data=cfubystrain_hd[cfubystrain_hd$strain=="A3" & !cfubystrain_hd$starvationdays %in% c(21,449),],type="b",pch=17,lwd=2,cex=2,lty=2)
points(log10(meanCFU)~ starvationdays, data=cfubystrain_hd[cfubystrain_hd$strain=="110" & !cfubystrain_hd$starvationdays%in% c(21,449),],type="b",pch=19,lwd=2,cex=2,lty=2)
#dev.off()


## graph of % cells in high-PHB group.
a3110$percentHighPHB <- a3110$meanpropHighcells*100
flowdatabystrainhd$percentHighPHB <- flowdatabystrainhd$meanpropHighcells*100
##initial time point is 100%:
a3110[a3110$sample==1,"percentHighPHB"] <- 100
flowdatabystrainhd[flowdatabystrainhd$sample==1,"percentHighPHB"] <- 100

# combine graph with % cells in high-PHB subgroup vs. population size over time--just use A3 so it's not so confusing:

plot.pch <- list(LD = 2, HD=17)
plot.col <- list(phigh = "black",cfu = "grey")
plot.lty <- list(phigh=1, cfu=2)
plot.specs <- list(pt.cex=2,lwd=2,xlim = c(0,30),cex.axis=1.4, cex.lab=1.3,xlab="starvation days",ylab1="% high-PHB cells",ylab2="population (log10 CFU/mL)")

#tiff(file="Fig4B_A3only.tiff",width=500,height=500,res=100)
f4bmar = par(mar=c(5,4,4,4)+0.1)
plot(percentHighPHB ~ starvationdays, data=a3110[a3110$starvationdays<50,],type="n",ylim=c(0,100),xlab="",ylab=plot.specs$ylab1,
     xlim=plot.specs$xlim,xaxt="n",cex.axis=plot.specs$cex.axis,cex.lab=plot.specs$cex.lab)
points(percentHighPHB ~ starvationdays, data=a3110[a3110$strain=="A3"  & a3110$starvationdays<50,],type="b",lwd=plot.specs$lwd,pch=plot.pch$LD,
       cex=plot.specs$pt.cex,col=plot.col$phigh)
points(percentHighPHB ~ starvationdays, data=flowdatabystrainhd[flowdatabystrainhd$strain=="A3" & flowdatabystrainhd$starvationdays<50,],type="b",
       lwd=plot.specs$lwd,pch=plot.pch$HD,cex=plot.specs$pt.cex,lty=plot.lty$phigh,col=plot.col$phigh)
par(new=T)
plot(log10(meanCFU) ~ starvationdays, data=flowdatabystrainhd,type="n",
     xlab="starvation days",cex.lab= plot.specs$cex.lab ,cex.axis=plot.specs$cex.axis,ylab="",xlim=plot.specs$xlim,ylim=c(5,7.5),yaxt="n")
points(log10(meanCFU) ~ starvationdays, data=a3110[a3110$strain=="A3"  & a3110$starvationdays<50,],type="b",pch=plot.pch$LD,
       lwd=plot.specs$lwd, cex=plot.specs$pt.cex,lty=plot.lty$cfu,col=plot.col$cfu)
points(log10(meanCFU)~ starvationdays, data=cfubystrain_hd[cfubystrain_hd$strain=="A3" & !cfubystrain_hd$starvationdays %in% c(21,449),],type="b",
       pch=plot.pch$HD,lwd=plot.specs$lwd,cex=plot.specs$pt.cex,lty=plot.lty$cfu,col=plot.col$cfu)
axis(side=4, labels=TRUE,cex.axis=plot.specs$cex.axis,col.axis=plot.col$cfu)
mtext(side = 4,text = plot.specs$ylab2,cex = plot.specs$cex.lab,line = 2,col=plot.col$cfu)
par(new=F)
f4bmar
#dev.off()

#tiff(file="starvingNoduleExtract_propHighLDvHD.tiff",width=450,height=450,res=100)
plot(percentHighPHB ~ starvationdays, data=a3110[a3110$starvationdays<50,],type="n",ylim=c(0,100),xlab="starvation days",ylab="% high-PHB cells",cex.lab=1.3,cex.axis=1.4)
points(percentHighPHB ~ starvationdays, data=a3110[a3110$strain=="A3"  & a3110$starvationdays<50,],type="b",lwd=2,pch=2,cex=2)
points(percentHighPHB ~ starvationdays, data=flowdatabystrainhd[flowdatabystrainhd$strain=="A3" & flowdatabystrainhd$starvationdays<50,],type="b",lwd=2,pch=17,cex=2,lty="dashed")
points(percentHighPHB ~ starvationdays, data=a3110[a3110$strain=="110" & a3110$starvationdays<50,],type="b",lwd=2,pch=1,cex=2)
points(percentHighPHB ~ starvationdays, data=flowdatabystrainhd[flowdatabystrainhd$strain=="110" & flowdatabystrainhd$starvationdays<50,],type="b",lwd=2,pch=19,cex=2,lty="dashed")
#dev.off()

##Look over entire range, and combine with plot of growth:
plot(percentHighPHB ~ starvationdays, data=a3110,type="n",ylim=c(0,100),xlab="starvation days",ylab="% high-PHB cells",cex.lab=1.3,cex.axis=1.4)
points(percentHighPHB ~ starvationdays, data=a3110[a3110$strain=="A3",],type="b",lwd=2,pch=2,cex=2)
points(percentHighPHB ~ starvationdays, data=a3110[a3110$strain=="110",],type="b",lwd=2,pch=1,cex=2)
par(new=TRUE)
plot(relativepop ~ starvationdays, data=a3110, type="n",yaxt="n",xaxt="n",xlab="",ylab="")
points(relativepop ~ starvationdays, data=a3110[a3110$strain=="A3",],cex=2,pch=2,lwd=2,col="grey",type="b")
points(relativepop ~ starvationdays, data=a3110[a3110$strain=="110",],cex=2,pch=1,lwd=2,col="grey",type="b")
par(new=FALSE)

#look at other strains too:
plot(meanpropHighcells ~ starvationdays, data=flowdatabystrain[flowdatabystrain$strain=="A3",],type="b",ylim=c(0,1))

for(i in 1:9){
points(meanpropHighcells ~ starvationdays, data=flowdatabystrain[flowdatabystrain$strain==strainz[i] &                                                                    flowdatabystrain$treatment=="LD",],type="b",col=rainbow(9)[i])
}

plot(relativepop ~ starvationdays, data=flowdatabystrain[flowdatabystrain$strain=="A3",],type="b",ylim=c(0,25))
for(i in 1:9){
  points(relativepop ~ starvationdays, data=flowdatabystrain[flowdatabystrain$strain==strainz[i] &                                                                    flowdatabystrain$treatment=="LD",],type="b",col=rainbow(9)[i])
}


with(flowdatabystrain[flowdatabystrain$sample==8,], plot(maxCFU/initialCFU,relativepop,ylab="15 month relative pop",xlab="maximum pop growth"))
for(i in 1:length(strainz)){
  flowdatabystrain[flowdatabystrain$strain==strainz[i] & flowdatabystrain$treatment=="LD","oneMonthrelativepop"] <-
    flowdatabystrain[flowdatabystrain$strain==strainz[i] & flowdatabystrain$treatment=="LD"
                     & flowdatabystrain$sample==6,"relativepop"]
}

plot(relativepop ~ starvationdays, data=flowdatabystrain[flowdatabystrain$strain=="A3",],type="b",ylim=c(0,25))
for(i in 1:9){
  points(relativepop ~ starvationdays, data=flowdatabystrain[flowdatabystrain$strain==strainz[i] &                                                                    flowdatabystrain$treatment=="LD",],type="b",col=rainbow(9)[i])
}
with(flowdatabystrain[flowdatabystrain$sample==8,], plot(oneMonthrelativepop,relativepop,ylab="15 month relative pop",xlab="one month pop growth"))

###############################################################
####3) Live-dead staining
# what % of gated cells absorbed propidium iodide?
# unfixed, no PHB stain:
livedead <- livedead[order(livedead$stain,livedead$sampleID),]
livedead$propdead <- with(livedead, countDead/countGated)
livedead[,c("stain","sampleID","propdead")]
mean_narm <- function(x){mean(x,na.rm=TRUE)}

with(livedead,tapply(propdead,list(stain,strain),function(x){mean(x,na.rm=TRUE)}))
with(livedead,tapply(propdead,list(stain,strain),function(x){min(x,na.rm=TRUE)}))
with(livedead,tapply(propdead,list(stain,strain),function(x){max(x,na.rm=TRUE)}))
livedeadplusphb$propdeadHigh <- livedeadplusphb$countDeadHigh/livedeadplusphb$countGatedHigh
livedeadplusphb[,c("propdeadHigh","countGatedHigh","sampleID","stain")]
##The problem here is that the PI staining without Bodipy includes more debris in the gated cell count,
# so the % dead cells is lower. What I really need to do is compare PI before and after fixing with Bodipy. 
## 


##############################
# How does initial PHB influence how much is left after a long starvation time?
plotPHB <- function(samp){
  sdays <- unique(flowdatabystrain[flowdatabystrain$sample==samp,"starvationdays"])
  plot(meanPHB_high ~ initialPHB, data=flowdatabystrain[flowdatabystrain$sample==samp,], main=paste(sdays,"days",sep=" "))
  points(meanPHB_high ~ initialPHB, data=flowdatabystrain[flowdatabystrain$sample==samp & flowdatabystrain$treatment=="HD",],pch=19,col="red")
  abline(h=0)
}
plotPHB("2")
plotPHB("3")
plotPHB("4")
plotPHB("5")
plotPHB("6")
plotPHB("8")

# Just out of curiosity, I want to see how my PHB measurements from the same strains in a different set of plants (SS3)
# compared to PHB used in plants going into the starvation experiment
load("SS3regatedEditedPHBdata.RData")
was14_starvingstrains <- ss3phbdata[ss3phbdata$strain %in% strainz,]
ss3means <- tapply(was14_starvingstrains$phbEst,was14_starvingstrains$strain,mean)
for(i in 1:length(ss3means)){
  flowdatabystrain[flowdatabystrain$strain==names(ss3means)[i],"meanPHBfromSS3"] <- ss3means[i]
}
plot(initialPHB ~ meanPHBfromSS3, data=flowdatabystrain[flowdatabystrain$sample==1 & flowdatabystrain$treatment=="LD",],ylim=c(0,0.45),xlim=c(0,0.45))
abline(0,1)
# most were pretty close, but some were a lot higher in the sample that went in starvation buffer. 
plot(relativepop ~ meanPHBfromSS3, data=flowdatabystrain[flowdatabystrain$sample==5,])
plot(relativepop ~ meanPHBfromSS3, data=flowdatabystrain[flowdatabystrain$sample==6,])
plot(relativepop ~ meanPHBfromSS3, data=flowdatabystrain[flowdatabystrain$sample==8,])

ggplot(flowdatabystrain[flowdatabystrain$treatment=="LD",],aes(starvationdays,phbused/initialPHB)) + geom_point(aes(col=strain)) + geom_line(aes(col=strain))


##see how PHB use for all strains compares to theoretical estimates:
firstmonth$phbused/29

flowdatabystrainld$propCeilingPHB <- flowdatabystrainld$dailyPHBuse/ceilingmaintenance
flowdatabystrainld$propFloorPHB <- flowdatabystrainld$dailyPHBuse/lowdormancy

ggplot(flowdatabystrainld, aes(x=starvationdays,y=propCeilingPHB)) + geom_point(aes(color=strain))

s1phbuse <- flowdatabystrainld[flowdatabystrainld$sample==1,c("strain","initialPHB","phbused","dailyPHBuse","propCeilingPHB")]
s2phbuse <- flowdatabystrainld[flowdatabystrainld$sample==2,c("strain","initialPHB","phbused","dailyPHBuse","propCeilingPHB")]
s5phbuse <- flowdatabystrainld[flowdatabystrainld$sample==5,c("strain","initialPHB","phbused","dailyPHBuse","propCeilingPHB","propFloorPHB")]
s6phbuse <- flowdatabystrainld[flowdatabystrainld$sample==6,c("strain","initialPHB","phbused","dailyPHBuse","propCeilingPHB","propFloorPHB")]
summary(s1phbuse$propCeilingPHB)
# Strain A4 had an increase in estimated mean PHB from sample 5 to sample 5--probably measurement error. 
s5phbuse

