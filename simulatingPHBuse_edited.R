##2018-06-11. Used to generate diagrams showing how long an individual rhizobia cell could support metabolism on a given amount of PHB. 
#ms1 <-getwd() # working directory for MS1
#load starvation survival functions:
source("functionsForEstimatingSurvivalOnPHB.R")

#first, make a simple graph with survival time on the x-axis and PHB on the y axis.
# get estimates for 1st year:
# the default start date is September 15th (day 258)
phb_high <- calculatePHB_1phase(gSB=0.5,rate = "high",temp=-999,t=365*100,d1=258)
phb_med10 <- calculatePHB_1phase(gSB=0.5,rate="med",t=365*100,midphb = 0.1,d1=258,temp=-999)
phb_med100 <- calculatePHB_1phase(gSB=0.5,rate="med",t=365*100,midphb=0.01,d1=258)
phb_med1000 <- calculatePHB_1phase(gSB=0.5,rate="med",t=365*1000,midphb=0.001,d1=258)
phb_low <- calculatePHB_1phase(gSB=0.5,rate="low",t=365*2000,d1=258)

# survival for 19 months (Sep 15th year 1 to may 15th, year 3 )
daterange <- as.Date(c("2018-05-15","2016-09-15"),format=c("%Y-%m-%d"))
daterange[1] - daterange[2] # 607 days
phb_high[607]
phb_low[607]
phb_med10[607]
phb_med100[607]
phb_med1000[607]

# how long could a cell with 0.9 or 0.09 pg survive on dormancy rates?
# make a function:
survivalestimate <- function(x){
    survivaltable <- matrix(NA,nrow=1,ncol=5,dimnames=list(c(x),c("lowD","medD","highD","lowM","highM")))
    survivaltable[1,1] <- min(which(phb_low>=x))
    survivaltable[1,2] <- min(which(phb_med1000>=x))
    survivaltable[1,3] <- min(which(phb_med100>=x))
    survivaltable[1,4] <- min(which(phb_med10>=x))
    survivaltable[1,5] <- min(which(phb_high>=x))
    survivaltable
  }
survivalestimate(0.58)
survivalestimate(0.58)/365

survivalestimate(0.082) # days
survivalestimate(0.082)/365 # days
survivalestimate(1) #days
survivalestimate(1)/365 #years

survivalestimate(0.5)/365 #years

survivalestimate(0.12) #0.12pg PHB = geometric mean from low PHB strain in Fig. 1
survivalestimate(0.12)/365 # hears. 
survivalestimate(1.16) #days
survivalestimate(1.16)/30 # months
##Will's S. meliloti:
survivalestimate(0.5)# late stationary phase
survivalestimate(0.5)/365
survivalestimate(0.3)# early stationary phase
#Those estimates are with temperature oscillations in MN. what about at room temp?
phb_high_roomtemp <- calculatePHB_1phase(rate="high",t=365*30,temp=26)
phb_low_roomtemp <- calculatePHB_1phase(rate="low",t=365*100,temp=26)
phb_med10_roomtemp <- calculatePHB_1phase(rate="med",midphb=0.1,t=365*30,temp=26)
phb_med100_roomtemp <- calculatePHB_1phase(rate="med",midphb=0.01,t=365*30,temp=26)
phb_med1000_roomtemp <- calculatePHB_1phase(rate="med",midphb=0.001,t=365*30,temp=26)

min(which(phb_high_roomtemp>=0.19))
min(which(phb_low_roomtemp>=0.19))/365
min(which(phb_med100_roomtemp>=0.19))/365
min(which(phb_med1000_roomtemp>=0.19))/365
min(which(phb_med10_roomtemp>=0.19))


##Get range of measurable PHB from senescent nodules:
load("senescentNoduleDataEditedA3B1.RData")
range(nr1$phbEst[!is.na(nr1$phbEst)])
summary(nr1$phbEst)
hist(nr1$phbEst)
hist(nr1[nr1$strain=="A3","phbEst"],border="red",xlim=c(-0.5,2),ylim=c(0,6))
hist(nr1[nr1$strain=="B1","phbEst"],add=TRUE,border="blue")
# This is the range from nodule means. Ford suggested I do histograms for all cells,
# regardless of which nodules they came from. What does this look like?
# this excludes one B1 nodule that was an outlier from everything else--possibly an A3 nodule introduced 
# by contamination
# what's the average of B1 nodule means if I exclude the outlier?
summary(nr1[nr1$strain=="B1" & nr1$nodID !="3B2","phbEst"])
#use 2B1 and 2A1 as illustrative nodules:

# actually, the counts aren't that different--B1 nodules have more cells by about 100. 
# put negative values at 0 PHB:
phist_expanded2$phbEst_zeroAdjusted <- phist_expanded2$phbEst
phist_expanded2[phist_expanded2$phbEst<0,"phbEst_zeroAdjusted"] <- 0


hist(phist_expanded2[phist_expanded2$strain=="A3","phbEst"],col="black",density=10,xlim=c(0,2),freq=FALSE)
hist(phist_expanded2[phist_expanded2$strain=="B1","phbEst"],add=TRUE,density=30,freq=FALSE)

hist(phist_expanded2[phist_expanded2$nodID=="2A1","phbEst"],col="black",density=10,xlim=c(0,2),freq=FALSE)
hist(phist_expanded2[phist_expanded2$nodID=="2B1","phbEst"],add=TRUE,density=30,freq=FALSE)


#What were the nodule averages for each strain? excluding the outlier of B1?
nr2 <- nr1[nr1$file!="NRnodules.016" & !is.na(nr1$phbEst), ]
# sort by strain:
nr2 <- nr2[order(nr2$strain),]
tapply(nr2$phbEst,nr2$strain, mean)
tapply(nr2$phbEst,nr2$strain, sd)
tapply(nr2$phbEst,nr2$strain, summary)

# I should use geometric mean from flowjo 
sort(tapply(nr2$phbEst,nr2$nodID, mean))
tapply(nr2$Count,nr2$nodID, mean)
##Pick a nodule that's representative (it's within-nodule mean is close to the among-nodule means)
# A3: median is 0.56, mean is 0.53
# B1: median is 0.13, mean is 0.15
##What's the mean from individual cell data?
tapply(phist_expanded2$phbEst,phist_expanded2$strain,mean)
#Considering the phbEst has lots of zeros and negatives, better
# to get mean based on fluorescence then use standards to calibrate
tapply(phist_expanded2$FL3,phist_expanded2$strain,mean)
tapply(phist_expanded2$FL3,phist_expanded2$strain,summary)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
gmmean_bystrain <- tapply(phist_expanded2$FL3,phist_expanded2$strain,gm_mean)
#convert to estimated PHB/cell (standards = standards from senescent A3/B1 nodules)
cal_a3b1 <- lm(FL3mean ~ phbPerCellPg, data=standards)
plot(FL3mean ~ phbPerCellPg, data=standards)
abline(cal_a3b1)
meanphb_bystrain <- (gmmean_bystrain - coef(cal_a3b1)[1])/coef(cal_a3b1)[2]
meanphb_bystrain
b1_1090 <- with(phist_expanded2[phist_expanded2$strain=="B1",],quantile(phbEst_zeroAdjusted,c(0.9)))
a3_1090 <- with(phist_expanded2[phist_expanded2$strain=="A3",],quantile(phbEst_zeroAdjusted,c(0.9)))

##What % of each strain could survive 20 months of active metabolism?
lowM20months <- phb_med10[20*30]
with(phist_expanded2[phist_expanded2$strain=="A3",],sum(phbEst_zeroAdjusted>=lowM20months)/length(phbEst_zeroAdjusted))
##77% of cells in the high-PHB strain have enough to last 20 months of low maintence
with(phist_expanded2[phist_expanded2$strain=="B1",],sum(phbEst_zeroAdjusted>=lowM20months)/length(phbEst_zeroAdjusted))
## 38% for the low-PHB strain.

# compare to estimate from 90th percentile fluorescence (just to double check)
gmp90_bystrain <- tapply(phist_expanded2$FL3,phist_expanded2$strain,function(x){quantile(x,c(0.9))})
p90phbbystrain <- (gmp90_bystrain - coef(cal_a3b1)[1])/coef(cal_a3b1)[2]
# It's the same--just wanted to double check. 

phb_high <- calculatePHB_1phase(gSB=0.5,rate = "high",temp=-999,t=365*100,d1=258)
phb_med10 <- calculatePHB_1phase(gSB=0.5,rate="med",t=365*100,midphb = 0.1,d1=258,temp=-999)
phb_med100 <- calculatePHB_1phase(gSB=0.5,rate="med",t=365*2000,midphb=0.01,d1=258)
phb_med1000 <- calculatePHB_1phase(gSB=0.5,rate="med",t=365*2000,midphb=0.001,d1=258)
phb_low <- calculatePHB_1phase(gSB=0.5,rate="low",t=365*3000,d1=258)


logdays <- log(c(1:(365*100)))
logdays2 <- log(c(1:(365*2000)))
logdays3 <- log(c(1:(365*3000)))
days <- c(1:(365*100))
#tiff(file="Supp_timeUntilPHBRunsOut_Sep15.tiff",width=640,height=470,units="px",res=100)
par(yaxt="n",xpd=FALSE)
plot(phb_low,logdays3,type="l",xlim=c(0,2),xlab="pg PHB per cell",ylab= "time until PHB runs out",
     ylim=c(0,log(365*2000)),cex.lab=1.3,cex.axis=1.2,lty=2,lwd=1.5)
points(phb_med10,logdays,type="l",lwd=1.5,lty=1)
points(phb_med100,logdays2,type="l",lwd=1.5,lty=2)
points(phb_med1000,logdays2,type="l",lwd=1.5,lty=2)
points(phb_high,logdays,type="l",lwd=1.5,lty=1)
par(yaxt="s",xpd=TRUE)
axis(2,at=c(log(c(7,30,30*6,365,365*5,365*10,365*50,365*100))),labels=c("1wk","1mo","6mo","1yr","5yr","10yr","50yr","100yr"),las=1)
#denstrip(ltcornSensecenceAdjusted,at=y_lab[1],adjust=0.5,
##Check histograms:
with(phist_expanded2[phist_expanded2$strain=="A3",],hist(phbEst_zeroAdjusted))
with(phist_expanded2[phist_expanded2$strain=="B1",],hist(phbEst_zeroAdjusted,add=TRUE))
with(phist_expanded2, tapply(phbEst_zeroAdjusted,strain,summary))
phist_expanded2[phist_expanded2$phbEst_zeroAdjusted>5,]
#There's one super high cell from B1. It's probably precipitate
phist_expanded2 <- phist_expanded2[phist_expanded2$phbEst_zeroAdjusted<10,]
hbreaks <- c(0,0.1,0.25,seq(0.5,8,by=0.25))

# Make a version of this graph that has histograms from senescent soybean nodules above the main graph
#tiff(file="Fig2survivalFromSept15PlusHistogram.tiff",width=500,height=700,res=100)
nf <- layout(matrix(c(1,2),2,1),heights=c(1,1.5))
#layout.show(nf)
defaultmargins <- c(5,4,4,2)+0.1
par(mar=c(0,4.1,2.1,2.1),xaxt="n")
#hist(phist_expanded2[phist_expanded2$strain=="A3","phbEst_zeroAdjusted"],add=TRUE,col="grey10",density=30,freq=TRUE,breaks=hbreaks)
hbreaks <- c(0,0.1,0.25,seq(0.5,8,by=0.25))
h1 <- hist(phist_expanded2[phist_expanded2$strain=="B1" ,"phbEst_zeroAdjusted"],freq=TRUE,main="",breaks=hbreaks,plot=F)
h1$counts <- (h1$counts/sum(h1$counts))*100
plot(h1,freq=TRUE,xlab="",main="",ylab="% of rhizobia",cex.lab=1.3,col="grey",xlim=c(0,2.5))
h2 <- hist(phist_expanded2[phist_expanded2$strain=="A3" ,"phbEst_zeroAdjusted"],freq=TRUE,main="",breaks=hbreaks,plot=F)
h2$counts <- 100*h2$counts/sum(h2$counts)
plot(h2,freq=TRUE,main="",add=TRUE,density=30)
b1_1090 <- with(phist_expanded2[phist_expanded2$strain=="B1",],quantile(phbEst_zeroAdjusted,c(0.9)))
a3_1090 <- with(phist_expanded2[phist_expanded2$strain=="A3",],quantile(phbEst_zeroAdjusted,c(0.9)))
legend("topright",bty="n",fill=c("grey","black"),density=c(10000,30),
       legend=c("low-PHB field strain","high-PHB field strain"),cex=1.3)
text(0.32,48.6,labels="A",cex=2.5)
par(xaxt="s",yaxt="n",xpd=FALSE,mar=c(4.1,4.1,0,2.1))
plot(phb_low,logdays3,type="l",xlim=c(0,2.5),ylim=c(0,log(365*2000)),xlab="initial PHB per cell (pg)",ylab= "",cex.lab=1.3,cex.axis=1.2,lty=1,lwd=1.5)
par(xpd=TRUE)
title(ylab="time until PHB runs out",cex.lab=1.3)
par(xpd=FALSE)
points(phb_med10,logdays,type="l",lwd=1.5,lty=2)
points(phb_med100,logdays2,type="l",lwd=1.5,lty=2)
points(phb_med1000,logdays2,type="l",lwd=1.5,lty=2)
points(phb_high,logdays,type="l",lwd=1.5,lty=1)
par(yaxt="s")
abline(v=c(b1_1090,a3_1090),lwd=2,col=c("grey","black"))
abline(v=c(meanphb_bystrain),lwd=2,col=c("black","grey")) #add geometric means (high strain is first)
axis(2,at=c(log(c(7,30,30*6,365,365*5,365*10,365*50,365*100))),labels=c("1wk","1mo","6mo","1yr","5yr","10yr","50yr","100yr"),las=1)
text(0.295,12.19,labels="B",cex=2.5)
#dev.off()


##Do the same thing, but this time use just my reference nodules
# Make a version of this graph that has histograms from senescent soybean nodules above the main graph
#tiff(file="survivalFromSept15PlusHistogram_1nodperstrain.tiff",width=500,height=700,res=100)
nf <- layout(matrix(c(1,2),2,1),heights=c(1,1.5))
#layout.show(nf)
defaultmargins <- c(5,4,4,2)+0.1
par(mar=c(0,4.1,2.1,2.1),xaxt="n")
#hist(phist_expanded2[phist_expanded2$nodID=="2A1","phbEst_zeroAdjusted"],add=FALSE,col="grey10",density=30,freq=TRUE,breaks=hbreaks)
h1 <- hist(phist_expanded2[phist_expanded2$nodID=="2B1" ,"phbEst_zeroAdjusted"],freq=TRUE,main="",breaks=hbreaks,plot=F)
h1$counts <- (h1$counts/sum(h1$counts))*100
plot(h1,freq=TRUE,xlab="",main="",ylab="% rhizobia",cex.lab=1.3,col="grey",xlim=c(0,2.5))
h2 <- hist(phist_expanded2[phist_expanded2$nodID=="2A1" ,"phbEst_zeroAdjusted"],freq=TRUE,main="",breaks=hbreaks,plot=F)
h2$counts <- 100*h2$counts/sum(h2$counts)
plot(h2,freq=TRUE,main="",add=TRUE,density=30)
b1_1090 <- with(phist_expanded2[phist_expanded2$nodID=="2B1",],quantile(phbEst_zeroAdjusted,c(0.9)))
a3_1090 <- with(phist_expanded2[phist_expanded2$nodID=="2A1",],quantile(phbEst_zeroAdjusted,c(0.9)))
legend("topright",bty="n",fill=c("grey","black"),density=c(10000,30),legend=c("low PHB field strain","high PHB field strain"),cex=1.3)

par(xaxt="s",yaxt="n",xpd=FALSE,mar=c(4.1,4.1,0,2.1))
plot(phb_low,logdays2,type="l",xlim=c(0,2.5),xlab="initial PHB per cell (pg)",ylab= "",cex.lab=1.3,cex.axis=1.2,lty=2,lwd=1.5)
par(xpd=TRUE)
title(ylab="time until PHB runs out",cex.lab=1.3)
par(xpd=FALSE)
points(phb_med10,logdays,type="l",lwd=1.5,lty=1)
points(phb_med100,logdays,type="l",lwd=1.5,lty=2)
points(phb_high,logdays,type="l",lwd=1.5,lty=1)
par(yaxt="s")
abline(v=c(b1_1090,a3_1090),lwd=2,col=c("grey","black"))
axis(2,at=c(log(c(7,30,30*6,365,365*5,365*10,365*50,365*100))),labels=c("1wk","1mo","6mo","1yr","5yr","10yr","50yr","100yr"),las=1)
par(mfcol=c(1,1))

#This version assumes start day of September 15th. Make a set of graphs with different start dates:
startDates <- c(228,258,289,320,350)
startMos <- c("Aug15","Sep15","Oct15","Nov15","Dec15")
for(i in 1:length(startDates)){
phb_high <- calculatePHB_1phase(rate="high",t=365*100,d1=startDates[i])
phb_med10 <- calculatePHB_1phase(rate="med",t=365*100,midphb = 0.1,d1=startDates[i])
phb_med100 <- calculatePHB_1phase(rate="med",t=365*2000,midphb=0.01,d1=startDates[i])
phb_med1000 <- calculatePHB_1phase(rate="med",t=365*3000,midphb=0.001,d1=startDates[i])
phb_low <- calculatePHB_1phase(rate="low",t=365*3000,d1=startDates[i])

#tiff(file=paste("FigA2_timeUntilPHBRunsOut_",startMos[i],".tiff",sep=""),width=500,height=500,units="px",res=100)
par(yaxt="n",xpd=FALSE)
plot(phb_low,logdays3,type="l",xlim=c(0,2),xlab="pg PHB per cell",ylab= "time until PHB runs out",cex.lab=1.3,cex.axis=1.2,lty=1,lwd=2)
points(phb_med10,logdays,type="l",lwd=2,lty=2)
points(phb_med100,logdays2,type="l",lwd=2,lty=2)
points(phb_med1000,logdays3,type="l",lwd=2,lty=2)
points(phb_high,logdays,type="l",lwd=2,lty=1)
title(main=startMos[i])
par(yaxt="s",xpd=TRUE)
axis(2,at=c(log(c(7,30,30*6,365,365*5,365*10,365*50,365*100))),labels=c("1wk","1mo","6mo","1yr","5yr","10yr","50yr","100yr"),las=1)
#dev.off()
}

