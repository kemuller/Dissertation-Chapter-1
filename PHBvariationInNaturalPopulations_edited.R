#####################################
#2018-06-11,
# Script for analyzing and presenting data on PHB variation in sympatric rhizobia populations. 
library(plyr)
library(ggplot2)
library(lme4)
library(sjstats)
library(lattice)
#################################
# Chamaecrista nodules collected from Greycloud dunes:
load("chamecristaCombinedNodData_edited.RDAtA")
cfas <- nod.data.combined
rm(nod.data.combined)

# Soybean nodules collected from Becker, MN:
load("senescentPlusJulyNodPHB_Becker.RDATA")
# There were issues with my original set of PHB standards that led me to overestimate PHB/cell. I re-stained and reran samples later 
# and figured out that my standards done a month later (phbEstfromNovStandards) were more accurate.
becker <- julynods
rm(julynods,julynodsremeasured,novStandards,octStandards,senescentnods)
# make sure PHB estimates are zero-bound:
becker$phbEst_zerobound <- becker$phbEstfromNovStandards
becker[becker$phbEstfromNovStandards<0,"phbEst_zerobound"] <- 0


# Waseca nodules collected in (2014), from the LTARN location, measured in trap plants in single inoculation:
load("SS3regatedEditedPHBdata.RData")
was14 <- ss3phbdata
rm(ss3phbdata,ss3standards)
#make sure PHB estimates are zero-bound:
was14$phbEst_zerobound <- was14$phbEst
was14[was14$phbEst<0, "phbEst_zerobound"] <- 0
# check for mistakes--duplicated plantID--and remove redundant measurements (nodules measured twice)
zz <- table(was14$plantID,was14$strain)
zz1 <- zz>0
strainsperplant <- apply(zz1,1,sum)
strainsperplant[strainsperplant>1]
# plant #285 is used twice
was14[was14$plantID==285,]
# I don't have any other data on these plant IDs, so I'll just change the names to distinguish them
was14[was14$plantID==285 & was14$strain=="E2","plantID"] <- "285a"
was14[was14$plantID==285 & was14$strain=="H6","plantID"] <- "285b"
# Now, look for double-measured nodules:
dupnods <- was14$nodID[duplicated(was14$nodID)]
dupdata <- was14[was14$nodID %in% dupnods,]
dupdataw <- reshape(dupdata,idvar="nodID",timevar="dilution", direction="wide")
plot(phbEst_zerobound.100 ~ phbEst_zerobound.10, data=dupdataw)
abline(a=0,b=1)
# there were a bunch of samples I measured at two dilutions. It looks like the measurements are
# pretty similar (based on the graph), so I'll just take the mean.
was14.raw <- was14
was14 <- ddply(was14, .(nodID), summarize, strain=unique(strain),plantID=unique(plantID),
               phbEst_zerobound=mean(phbEst_zerobound,na.rm=TRUE))
dim(was14.raw)
dim(was14) # double checking that it had the expected number of dimensions.
dim(dupdataw)
# I know I need to exclude strain B6 because it was contaminated with dsRed-marked USDA 110
was14 <- was14[was14$strain !="B6",]

#take within-strain mean (since I want frequency by strain):
was14bystrain <- ddply(was14, .(strain),summarize, meanPHB = mean(phbEst_zerobound, na.rm=TRUE),sdPHB= sd(phbEst_zerobound,na.rm=TRUE))
# also take within-plant mean:
was14byplant <- ddply(was14, .(plantID),summarize,strain=unique(strain) ,meanPHB = mean(phbEst_zerobound, na.rm=TRUE),sdPHB= sd(phbEst_zerobound,na.rm=TRUE))


#singly inoculated plants with two strains from waseca (A3 and B5)
a3b5 <- read.csv("flowDataA3B5nod_bdpy.csv",stringsAsFactors = F) 
a3b5standard <- read.csv("PHBstandards_A3B5.csv",stringsAsFactors = F)
a3b5standard <- a3b5standard[a3b5standard$file !="5Mar15.001",]
plot(FL1mean ~ nodAgeDays, data=a3b5)
cal <- lm(FLmean ~ PHBperCellpg, data=a3b5standard[a3b5standard$stain=="bodipy",])
a3b5$phbEst <- with(a3b5,(FL1mean - coef(cal)[1])/coef(cal)[2])
a3b5$phbEst_zeroscaled <- a3b5$phbEst
a3b5[a3b5$phbEst<0,"phbEst_zeroscaled"] <- 0
a3b5[grep("A3",a3b5$plantID),"strain"] <- "A3"
a3b5[grep("B5",a3b5$plantID),"strain"] <- "B5"
a3b5 <- a3b5[a3b5$plantID !="BodipyControl",]
# any repeated measurements?
a3b5$nodID <- with(a3b5, paste(plantID,nodAgeDays,nodLet,sep="."))
repeats <- a3b5$nodID[duplicated(a3b5$nodID)]
repdata <- a3b5[a3b5$nodID %in% repeats,]
repdata <- repdata[order(repdata$nodID),]
# I already went through my flow cytometry data and marked which files to delete 
#(based on low cell count/samples being too dilute). There was a malfunction in the flow cytometer that I didn't learn about until later. 
a3b5 <- a3b5[is.na(a3b5$delete),]
repeats <- a3b5$nodID[duplicated(a3b5$nodID)]
length(repeats) # that takes care of all of them.

dma3b5 <- read.csv("20150206biomassA3B5.csv",stringsAsFactors = F,skip=1)
a3b5 <- merge(a3b5,dma3b5[,names(dma3b5)[names(dma3b5)!="strain"]],by="plantID",all.x=TRUE)
table(a3b5$plantID,a3b5$nodulesRmvd) # based on my lab notebook, nodulesRmvd is the number I wrote down for 
# nodules that were measured for PHB. For the most part, it's reasonably accurate (off by one or two nodules for a couple plants)
# noduleNum is the number of nodules that were weighed.

# Waseca nodules from trap plants in the greenhouse (soil collected in 2016)
load("phbplusplantdata_final.RDATA")
# use the LTARN experiment only (I don't want to use too much of this, since it's going in another paper)
was17 <- phbdata
was17$rotation2 <- factor(was17$rotation2)
rm(mixedphbdata, pd, phbdata, standards, uninoculatedplants)
was17$rotplot <- paste(was17$rotation2,was17$plotrep,sep=".")


##Split root plants inoculated with a subset of trap-plant nodule isolates:
##compare phenotypic variation when I know it's the same strain vs. when I don't know:
load("editedPHBData_WREsplitroot1.RDATA")
load("noduleMassWithH2assays.RDATA")
splitroot <- phbdata
rm(phbdata,standards)
# get rid of NA's and repeat measurements
# This data has columns for PHB fluorescence using arithmetic mean (FL1mean) and geometric mean (FL1geomean). I double checked the data-editing
# script and confirmed that the PHB estimates use geometric mean. 
splitroot <- splitroot[!is.na(splitroot$FL1geomean),]
splitroot$nodID <- paste(splitroot$plate,splitroot$column,splitroot$row,sep=".")
repmeasures <- unique(splitroot$nodID[duplicated(splitroot$nod)])
repdata <- splitroot[splitroot$nodID %in% repmeasures,]
repdata #there's only one nodule with a repeated measurement. One of them has 'discard' in the notes column.
splitroot <- splitroot[splitroot$file != "wresplitroot1.234",]
splitroot$focal <- "focal"
splitroot[splitroot$strain=="A3","focal"] <- "ref"
splitroot_byplant <- ddply(splitroot, .(plantID,focal), summarize,meanPHB= mean(phbEst, na.rm=TRUE), sdPHB= sd(phbEst, na.rm=TRUE),strain=unique(strain))
splitroot_byplantw <- reshape(splitroot_byplant, idvar="plantID",timevar = "focal",direction="wide")
splitroot_bystrain <- ddply(splitroot, .(strain),summarize, meanPHB = mean(phbEst, na.rm=TRUE), sdPHB = sd(phbEst, na.rm=TRUE),iqr=IQR(phbEst, na.rm=TRUE))
splitroot_bystrainw <- ddply(splitroot_byplantw,.(strain.focal), summarize, meanPHB.focal = mean(meanPHB.focal,na.rm=TRUE), sdPHB.focal = sd(meanPHB.focal, na.rm=TRUE),
                             meanPHB.ref = mean(meanPHB.ref, na.rm=TRUE), sdPHB.ref = sd(meanPHB.ref, na.rm=TRUE))

##show distribution of trap plant PHB and mark isolates from split root plants on graph:
##Need to match up isolate codes with nodule locations in trap plant data:
isolatecodes <- read.csv("isolateCodesForFreezerwithfilenames_WRE.csv",stringsAsFactors = F)

###Plots for Fig. 1, Phenotypic variation among nodules with sympatric rhizobia.
#### Fig. 1A: distribution in field-collected nodules:
##Put field-collected nodules on the same plot (Becker and Chamaecrista)
xrange <- c(0,2.7)
x <- becker$phbEst_zerobound
dens <- density(x,from=0)
x1 <- cfas$phbEst
dens1 <- density(x1,from=0)

#tiff(file="Fig1AVp_field.tiff", width=500,height=450,res=100)
plot(dens,type="l",lwd=2,main="field-collected nodules",cex.lab=1.4,cex.axis=1.3,xlab="mean PHB/cell (pg)",xlim=xrange,
     ylim=c(0,max(dens$y,dens1$y)))
points(dens1, type="l", lwd=2,col=grey(0.6))
lab1 <- paste("soybean","(n=",dens$n,")",sep="")
lab2 <- paste("Chamaecrista","(n=",dens1$n,")",sep="")
rug(x,ticksize=0.07)
rug(x1,pos=0.1,ticksize=0.07,col=grey(0.6))
text(1,2.5,labels=lab1,col="black",adj=c(0,0),cex=1.4)
text(1,2.2,labels=lab2,col=grey(0.6),adj=c(0,0),cex=1.4)
#dev.off()

####Fig. 1B. trap plant nodules:
#tiff(file="Fig1BVp_trapplants.tiff",width=500,height=450,res=100)
x2 <- was17$phbEst_zerobound
dens2 <- density(x2,from=0)
x3 <- isolatecodes$phbEst
plot(dens2,xlim=xrange,lwd=2,main="trap plant nodules",xlab="mean PHB/cell (pg)",cex.lab=1.4, cex.axis=1.3)
rug(x2, ticksize=0.05)
rug(x3,pos=0.1,col=grey(0),lwd=1.5)
text(1,1.5,labels=paste("soybean (n=",dens2$n,")",sep=""),cex=1.4,adj=c(0,0))
#dev.off()
##Plot strain means with 95% confidence interval for reference strain (holding genetics constant)
# and for focal strain (varying genetics).

###Fig. 1C: split-root plants:
splitroot_byplant$focal <- factor(splitroot_byplant$focal, levels=c("ref","focal"))
# 95% confidence interval of PHB measurements across all individual nodules (captures within/between plant variation as well as strain variation)
totCIbynod <- tapply(splitroot$phbEst,splitroot$focal,function(x){quantile(x,na.rm=TRUE,probs=c(0.025,0.975))})

#tiff(file="Fig1CsplitrootVp.tiff",width=500,height=500,res=100)
plot(meanPHB.focal ~ meanPHB.ref, data=splitroot_bystrainw, ylim=c(0,0.6), xlim=c(0,0.6),
     xlab="mean PHB/cell of shared reference strain (pg)",ylab="mean PHB/cell of focal isolate (pg)", cex.lab=1.4,cex.axis=1.3,type="n", 
     main="split-root plants")
rect(xleft = totCIbynod$ref[1],ybottom = totCIbynod$focal[1],xright=totCIbynod$ref[2],ytop=totCIbynod$focal[2],col = grey(0.5,alpha=0.5))
points(meanPHB.focal ~ meanPHB.ref, data=splitroot_bystrainw,pch=1,cex=1.5,lwd=2)
text(0.45,0.02,labels="n=21 nodule isolates",cex=1.4)
#dev.off()

# how does the breadth of varition (width of the confidence interval)
# compare between focal and reference strains?
ciwidth_foc <- totCIbynod$focal[2]-totCIbynod$focal[1]
ciwidth_ref <- totCIbynod$ref[2]-totCIbynod$ref[1]
ciwidth_foc/ciwidth_ref

#############Paritioning environmental vs. genetic variation:
#####Do differences in plant size contribute to variation between plants?
was17byplant <- ddply(was17, .(plantID),summarize, meanPHB=mean(phbEst_zerobound,na.rm=TRUE),
                      shootMassG=mean(shootMassG,na.rm=TRUE),nodMassG=mean(nodMassG,na.rm=TRUE),
                      plot= unique(plot),rotation=unique(rotation2),nodCount=mean(nodCount,na.rm=TRUE))
                      
chamabyplant <- ddply(cfas, .(plantID), summarize, meanPHB= mean(phbEst, na.rm=TRUE),
                      shootMassG = mean(shootDWg, na.rm=TRUE))                                                    

a3b5byplant <- ddply(a3b5,.(plantID), summarize, meanPHB=mean(phbEst,na.rm=TRUE),shootMassG= mean(shootMassg, na.rm=TRUE),strain=unique(strain))
#put shoot vs. mean PHB across nodules in graphs:
#tiff("FigC3chamaecrista_shootVPHB.tiff", width=550, height=500, res=100)
par(mar=c(5,5,2,2)+0.1)
plot(phbEst ~ shootDWg, data=cfas,lwd=2,xlab="shoot dry mass (g)"
     ,ylab="PHB/cell (pg)",cex.lab=1.4,cex.axis=1.3,type="n")
# put a line connecting nodules from the same plant
multiplants <- names(table(cfas$plantID)[table(cfas$plantID)>1])
zz <- matrix(NA,ncol=length(multiplants),nrow=3,dimnames = list(c("shootDWg","minPHB","maxPHB"),c(multiplants)))
zz[1,] <- with( cfas[cfas$plantID %in% multiplants,] , tapply(shootDWg, plantID, function(x){mean(x,na.rm=T)}))
zz[2,] <- with( cfas[cfas$plantID %in% multiplants,] , tapply(phbEst, plantID, function(x){min(x,na.rm=T)}))
zz[3,] <- with( cfas[cfas$plantID %in% multiplants,] , tapply(phbEst, plantID, function(x){max(x,na.rm=T)}))
rect(xleft=zz[1,]-0.008,ybottom=zz[2,]-0.02,xright = zz[1,]+0.008,ytop = zz[3,]+0.02,lwd=2)
points(phbEst ~ shootDWg, data=cfas[cfas$growthnamesStage=="veg",],pch=21, lwd=2,cex=1.5)
points(phbEst ~ shootDWg, data=cfas[cfas$growthStage=="bud",],pch=21, bg="grey",lwd=2,cex=1.5)
points(phbEst ~ shootDWg, data=cfas[cfas$growthStage=="flw",],pch=25,lwd=2,cex=1.5)
points(phbEst ~ shootDWg, data=cfas[cfas$growthStage=="earlypod",],pch=25, bg="grey",lwd=2,cex=1.5)
legend("topright",legend = c("vegetative","buds present","flowering","fruiting"),
       pch=c(21,21,25,25),pt.bg=rep(c("white","grey"),2),pt.lwd=2,pt.cex=1.5,box.lty=1,cex=1.4)
par(mar=c(5,4,4,2)+0.1)
#dev.off()


#tiff("FigC3WRE_shootVPHB.tiff", width=500, height=500, res=100)
plot(meanPHB ~ shootMassG, data=was17byplant,lwd=2,xlab="shoot dry mass (g)",ylab="PHB/cell (pg)",cex.lab=1.4,cex.axis=1.3,cex=1.5)
#dev.off()

##This one doesn't have enough plant replicates to be meaningful
#tiff("FigC3A3B5_shootVPHB.tiff", width=500, height=500, res=100)
plot(meanPHB ~ shootMassG, data=a3b5byplant,lwd=2,xlab="shoot dry mass (g)",ylab="PHB/cell (pg)",cex.lab=1.4,cex.axis=1.3,type="n")
points(meanPHB ~ shootMassG, data=a3b5byplant[a3b5byplant$strain=="A3",],lwd=2,cex=1.5,pch=19)
points(meanPHB ~ shootMassG, data=a3b5byplant[a3b5byplant$strain=="B5",],lwd=2,cex=1.5,pch=21)
legend(0.3,0.3,pch=c(19,21),legend=c("high-PHB isolate", "low-PHB isolate"),
       pt.cex=1.5,pt.lwd=2,cex=1.3)
#dev.off()

ggplot(was17byplant, aes(nodMassG, shootMassG)) + geom_point(aes(size=meanPHB))
ggplot(was17byplant, aes(nodCount, shootMassG)) + geom_point(aes(size=meanPHB))

##Within-plant variation due to nodule age:
#tiff("FigC2PHBvsNodAge.tiff", width=500,height=500,res=100)
plot(jitter(a3b5$nodAgeDays,1),a3b5$phbEst_zeroscaled,xlab="nodule age (days)", ylab="PHB/cell(pg)",cex.lab=1.4,cex.axis=1.3,type="n")
points(phbEst_zeroscaled ~ jitter(nodAgeDays,1), data=a3b5[a3b5$strain=="B5",],pch=19,cex=1.5,col=grey(0.5))
points(phbEst_zeroscaled ~ jitter(nodAgeDays,1), data=a3b5[a3b5$strain=="A3",],pch=1,cex=1.5)
#dev.off()

# any statistical assocaition between PHB and nodule age? (use B5 nodules only, since they span a wider range)
b5 <- a3b5[a3b5$strain=="B5",]
nodagevphb <- lm(phbEst_zeroscaled ~ nodAgeDays, data=b5)
plot(phbEst_zeroscaled ~ nodAgeDays, data=b5)
summary(nodagevphb)
# no significant association between PHB and nodule age--for nodules differing in age by about 2 weeks. 

##variance partitioning: two strains only, single inoculation:
va3b5 <- lmer(phbEst_zeroscaled ~ 1+ (1|strain) + (1|plantID), data=a3b5)
varparts_a3b5 <- as.data.frame(VarCorr(va3b5))

varparts_a3b5['betweenStrain','vcov']/sum(varparts_a3b5$vcov) # 93% of variation is between strain
icc(va3b5) # this does the same thing (just wanted to double check I interpreted the lmer output correctly)

##Variance partitioning with lmer (random effects)
# proportion of variance in nodule PHB explained by within-plant vs. between-plant differences:
rbeck <- lmer(phbEst_zerobound ~ 1 + (1|plantID), data=becker) # field-collected soybean nodules from Becker
summary(rbeck)
pvbeck <- data.frame(VarCorr(rbeck))
pvbeck[pvbeck$grp=="Residual","vcov"]/sum(pvbeck$vcov)
icc(rbeck)
1-icc(rbeck)

# trap plant soybean nodules (Waseca rotation experiment):
rwas17 <- lmer(phbEst_zerobound ~ 1 + (1|plantID) , data=was17)
pvwas17 <- data.frame(VarCorr(rwas17))
icc(rwas17)
pvwas17[pvwas17$grp=="Residual","vcov"]/sum(pvwas17$vcov)


####For the split root plants: Is a significant proportion of PHB variance explained by strain?
# hypothesis test: generate a null distribution of % variance explained by strain by randomly reassigning plant replicates to strains
# (permutation test) and calculating intraclass Correlation for strain (proportion of variance explained by strain). 
# Variation among strains is signficant if the proportion of variance explained by strain is larger than 95% of the null values (1-sided), then 
# make a dataset that excludes that common reference strain, but still has focal strain:
splitroot_focal <- splitroot[splitroot$strain !="A3",]
splitroot_focal$strain <- factor(splitroot_focal$strain)
focalplants <- splitroot_byplantw[!is.na(splitroot_byplantw$strain.focal),"plantID"] # unique plant IDs, excluding plants that
# had the references strain only (growth on only one side of the pouch)
# get vector of strains corresponding to each focal plant ()
focalstrains <-splitroot_byplantw[!is.na(splitroot_byplantw$strain.focal),"strain.focal"] 
splitroot_byplantw_focal <- splitroot_byplantw[splitroot_byplantw$plantID %in% focalplants,]

vobs <- lmer(phbEst ~ 1 + (1|strain), data=splitroot_focal)
pvstrain_obs <- as.numeric(icc(vobs))

nboot <- 1e4
# make a matrix of randomly permuted plantIDs
nullpla <- c()
set.seed(171103)
for(i in 1:nboot){
  nullpla <- cbind(nullpla,sample(focalplants,replace=F))
}

nullpvstrain <- apply(nullpla,2,function(x){
  pla_strain <- data.frame(plantID= x, strain= focalstrains)
  null.data <- merge(pla_strain,subset(splitroot_focal,select=-c(strain)))
  v.null <- lmer(phbEst ~ 1+  (1|strain), data= null.data)
  return(as.numeric(icc(v.null)))
})

hist(nullpvstrain,main="proportion of variance explained by strain (null)",xlim=c(0,0.7))  
q95null <- quantile(nullpvstrain,prob=c(0.95))
abline(v=pvstrain_obs)
abline(v=q95null,col="red")
# p< 0.0001 (none of the 10000 permutations produced higher proportion of variance explained by strain)
# instead of using the 95% confidence interval, it would make more sense here 
# to see if any of the permutations produced heritability estimates that were as 
#high or higher than observed--but it should be including the 95% confidence interval


##Use bootstrapping to generate a 95% confidence interval for heritability 
#(randomly resample plants with replacement)
bs.pla <- c()
nboot <- 1e4
set.seed(1711003.1)
for(i in 1:nboot){
  bs.pla <- cbind(bs.pla, sample(focalplants,replace=T))
}

bspvstrain <- apply(bs.pla,2,function(x){
  pla_strain <- data.frame(plantID= x)
  bs.data <- merge(pla_strain,splitroot_focal,by="plantID")
  v.bs <- lmer(phbEst ~ 1+  (1|strain), data= bs.data)
  return(as.numeric(icc(v.bs)))
})

bs.ci95 <- quantile(bspvstrain,prob=c(0.025,0.975))
bs.ci95
pvstrain_obs -bs.ci95


sum(nullpvstrain >= pvstrain_obs) # none of the null estimates overlapped
max(nullpvstrain)

#what proportion of null values are within the 95% confidence interval?
sum(nullpvstrain>=(bs.ci95[1]))/nboot # still <0.0001

############
# Another important question for split-root plants is whether PHB in the reference strain varies more when paired with the same 
# field isolate vs. different field isolates. 
splitroot_ref <- splitroot[splitroot$focal=="ref",]
# need to make a column for the strain that is on the other side of the plant (and maybe its mean PHB too).
plantz <- splitroot_byplantw$plantID
for(i in 1:length(plantz)){
  splitroot_ref[splitroot_ref$plantID==plantz[i], "strain.focal"] <- splitroot_byplantw$strain.focal[i]
  splitroot_ref[splitroot_ref$plantID==plantz[i], "meanPHB.focal"] <- splitroot_byplantw$meanPHB.focal[i]
}
# double check that I assigned the correct strains to the correct plants:
unique(splitroot_ref[splitroot_ref$strain.focal=="A1","plantID"])
splitroot_byplantw[splitroot_byplantw$strain.focal=="A1","plantID"]
# any missing?
nofoc <- unique(splitroot_ref[is.na(splitroot_ref$strain.focal),"plantID"])
splitroot_byplantw[splitroot_byplantw$plantID %in% nofoc,] # three plants had root growth on only the A3 side of the pouch. 
# (singly inoculated)
splitroot_ref[splitroot_ref$plantID %in% nofoc,"strain.focal"] <- "none"
splitroot_refbyplant <- ddply(splitroot_ref,.(plantID),summarize, strain.focal=unique(strain.focal), meanPHB.ref = mean(phbEst,na.rm=TRUE), 
                              meanPHB.foc = mean(meanPHB.focal,na.rm=TRUE))
# remove focal strains that have fewer than 3 plant replicates 
zz <- table(splitroot_refbyplant$strain.focal)
zz[zz>=3] # includes only 12 strains instead of 19, but that's OK
splitroot_ref3reps <- splitroot_ref[splitroot_ref$strain.focal %in% names(zz[zz>=3]),]

## View data on reference strain paired with different field isolates:
plot(phbEst ~ factor(strain.focal), data=splitroot_ref)
plot(phbEst ~ factor(strain.focal), data=splitroot_ref3reps)
# now, partition variance:
vref <- lmer(phbEst ~ 1 + (1|strain.focal), data =splitroot_ref3reps)
summary(vref)
pvref <- as.numeric(icc(vref))
# variance between focal strains (for the same reference) accounts for about 3% of total variance.
1-pvref

# what about variance within and among plants for reference strain--regardless of the focal strain?
vref2 <- lmer(phbEst ~ 1 + (1|plantID), data =splitroot_ref3reps)
icc(vref2)

vref2 <- lmer(phbEst ~ 1 + (1|plantID) + (1|strain.focal), data =splitroot_ref3reps)
icc(vref2)

# since PHB distribution is more-or less normal for the reference strain, I can use a t-test to compare 
# between focal strains:
ref_pairwise <- with(splitroot_ref3reps, pairwise.t.test(phbEst, factor(strain.focal), pool.sd = F, p.adjust.method = "none"))
hist(ref_pairwise$p.value) # not sure an FDR test correction works. 
ref_pairwiseBH <- with(splitroot_ref3reps, pairwise.t.test(phbEst, factor(strain.focal), pool.sd = F, p.adjust.method = "BH"))
pvals <- ref_pairwiseBH$p.value
pvals[pvals<0.05] # none of the comparisons are significantly different after multiple test correction. 

# I also have data on nodule occupancy:
nodMassw[is.na(nodMassw$nodcount.foc),"nodcount.foc"] <- 0
nodMassw[is.na(nodMassw$nodcount.ref),"nodcount.ref"] <- 0
nodMassw$nodcount.total <- nodMassw$nodcount.foc + nodMassw$nodcount.ref
nodMassw$nodOccupancy.ref <- nodMassw$nodcount.ref/nodMassw$nodcount.total
splitroot_ref3reps <- merge(splitroot_ref3reps, nodMassw[,c("plantID","nodOccupancy.ref","nodcount.total")])
splitroot_ref <- merge(splitroot_ref, nodMassw[,c("plantID","nodOccupancy.ref","nodcount.total")])
splitroot_refbyplant <- merge(splitroot_refbyplant,nodMassw[,c("plantID","nodOccupancy.ref","nodcount.total")] )

plot(splitroot_ref3reps$nodOccupancy.ref,splitroot_ref3reps$phbEst)
plot(splitroot_refbyplant$nodOccupancy.ref,splitroot_refbyplant$meanPHB.ref)
xyplot(phbEst ~ nodOccupancy.ref | factor(strain.focal), data=splitroot_ref3reps)
### Get bootstrap CI and do a hypothesis test: is the observed between-focal strain variance different than what you'd get from randomly 
# reassigning plants to focal strains?
nboot <- 1e4
# make a matrix of randomly permuted plantIDs
 plantz <- unique(splitroot_ref3reps$plantID)
 splitroot_refbyplant3reps <- splitroot_refbyplant[splitroot_refbyplant$plantID %in% plantz,]

nullpla <- c()
set.seed(180108.1)
for(i in 1:nboot){
  nullpla <- cbind(nullpla,sample(splitroot_refbyplant3reps$plantID, replace=F))
  }

nullpvref <- apply(nullpla,2,function(x){
  pla_strain <- data.frame(plantID = x, strain.focal = splitroot_refbyplant3reps$strain.focal)
  null.data <- merge(pla_strain,subset(splitroot_ref3reps,select=-c(strain.focal)))
  v.null <- lmer(phbEst ~ 1+  (1|strain.focal), data= null.data)
  return(as.numeric(icc(v.null)))
})

# what proportion of null icc values (proportion of variance in the reference strain due to variation between focal strain) 
#are greater than or equal to the observed icc?
sum(nullpvref >= pvref)/length(nullpvref) # over 80%
hist(nullpvref)
abline(v=pvref, col="red")
abline(v=quantile(nullpvref,0.95),col="blue")

#############

##Do a similar variance analysis for singly inoculated plants (waseca 14)
plot(phbEst_zerobound ~ factor(strain), data=was14)

vwas14 <- lmer(phbEst_zerobound ~ 1 + (1|strain), data=was14)
summary(vwas14)
pvstrainwas14_obs <- as.numeric(icc(vwas14))
pvstrainwas14_obs
# generate a bootstrap 95% confidence interval (randomly resample plantIDs with replacement)
focalplants_was14 <- was14byplant$plantID
bs.pla <- c()
nboot <- 1e4
set.seed(1711004)

for(i in 1:nboot){
  bs.pla <- cbind(bs.pla, sample(focalplants_was14,replace=T))
}

bspvstrain_was14 <- apply(bs.pla,2,function(x){
  pla_strain <- data.frame(plantID= x)
  bs.data <- merge(pla_strain,was14,by="plantID")
  v.bs <- lmer(phbEst_zerobound ~ 1+  (1|strain), data= bs.data)
  return(as.numeric(icc(v.bs)))
})
hist(bspvstrain_was14)
abline(v=pvstrainwas14_obs)
bs.ci95_was14 <- quantile(bspvstrain_was14,prob=c(0.025,0.975))
# number to present in text (estimated broad-sense heritability in singly inoculated plants,28 strains)
pvstrainwas14_obs
pvstrainwas14_obs - bs.ci95_was14

#hypothesis test for single inoculation (10,000 permutations)
nboot <- 1e4
# make a matrix of randomly permuted plantIDs
nullpla <- c()
set.seed(171104.1)
focalstrains_was14 <- was14byplant$strain
for(i in 1:nboot){
  nullpla <- cbind(nullpla,sample(focalplants_was14,replace=F))
}

nullpvstrain_was14 <- apply(nullpla,2,function(x){
  pla_strain <- data.frame(plantID= x, strain= focalstrains_was14)
  null.data <- merge(pla_strain,subset(was14,select=-c(strain)),by="plantID")
  v.null <- lmer(phbEst_zerobound ~ 1+  (1|strain), data= null.data)
  return(as.numeric(icc(v.null)))
})

hist(nullpvstrain_was14,main="proportion of variance explained by strain (null)",xlim=c(0,0.9))  
q95null <- quantile(nullpvstrain_was14,prob=c(0.95))
abline(v=pvstrainwas14_obs)
abline(v=q95null,col="red")
sum(nullpvstrain_was14>=pvstrainwas14_obs)/nboot
# none of the null estimates overlap with the observed--but what about
# the bottom of the 95% confidence interval for the observed?
# same story. The max of the null is 0.71, the observed is between 0.78 and 0.92
sum(nullpvstrain_was14 >=bs.ci95_was14[1])
max(nullpvstrain_was14)

# same story, p<0.0001


##Appendix D: 
# Show phenotypic variability in PHB/cell (nodule average) within strains in single inoculation.
was14bystrain <- was14bystrain[order(was14bystrain$meanPHB),]
# order factor levels by mean PHB (order by strain)
wstrainsmeans <- sort(tapply(was14$phbEst_zerobound, was14$strain,function(x){mean(x,na.rm=TRUE)}))
wstrainorder <- names(wstrainsmeans)
was14$strain_byphb <- factor(was14$strain, levels=wstrainorder)

#tiff(file="FigD1_Waseca_withinstrainVp.tiff", width=500, height=500, res=100)
boxplot(phbEst_zerobound~ strain, data=was14,horizontal=TRUE,yaxt="n",cex.axis=1.3,
        ylab="strain",cex.lab=1.4, xlab="mean PHB/cell (pg)")
points(wstrainsmeans, c(1:27),pch=4,col=grey(0.5), lwd=2)
axis(2,at=1:27, labels=FALSE)
# mark USDA 110 with an arrow.
arrows(0.3, which(levels(was14$strain_byphb)=="110") , 0.2, which(levels(was14$strain_byphb)=="110"))
#dev.off()

plotdist <- function(x,title="title",bw=1,nodorstrain="nodules"){
  dens <- density(x, from=0, adjust=bw)
  plot(dens,type="l",lwd=2,main=title,cex.lab=1.4,cex.axis=1.3,xlab="mean PHB/cell (pg)")
  legend("topright",legend=c( paste("n=",dens$n," ",nodorstrain,sep=""), paste("bandwidth=",round(dens$bw,5),sep="")),bty="n",cex=1.4)
  rug(x,ticksize=0.05)
}

#tiff(file="Fig1_dist_was14bystrain.tiff",width=500,height=500,res=100)
plotdist(was14bystrain$meanPHB,title="Waseca strains in growth pouches",nodorstrain="strains")
usda110phb <- was14bystrain[was14bystrain$strain=="110","meanPHB"]
arrows(usda110phb, 4, usda110phb, 1)
#dev.off()

#Becker soybean (field)
bpla <- sort(tapply(becker$phbEst_zerobound, becker$plantID,mean))
bids <- names(bpla)
becker$plantID_byphb <- factor(becker$plantID, levels=bids)

#tiff(file="becker_withinplantVp.tiff", width=500, height=500,res=100)
plot(phbEst_zerobound ~ plantID_byphb, data=becker,horizontal=TRUE,yaxt="n",
     xlab="host plant",ylab="mean PHB/cell (pg)",cex.lab=1.4, cex.axis=1.3)
axis(2, labels=FALSE, at=1:length(bpla))
par(yaxt="s")
#dev.off()

# LTARN trap plants
ltarnpla <- sort(tapply(was17$phbEst_zerobound, was17$plantID, function(x){mean(x,na.rm=TRUE)}))
was17$plantID_byphb <- factor(was17$plantID, levels=names(ltarnpla))
#tiff("ltarn_witihinplantVp.tiff", width=500, height=500, res=100)
plot(phbEst_zerobound ~ plantID_byphb, data=was17,horizontal=TRUE,yaxt="n", xlab="host plant",ylab="mean PHB/cell (pg)",
     cex.lab=1.4, cex.axis=1.3)
axis(2,labels=FALSE, at=1:length(ltarnpla))
par(yaxt="s")
#dev.off()

##Taxonomic differences in PHB accumulation: based on tolerance to BG:
####Chamaecrista fasciculata: calculate bootstrapped means and standard error for comparing PHB of nodules with 
# BG-sensitive vs. tolerant rhizobia
bs.sensitive <- c()
bs.tolerant <- c()
nboot <- 1e4
set.seed(171104.3)
for(i in 1:nboot){
  bs.sensitive <- cbind(bs.sensitive,sample(cfas[cfas$growsWithBG==FALSE,"phbEst"],replace=TRUE))
  bs.tolerant <- cbind(bs.tolerant, sample(cfas[cfas$growsWithBG==TRUE,"phbEst"],replace=TRUE))
}
sensitive.mean <- mean(apply(bs.sensitive,2,mean))
sensitive.se <- sd(apply(bs.sensitive,2,mean))
tolerant.mean <- mean(apply(bs.tolerant,2,mean))
tolerant.se <- sd(apply(bs.tolerant,2,mean))
print(
  rbind(paste("BGsensitivemeanPHB =",round(sensitive.mean,3),"+/-",round(sensitive.se,3),sep=" "),
        paste("BGtolerantmeanPHB =",round(tolerant.mean,3),"+/-",round(tolerant.se,3),sep=" "))
)


cfas$BGtolerance <- "tolerant"
cfas[cfas$growsWithBG==FALSE,"BGtolerance"] <- "sensitive"

#tiff("BGvsPHB.tiff", width=500, height=500,res=100)
plot(cfas$phbEst ~ factor(cfas$BGtolerance),xlab = "BG tolerance",
     ylab="nodule mean PHB/cell (pg)", cex.lab=1.4 ,cex.axis=1.3)
text(1,0.2,labels="*", cex=2)
#dev.off()

t.test(log(phbEst) ~ growsWithBG, data=cfas)
plot(log(cfas$phbEst) ~ factor(cfas$BGtolerance),xlab = "BG tolerance",
     ylab="nodule mean PHB/cell (pg)", cex.lab=1.4 ,cex.axis=1.3)
# a problem with this approach is that 9 of my BG-sensitive nodules came from the same
# plant. So this could be pseudoreplication. Do I still get a difference if I collapse nodules
# from the same plant (with the same tolerance) into a single mean?
cfas1 <- ddply(cfas, .(plantID, BGtolerance), summarize,meanPHB = mean(phbEst, na.rm=TRUE))
table(cfas1$plantID)
table(cfas1$BGtolerance)
#that's only 5 sensitive and 9 tolerant. 

plot(meanPHB ~ factor(BGtolerance), data=cfas1 )
plot(log(meanPHB) ~ factor(BGtolerance), data=cfas1 )
t.test(log(meanPHB) ~ factor(BGtolerance), data=cfas1)

# distribution doesn't look good either way. Better to use a bootstrap test:
set.seed(42)
nboot <- 1e5
mean.diffs <- numeric(nboot)
for(i in 1:nboot){
  bstrapS = sample(cfas1[cfas1$BGtolerance=="sensitive","meanPHB"],replace=T)
  bstrapT = sample(cfas1[cfas1$BGtolerance=="tolerant","meanPHB"],replace=T)
  mean.diffs[i] <- mean(bstrapT )- mean(bstrapS)
}
perm1e4 <- c()
for(i in 1:1e4){
  perm1e4 <- cbind(perm1e4, sample(as.character(cfas1$BGtolerance),replace=F))
}
mean.diffs <- apply(perm1e4,2 ,function(x){
  means.perm = tapply(cfas1$meanPHB, factor(x), mean)
  perm.means.diff = means.perm[2] - means.perm[1]
  return(perm.means.diff)
})
obs.means <- tapply(cfas1$meanPHB, cfas1$BGtolerance,mean)
obs.diff <- obs.means[2] - obs.means[1]
plot(density(mean.diffs))
abline(v=c(obs.diff,-obs.diff))
2* sum(mean.diffs < -obs.diff | 
      mean.diffs > obs.diff)/1e4
 #p=0.0122. That seems more reasonable. 

#lumping together nodules from the same plant doesn't change results--means are still significantly different. 

#####################################
# since I have multiple factors, 
###real quick: just look into any differences between inoculation treatments at becker:
becker[grep("A",becker$plantID),"inoc"] <- "A3"
becker[grep("B",becker$plantID),"inoc"] <- "B5"
becker[grep("110",becker$plantID),"inoc"] <- "110"
becker[grep("U",becker$plantID),"inoc"] <- "none"
with(becker, table(inoc, plantID))
becker$inoc <- factor(becker$inoc, levels=c("none","B5","110","A3"))
plot(becker$phbEst_zerobound ~ becker$inoc)
plot((log(becker$phbEst_zerobound+0.001)) ~ becker$inoc)
hist(becker$phbEst_zerobound)
hist(log(becker$phbEst_zerobound+0.001))

becker[becker$plantID %in% c("1101","A1","B1","U1"),"plantRep"] <- "pla1"
becker[becker$plantID %in% c("1102","A2","B2","U2"),"plantRep"] <- "pla2"
becker[becker$plantID %in% c("1103","A3","B3","U3"),"plantRep"] <- "pla3"
becker[becker$plantID %in% c("1104","A4","B4","U4"),"plantRep"] <- "pla4"
becker$plantRep <- table(becker$plantRep)
table(becker$plantID)
lm.inoc <- lm(phbEst_zerobound ~ inoc + inoc:plantRep, data=becker)
lm.inoc2 <- lm(phbEst_zerobound ~ inoc            , data=becker)
lm.inocnull <- lm(phbEst_zerobound ~ 1           , data=becker)
anova(lm.inoc2, lm.inoc)
anova(lm.inocnull, lm.inoc2)
aov(phbEst_zerobound ~ inoc + inoc:plantRep, data=becker)
lm.inoc <- lmer(phbEst_zerobound ~  (1|plot) + (1|plantID), data=becker)
summary(lm.inoc)
plantv <- 1.309e-02
plotv <- 3.478e-18
residv <- 3.799e-02
residv/(plantv+plotv+residv)
plantv/(plantv+plotv+residv)
plotv/(plantv+plotv+residv)
