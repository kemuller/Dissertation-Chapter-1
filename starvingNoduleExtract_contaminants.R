###This is a script to test whether excluding replicates with contaminants 
# has any impact on results (or if I even have enough data)
# Something is wonky with my initial PHB data--I need to check this. 
load("starvingNoduleExtractPHBuse.RDATA")
rawcont <- read.csv("contaminationNST1.csv",stringsAsFactors = F,skip=1)

# Eventually I gave up counting contaminants and just scored them on quality.
# but I need some way to score all of them
#Are there any from samples 5-8 that don't have a "quality" score?
s5t8 <- ld[ld$sample>=5,]
sum(is.na(s5t8$quality))
s5t8$quality # nope.

# It may not be necessary to look at samples 2-4. maybe just stick to 
# initial, then 1, 5, and 15 months. 
ldplusfd <- flowdatacombined[flowdatacombined$treatment=="LD" & flowdatacombined$sample %in% c(1,5:8),]

# look at 1 month population growth--is the relationship between 
# reproduction and PHB use different if you exclude samples with contamination?
cfubystrain_ld2 <- cfubystrain_ld[!is.na(cfubystrain_ld$meanCFU),]
initialLD <- flowdatabystrain[flowdatabystrain$sample==1,"initialCFU"]
initialPHB <- flowdatabystrain[flowdatabystrain$sample==1,"initialPHB"]
strainz <- flowdatabystrain$strain

# clean samples: quality 1 or 2:
ldclean <- ldplusfd[ldplusfd$quality<3,]
ldclean <- ldclean[!is.na(ldclean$sampleID),]
with(ldclean[ldclean$sample==5,], plot(initialPHB,phbused))
with(ldplusfd[ldplusfd$sample==5,], plot(initialPHB,phbused))
# results are qualitatively similar--there are just fewer samples.

# look at population growth vs. PHB used:
with(ldclean[ldclean$sample==5,], plot(phbused,CFUchange))
with(ldplusfd[ldplusfd$sample==5,], plot(phbused,CFUchange))
with(ldclean[ldclean$sample==5,], points(phbused,CFUchange,col="red"))
ldclean$relativePop <- ldclean$meancellsperml/ldclean$initialCFU
# take average by strain using only clean samples:
ldclean_bystrain <- ddply(ldclean, .(strain, sample), meanPHBhigh = 
                            mean(phbEst_high2, na.rm=T), meanPHBlow=mean(phbEst_low, na.rm=T),
                          meanCFU = mean(meancellsperml, na.rm=TRUE),
                          CFUchange = mean(CFUchange, na.rm=TRUE), starvationdays=mean(starvationdays),
                          phbused = mean(phbused, na.rm=TRUE),
                          initialPHB=mean(initialPHB),
                          initialCFU = mean(initialCFU),
                          relativePop = mean(relativePop, na.rm=TRUE),summarize)
with(ldclean_bystrain[ldclean_bystrain$sample==5,], plot(phbused,CFUchange))                          
with(ldclean_bystrain[ldclean_bystrain$sample==5,], plot(initialPHB,phbused))

# Look at relative pop size. 
plot(relativePop ~ initialPHB,ldclean[ldclean$sample==6,])
plot(relativePop ~ initialPHB,ldclean[ldclean$sample==8,])

# Make graphs corresponding to Fig. 3A and B, using only "clean" samples.
#tiff("FigE2_3Acleanonly.tiff", width=500, height=500, res=100)
par(mar=c(5,5,3,2)+0.1)
plot(relativePop ~ initialPHB,ldclean_bystrain[ldclean_bystrain$sample==6,],pch=19, cex=2, 
     ylab= expression("relative population "~ (N[t] / N[0])),xlab="initial PHB per cell (pg)",
     cex.lab=1.4, cex.axis=1.3)
points(relativePop ~ initialPHB,ldclean_bystrain[ldclean_bystrain$sample==8,],pch=1,cex=2,lwd=2)
xx <- ldclean_bystrain[ldclean_bystrain$sample==6,"initialPHB"]
names(xx) <- ldclean_bystrain[ldclean_bystrain$sample==6,"strain"]
y_0 <- ldclean_bystrain[ldclean_bystrain$sample==6,"relativePop"]
names(y_0) <- ldclean_bystrain[ldclean_bystrain$sample==6,"strain"]
y_1 <- ldclean_bystrain[ldclean_bystrain$sample==8,"relativePop"]
names(y_1) <- ldclean_bystrain[ldclean_bystrain$sample==8,"strain"]
y_1 <- y_1[intersect(names(y_0),names(y_1))]
y_0 <-y_0[names(y_1)]
xx <- xx[names(y_1)]
arrows(xx,y_0, xx, (y_1+0.2), length=0.1, lwd=2, angle=20, col=grey(0.5))
abline(h=1, lty=2, lwd=2)
legend("topleft", pt.cex=2, pch=c(19,1), pt.lwd=2, legend=c("127 days", "450 days"), title="time in starvation culture")
par(mar=c(5,4,4,2)+0.1)
#dev.off()

#tiff("FigE2_3Bcleanonly.tiff", width = 500, height= 500, res=100)
nf <- layout(matrix(c(1,2),2,1),heights=c(1,0.5))
par(mar=c(0,5.1,3.1,2.1))
growers <- ldclean_bystrain[ldclean_bystrain$sample==5 & ldclean_bystrain$CFUchange>0,]
shrinkers <- ldclean_bystrain[ldclean_bystrain$sample==5 & ldclean_bystrain$CFUchange<0,]
xspan <- range(c(growers$phbused, shrinkers$phbused))
with(growers, plot(phbused,CFUchange,pch=19,cex=2,xlab="1 month PHB use (pg/cell)",
     ylab="",cex.lab=1.4,cex.axis=1.3,xaxt="n",yaxt="n",xlim=xspan +c(-0.01,0.01)))
axis(2,at=seq(2e6,5e6,by=0.5e6), labels=seq(2,5,by=0.5),cex.axis=1.3)
phbvrep_clean <- lm(CFUchange ~ phbused, data=growers)
abline(phbvrep_clean,lwd=2)
lm_coef <- round(coef(phbvrep_clean)/10^6,2)
r2 <- round(summary(phbvrep_clean)$r.squared,2)
eq <- substitute(paste("y = ","(",m, "*",10^6,")",
                       " x"," + ", b,"*", 10^6 ), list(m=lm_coef[2], b=lm_coef[1]))
eqr2 <- substitute(paste("R"^{2}," = ", rsquared),list(rsquared=r2))
text(0.06, 3.7e6, labels = eq, cex=1.4, pos=4)
text(0.06, 3.3e6, labels = eqr2, cex=1.4, pos=4)
par(mar=c(5.1,5.1,1,2.1))
plot(range(c(growers$phbused,shrinkers$phbused)),c(-3.5e5,-4.5e5),cex.axis=1.4, xaxt="s",ylab="",xlab="PHB used in 1st month",
     cex.lab=1.4, cex.axis=1.3, yaxt="n",type="n" , xlim=xspan +c(-0.01,0.01))
points(CFUchange ~ phbused, data=shrinkers, cex=2, pch=19)
axis(2,at=c(-3.5e5,-4e5),cex.axis=1.4,labels=c(-3.5, -4))
mtext(side=2,text= expression("net population growth "~(10^6~ "cells/mL"))
      ,xpd=NA ,line=2.5,at=-4,cex=1.4)
layout(1)
par(mar=c(5,4,4,2)+0.1,yaxt="s",xaxt="s",xpd=FALSE)
#dev.off()