###############'
#' How long could an individual rhizobial cell in the soil survive on endogenous PHB?


##################

soilTemp <- read.csv("averagemonthlysoiltemp_waseca.csv")
date <- paste(soilTemp$year,soilTemp$month,"15", sep="-")
soilTemp$date <- as.Date(date, format=c("%Y-%m-%d"))

##Average monthly soil temperature data from Waseca from 2010-2014 (daily data is available, but would take a long time to enter by hand)
##Useful matrix for converting monthly data to daily (ignore leap day)
month2day <- matrix(NA,nrow=4,ncol=12,byrow=TRUE, dimnames=list(c("imonth","Jd1ofm","dayspermonth","Jmidm"),
                                                                c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")))
month2day[1,] <- c(1:12)
month2day[2,] <- c(1,32,60,91,121,152,182,213,244,274,305,335)
month2day[3,] <- c(31,28,31,30,31,30,31,31,30,31,30,31)
month2day[4,] <- c(16,45.5,75,105.5,136,166.5,197,228,258.5,289,319.5,350)

for(i in 1:12){
  soilTemp[soilTemp$month==c(1:12)[i],"lastdayofmonth"] <- month2day[3,][i]
  soilTemp[soilTemp$month==c(1:12)[i],"indexofday1"] <- month2day[2,][i]
}


dailysoiltemp <- soilTemp[rep(row.names(soilTemp),soilTemp$lastdayofmonth),] #expand to have one row per day
#add in column for day of month
for(i in 1:12){
  dailysoiltemp[dailysoiltemp$month==c(1:12)[i],"day"] <- c(1:month2day[3,][i])
}
#make a column for calendar date:
date <- with(dailysoiltemp, paste(year,month,day,sep="-"))
dailysoiltemp$date <- as.Date(date, format=c("%Y-%m-%d"))

##View oscilation:
plot(average4inchtempC ~ date, data=soilTemp,ylab="soil temperature (degC)",main="soil temperature measured at Waseca SROC")
points(average8inchtempC ~ date, data=soilTemp, col="red")
legend("bottomright",legend=c("4 inches","8 inches"),fill=c("black","red"),title="depth")

###Use the 8-inch temperature, since it's a bit more stable (the difference is pretty negligible in the warm season, 
# which is where it would have the biggest effects on survival estimates Get a vector of hourly data:
##To make life easier, fit a mathematical model to the temperature data:
soilTemp <- soilTemp[order(soilTemp$date),]
plot(average8inchtempC~indexofday1, data=soilTemp[soilTemp$year=="2010",],type="s") ##
abline(v=c(91,273))

##Use a 1st order Fourier equation (from p. 262 of Thornley and France)
##Maximum and minimum yearly soil temp: use average of 6 years:
Tmax <- mean(with(soilTemp,tapply(average8inchtempC,year,max)))
Tmax <- round(Tmax) #round to nearest degree
##2012 is a major outlier--exclude from mean
mintemp <- with(soilTemp,tapply(average8inchtempC,year,min))
Tmin <- mean(mintemp[row.names(mintemp)!="2012"])
Tmin <- round(Tmin)
Tmn = 0.5 * (Tmax + Tmin)
Tvr = 0.5 * (Tmax - Tmin)
de = 121 #day on which temp reaches Tmn (or de +-182)--May 1st and Nov 1st
soilTemp$sinusoidTemp <- Tmn + Tvr*sin(2*pi*(soilTemp$indexofday1-de)/365)

#format for paper
#tiff(file="FigA1tempOscillation.tiff",width=607,height=350,res=100)
plot(average8inchtempC~date, data=soilTemp,type="s",lwd=2,col=grey(0),
     ylab=c('8in soil temp(deg.C)'),cex.lab=1.2,cex.axis=1.3,xlab="") ##
points(sinusoidTemp ~ date, data=soilTemp,type="l",lwd=2,lty=1,col=grey(0.4))
#dev.off()
##Calculate calendar day as a function of time after starting day
# let's say rhizobia are released on day 271 (October 1st)
d1 = 271
x = 1:100
z = d1 + x
d = z - 365*(floor(z/365) + ceiling(z/365) - floor(z/365) - 1)
##This is more succinct:
d = z - 365*(ceiling(z/365) - 1)

#How long does a cell with given PHB last before dying?
##Initial assumptions: a cell with 2pg of structural biomass and 2 pg PHB goes into dormancy on October 1st (day 274).
# Solve iteratively--should be doable analytically, but my brain can't quite wrap around it.
#Stoichiometry
#PHB monomer = C4H6O2
#Structural biomass = C4H8O2N
gPHBtogC = 48/86.09
gSBtogC = 48/102.11
PHB = 2
PHBgC <- PHB * gPHBtogC
SBgc <- 2 * gSBtogC
d1 = 274 #day 1 = October 1st

#How much C used in 1 year:
t= 24*365 #time in hours
z = d1 + 1:t/24
d = z - 365*(ceiling(z/365) - 1)
Tc = 11 + 14*sin(2*pi*(d-121)/365)
mt = 10^(0.084*Tc - 8.17)
Cused = sum(mt * SBgc)
PHBused = Cused/gPHBtogC

##Calculate PHB used in 30 years for dormant, maintenance, and growing metabolism
ty = 30*365*24 #number of hours in 30 years (excluding leap year)
yearsdormant <- (1:ty)/(24*365)

##Dormant (aspartic acid racemization)
z = d1 + (1:ty)/24
d = z - 365*(ceiling(z/365) - 1)
Tc = 11 + 14*sin(2*pi*(d-121)/365)
mt = 10^(0.084*Tc - 8.17)
Cused = cumsum(mt * SBgc)
PHBuseddormant = Cused/gPHBtogC

##Maintenance metabolism parameters:
z = d1 + (1:ty)/24
d = z - 365*(ceiling(z/365) - 1)
Tc = 11 + 14*sin(2*pi*(d-121)/365)
mt = 10^(0.099*Tc - 5.14)
Cused = cumsum(mt * SBgc)
PHBusedmaint = Cused/gPHBtogC

##Growth metabolism parameters
z = d1 + (1:ty)/24
d = z - 365*(ceiling(z/365) - 1)
Tc = 11 + 14*sin(2*pi*(d-121)/365)
mt = 10^(0.0987*Tc - 2.161)
Cused = cumsum(mt * SBgc)
PHBusedgrowing = Cused/gPHBtogC


yearsdormant <- ty/(365*24)
plot(PHBuseddormant~yearsdormant, type="l",ylim=c(0,2))
points(PHBusedmaint~yearsdormant, type="l", col="red")
points(PHBusedgrowing~yearsdormant, type="l", col="blue")

##How long does it take to use 0.1pg PHB?
dormant.1pghours <- min(which(round(PHBuseddormant,2)==0.1))
dormant.1pghours/(24*365)
##0.1 pg PHB lasts about 30 years with dormant metabolism rates
maint.1pghours <- min(which(round(PHBusedmaint,2)==0.1))
maint.1pghours/24
##0.1pg PHB lasts about a week with maintenance metabolism rates
# 1pg PHB would last 255.6 days, or nearly a year. 2 pg PHB would last 277 days.
growth.1pghours <- min(which(round(PHBusedgrowing)==1)) #Growing, 0.1 pg PHB would be used in < 1h
#0.1pg PHB lasts < 1 hour with growing metabolism rates. 
# 1pg PHB would last 1-2 hours. 2pg PHB would last about 3 hours.

##That is 0.0018 gPHB carbon in a year, or 0.00324 g PHB per cell. 
## A cell with 0.3pg PHB could last for nearly 100 years!
## A cell with 1pg PHB would last over 300 years!
##let's look at estimated years of survival as a function of initial PHB (a year being October1st-October1st)
phbrange = seq(0,1,by=0.01)
survivalyears = phbrange/PHBused
plot(phbrange,survivalyears)
abline(h=30)
abline(h=5)

##In that case, a difference of 30 years wouldn't make much difference--unless rhizobia use up all but a fraction of nodule PHB during senescence
# (which is possible)
which(round(survivalyears)==30)
phbsurvival <- cbind(phbrange,survivalyears)

Cremaining <- numeric(t)
Cremaining[1] <- PHBgC
for(i in 2:t){
Cremaining[i] <- Cremaining[i-1] - mt[i]*SBgc  
}
plot(Cremaining) #shows slight oscillation--probably some complicated Fourier equation, but I don't really need to do it. 


