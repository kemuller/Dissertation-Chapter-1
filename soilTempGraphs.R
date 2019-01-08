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

