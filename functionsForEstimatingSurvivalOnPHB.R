##These functions were originally written in the script survivalwithSoilTempAndNodSenescence.R. They are copied and revised here so this
# script can be run to load functions for estimating PHB used under different metabolic scenarios. 
# The analysis used to figure out functions for temperature and metabolic rates is in the script soilTempGraphs.R

##Constants:
#Stoichiometry: From Bormann et al. 2000
#PHB monomer = C4H6O2
#Structural biomass = C4H8O2N
gPHBtogC = 48/86.09
gSBtogC = 48/102.11
PHB = 2
PHBgC <- PHB * gPHBtogC
SBgc <- 0.5 * gSBtogC # I originally used 2pg of structural biomass. But based on measurements of 
# bacteroids from soybean nodules (Bergersen 1985), 0.3pg is more realistic. I'll go with 
# 0.5 to be consistent with published measurements of E. coli. 
d1 = 258 #day 1 = mid-September (around harvest)

##I found a different source (Roels 1980) that has a more reliable stoichiometry for bacterial biomass,
# but it differs from that presented in Bormann et al.:
# Roels' stoichiometry was C5H9O2.5N
C = 5
H = 9
O= 2.5
N = 1
mwSB = C*12 + H*1 + O*16 + N*14
gSBtogC2 = C*12/mwSB
##This stoichiometry is 49% C in bacterial biomass, compared to 47% C for the Bormann estimate
# I think I should stick to my original. 
rm(C,H,O,mwSB,gSBtogC2)

##Phase 1 can be high or med, phase 2 can be low or med.
calculatePHB_2phase <- function(gSB=0.5,
  phase1 = "high",phase2 = "low",t=365*30,senescenceperiod=30, d1=258, midphb1=0.1, midphb2=0.01){ #coefficient is for middle metabolic rate (proportion of high)
  z = d1 + 1:t
  d = z - 365*(ceiling(z/365) - 1)
  Tc = 11 + 14*sin(2*pi*(d-121)/365)
  gCpergPHB = 48/86.09 # stoichiometry of PHB and structural biomass (SB)
  gCpergSB = 48/102.11
  
  if(phase1=="high"){
    mphase1_hourly = 10^(0.099*Tc[1:senescenceperiod] - 5.14) #high maintenance rate
  }else{
    mphase1_hourly = midphb1*10^(0.099*Tc[1:senescenceperiod] - 5.14)} #10 % of high maintenance rate (med)
  
  if(phase2=="low"){
    mphase2_hourly = 10^(0.084*Tc[(senescenceperiod+1):t] - 8.17) # dormant metabolic rate per hour. make it per day:
  }else{
    mphase2_hourly = midphb2*10^(0.099*Tc[(senescenceperiod+1):t] - 5.14)} # 10% of high maintenance rate
  
  mphase1_daily = mphase1_hourly*24
  mphase2_daily = mphase2_hourly*24
  SBgC = gSB * gCpergSB
  PHBused_phase1 = (mphase1_daily * SBgC)/gCpergPHB
  PHBused_phase2 = (mphase2_daily * SBgC)/gCpergPHB
  PHBused <- cumsum(c(PHBused_phase1,PHBused_phase2))
  PHBused}

##Model with constant metabolic rate (plus an argument for constant temperature)
calculatePHB_1phase <- function(gSB=0.5,rate = "high",t=365*30, d1=258, midphb = 0.1,temp = -999){ #midphb coefficient is what you multiply high rate by to get middle rate
  z = d1 + 1:t
  d = z - 365*(ceiling(z/365) - 1) # function to run through the days of the year (1-365)
  Tc = 11 + 14*sin(2*pi*(d-121)/365) # temperature oscilation (in celsius)
  gCpergPHB = 48/86.09 # stoichiometry of PHB and structural biomass (SB)
  gCpergSB = 48/102.11
  
  if(temp== -999){ # this is set up so that if you don't specify temperature, it uses the oscillating temperature based on the field measurements.
    # otherwise, it will use a constant temperature (in celsius)
  mdormant_hourly = 10^(0.084*Tc[1:t] - 8.17)##Lowest rate.# dormant metabolic rate per hour. 
  mmaintenance_hourly = 10^(0.099*Tc[1:t] - 5.14) # highest rate (hourly)
  mmaintenance10_hourly = midphb* (10^(0.099*Tc[1:t] - 5.14)) # 10 % of highest rate
  }else{
    mdormant_hourly =10^(0.084*rep(temp,t) - 8.17)
    mmaintenance_hourly = 10^(0.099*rep(temp,t) - 5.14) # highest rate (hourly)
    mmaintenance10_hourly = midphb*(10^(0.099*rep(temp,t) - 5.14)) # 10 % of highest rate
    }
  mdormant_daily = mdormant_hourly*24 # multiply hourly rate by 24 to get daily rate (soil temperature is assumed constant on a daily basis)
  mmaintenance_daily = mmaintenance_hourly*24
  mmaintenance10_daily = mmaintenance10_hourly*24
  
  if(rate == "high"){
    mrate <- mmaintenance_daily
  } else{
    if (rate =="low"){
      mrate <- mdormant_daily }
    else{mrate <- mmaintenance10_daily}}
  
  gCinSB = gSB * gCpergSB # carbon in structural biomass (units are pg).
  PHBusedperday = (mrate * gCinSB)/gCpergPHB # divides carbon used per day by carbon per g PHB (to get pg PHB from pg C).
  PHBused <- cumsum(c(PHBusedperday))
  PHBused}

calculatePHB_hourly <- function(gSB=0.5,rate = "highmaintenance",t=24, midphb = 0.1,temp=26){ #midphb coefficient is what you multiply high rate by to get middle rate
  ##Calculate hourly PHB use for time periods < 1 day. d1 chooses which day of year you're at
  mgrowth_hourly = 10^(0.0987*(rep(temp,t)) - 2.161)# growth
  mmaintenance_hourly = 10^(0.099*(rep(temp,t)) - 5.14) # highest rate (hourly)
  mmaintenance10_hourly = midphb* (10^(0.099*rep(temp,t) - 5.14)) # 10 % of highest rate
  gCpergPHB = 48/86.09 # stoichiometry of PHB and structural biomass (SB)
  gCpergSB = 48/102.11
  if(rate == "highmaintenance"){
    mrate <- mmaintenance_hourly
  } else{
    if (rate =="lowmaintenance"){
      mrate <- mmaintenance10_hourly }
    else{mrate <- mgrowth_hourly}}
  gCinSB = gSB * gCpergSB
  PHBusedperhour = (mrate * SBgc)/gCpergPHB
  PHBused <- cumsum(c(PHBusedperhour))
  PHBused}

# calculatePHB_2phase <- function(phase1 = "high",phase2 = "low",t=365*30,senescenceperiod=30, d1=258, midphb1=0.1, midphb2=0.01){ #coefficient is for middle metabolic rate (proportion of high)
#   z = d1 + 1:t
#   d = z - 365*(ceiling(z/365) - 1)
#   Tc = 11 + 14*sin(2*pi*(d-121)/365)
#   
#   if(phase1=="high"){
#     mphase1_hourly = 10^(0.099*Tc[1:senescenceperiod] - 5.14) #high maintenance rate
#   }else{
#     mphase1_hourly = midphb1*10^(0.099*Tc[1:senescenceperiod] - 5.14)} #10 % of high maintenance rate (med)
#   
#   if(phase2=="low"){
#     mphase2_hourly = 10^(0.084*Tc[(senescenceperiod+1):t] - 8.17) # dormant metabolic rate per hour. make it per day:
#   }else{
#     mphase2_hourly = midphb2*10^(0.099*Tc[(senescenceperiod+1):t] - 5.14)} # 10% of high maintenance rate
#   
#   mphase1_daily = mphase1_hourly*24
#   mphase2_daily = mphase2_hourly*24
#   PHBused_phase1 = (mphase1_daily * SBgc)/gPHBtogC
#   PHBused_phase2 = (mphase2_daily * SBgc)/gPHBtogC
#   PHBused <- cumsum(c(PHBused_phase1,PHBused_phase2))
#   PHBused}
# 
# # assumes wakeup duration is 7 days.
# PHBuse_setwakeup <- function(sleeping = "low", t=365*2, d1=258,
#                              wakeupstart = c(258,91), midphb_dormant=0.001, midphb_active=1){ #coefficient is for middle metabolic rate (proportion of high)
#   z = d1 + 1:t
#   d = z - 365*(ceiling(z/365) - 1)
#   Tc = 11 + 14*sin(2*pi*(d-121)/365)
#   time <- data.frame(day=d, Tc = Tc)
#   time$state <- "dormant"
#   wakeupdays <- c()
#   for(i in 1:length(wakeupstart)){
#     wakeupdays <- c(wakeupdays,wakeupstart[i] + c(0:6))
#   }
#   time[time$day %in% wakeupdays,"state"] <- "active"
#   
#   
#   if(sleeping=="low"){
#     mdormant_daily = 24*10^(0.084*Tc[(senescenceperiod+1):t] - 8.17) # dormant metabolic rate per hour, times 24 to make it per day:
#   }else{mdormant_hourly <- midphb_dormant*10^(0.099*Tc[(senescenceperiod+1):t] - 5.14)} # default is mid-dormant (0.1% of ceiling)
# 
#     mactive_hourly = midphb_active*10^(0.099*Tc[1:senescenceperiod] - 5.14)} #10 % of high maintenance rate (med)
#   
#   
#   mphase1_daily = mphase1_hourly*24
#   mphase2_daily = mphase2_hourly*24
#   PHBused_phase1 = (mphase1_daily * SBgc)/gPHBtogC
#   PHBused_phase2 = (mphase2_daily * SBgc)/gPHBtogC
#   PHBused <- cumsum(c(PHBused_phase1,PHBused_phase2))
#   PHBused}
