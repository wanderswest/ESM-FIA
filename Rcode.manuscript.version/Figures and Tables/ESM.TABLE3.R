# EMPIRICAL SUCCESSION MAPPING Feb. 2015
# Climate/disturbance correlation table
# TRAVIS ANDREWS - tda2@me.com

 Do WITH wighted comparison correction for consistency!
 
 
#Load library, set options, simple fuctions
library(RColorBrewer)
library(sp)
library(data.table)
library(spam)
library(nlme)
library(maps)
library(ggplot2)
library(modeest)
library(effects)
library(plyr)
library(reshape2)
library(ggmap)
library(boot)
library(gplots)
library(shape)
library(SDMTools)
library(grid)
library(gridExtra)
library(plotrix)
library(foreach)
library(doParallel)
library(doMC)
registerDoMC(cores=2)
library(quantreg)
library(fields)
library(grDevices)
library(RCurl)
options(scipen=999) #remove sci notation
meanFunc <- function(x,i){mean(x[i], na.rm=TRUE)}


##Load data from GITHUB
 data.url <- "https://raw.githubusercontent.com/wanderswest/ESM-FIA/master/ESM.data.csv"
 
    ESM.data <- getURL(data.url)                
    ESM.data <- read.csv(textConnection(ESM.data), header = TRUE, sep = ",", quote="\"", dec=".")
RM5mergedfullcut1991.11.17.2014  <- as.data.table(ESM.data)

#or load filtered data from local drive
#RM5mergedfullcut1991.11.17.2014  <- as.data.table(read.csv("/Users/travis/GitHub/ESM-FIA/ESM.data.csv", header = TRUE, sep = ",", quote="\"", dec="."))

RM5mergedfullcut1991.11.17.2014 <-  subset(RM5mergedfullcut1991.11.17.2014, STDORGCD==0 ) # remove planted forests
S4swPREC <-  subset(RM5mergedfullcut1991.11.17.2014, STDORGCD==0 & cut==0 ) #exclude planted and harvested forests
F2 <- subset(S4swPREC)

print(c(round(sum(RM5mergedfullcut1991.11.17.2014$endsumTPA)/(1/0.166166884),0), "Trees surveyed"))




###############################################################################
# Extended data TABLE 3 - CLimate stats #################################################
#########
# Forests compared across same LAT LON grid
# Stats adjusted to account for minor differences in average tree densities. 
# Simplistic methods used because disturbance stats are inherently imperfect


#ratio of mean conditions with Null weighted spatially by Condition
meanFunc2 <- function(x,i){(mean(x$variCOND[i], na.rm=TRUE)/mean((x$variNULL[i]), na.rm=TRUE))} 

F2 <- subset(S4swPREC)  

#Convert variables to standard names and metric units
F2$DIAmean <- (F2$DIAendmean* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean <- (F2$DIAbeginmean *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum <- (F2$endsumTPA /0.404686)			#number of trees per hectare by plot at resurvey
F2$PREV_TPAsum <- (F2$startsumTPA /0.404686)	#number of trees per hectare by plot at initial survey (cm)
F2$PREV_STOCKINGmid <- F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey
 
#set grid cell size
a <- b <- 2  #2=0.5degrees, 0.5=2degrees
	F2$iLAT <- ((ceiling(F2$LAT*a))/a)-0.5 #+0.25  #round down to integer then add .25 ; switched /* removed ).25 addition
	F2$iLON <- ((floor(F2$LON*b))/b)+0.5 #-0.25  #round down to integer then subtract .25; switched /*  ).25 addition
	F2$LATLONYR <- (F2$iLAT*10000000000000) +(F2$iLON*(-100000000)) +F2$INVYR
	F2$LATLON <- as.factor(F2$iLAT*100000000+ F2$iLON*-1000) #LATLON identifyier

#Find climate variable mean and range for each grid cell
F2$TEMPminLL <- ave(F2$temp5, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$TEMPgrowminLL <- ave(F2$temp5grow, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$TEMP5chg <- F2$temp5-F2$TEMPminLL
F2$TEMP5growchg <- F2$temp5grow-F2$TEMPgrowminLL

F2$PRECminLL <- ave(F2$PREC5, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$PRECgrowminLL <- ave(F2$PREC5grow, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$PREC5chg <- F2$PREC5-F2$PRECminLL
F2$PREC5growchg <- F2$PREC5grow-F2$PRECgrowminLL

F2$SWdroughtLL <- ave(F2$SWdrought5, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$SWpluvialLL <- ave(F2$SWpluvial5, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$SWdrought5chg <- F2$SWdrought5-F2$SWdroughtLL
F2$SWpluvial5chg <- F2$SWpluvial5-F2$SWpluvialLL
F2$nsampTOT <- ave(F2$DIAmean>0,F2$LATLON, FUN=function(x) sum(x, na.rm=TRUE))

F2$Drysummer <- ifelse(F2$SWdrought5chg <= quantile(F2$SWdrought5chg, 0.17, na.rm=T), F2$startsumTPA/6.018,0)
F2$Wetsummer <- ifelse(F2$SWpluvial5chg >= quantile(F2$SWpluvial5chg, 0.83, na.rm=T), F2$startsumTPA/6.018,0)
F2$Dryperiod <- ifelse(F2$PREC5growchg <= quantile(F2$PREC5growchg, 0.17, na.rm=T), F2$startsumTPA/6.018,0)
F2$Wetperiod <- ifelse(F2$PREC5growchg >= quantile(F2$PREC5growchg, 0.83, na.rm=T), F2$startsumTPA/6.018,0)
F2$SPCDdiversity2 <- F2$SPCDdiversity-F2$PREV_SPCDdiversity

#set disturbance conditions to binary
#F2$fire <- ifelse(F2$fire>0, 1,0) # F2$carbon1sum-F2$PREVcarbon1sum,0) #ifelse(F2$fire>0, 1,0)
#F2$insect <- ifelse(F2$insect>0, 1,0) # F2$carbon1sum-F2$PREVcarbon1sum,0) # 1,0)
#F2$disease <- ifelse(F2$disease>0, 1,0) # F2$carbon1sum-F2$PREVcarbon1sum,0) # 1,0)
#F2$animal <- ifelse(F2$animal>0,  1,0) #F2$carbon1sum-F2$PREVcarbon1sum,0) # 1,0)
#F2$weather <- ifelse(F2$weather>0,  1,0) #F2$carbon1sum-F2$PREVcarbon1sum,0) # 1,0)
#F2$SPCDdiversity2 <- ifelse(F2$SPCDdiversity2>0, 1,0)

vari <- c("SWdrought5chg",  "SWpluvial5chg", "PREC5growchg", "TEMP5growchg") 
vari1 <- subset(F2, select=c(SWdrought5chg, SWpluvial5chg, PREC5growchg, TEMP5growchg)) 
vari2 <- c("fire", "insect", "disease", "animal", "weather", "Drysummer", "Wetsummer", "Dryperiod", "Wetperiod", "SPCDdiversity2")

###ADJUST TOP/bottom QUANTILE HERE
quantile1<-0.17 #0.03 #0.17

F5b <- NULL
for(v2 in 1:5)
{
for(v in 1:4)
{
for(c in 1:2)
{
if(c==1 & vari[v]=="SWdrought5chg"){
	f22 <- subset(F2, SWdrought5chg < quantile(F2$SWdrought5chg, quantile1, na.rm=TRUE) & SWpluvial5chg<= quantile(F2$SWpluvial5chg, 0.85, na.rm=TRUE))
	f23 <- subset(F2, SWdrought5chg <= quantile(F2$SWdrought5chg, 0.7, na.rm=TRUE) & SWdrought5chg >= quantile(F2$SWdrought5chg, 0.3, na.rm=TRUE) & SWpluvial5chg<= quantile(F2$SWpluvial5chg, 0.83, na.rm=TRUE))
	 }	#
if(c==2 &  vari[v]=="SWpluvial5chg"){
	f22 <- subset(F2, SWpluvial5chg > quantile(F2$SWpluvial5chg, 1-quantile1, na.rm=TRUE) & SWdrought5chg>= quantile(F2$SWdrought5chg, 0.15, na.rm=TRUE)) 
	f23 <- subset(F2, SWpluvial5chg <= quantile(F2$SWpluvial5chg, 0.7, na.rm=TRUE) & SWpluvial5chg >= quantile(F2$SWpluvial5chg, 0.3, na.rm=TRUE)& SWdrought5chg>= quantile(F2$SWdrought5chg, 0.17, na.rm=TRUE)) 
	}	

if(c==1 & vari[v]=="SWpluvial5chg"){next}
if(c==2 & vari[v]=="SWdrought5chg"){next}

if(c==1 & v>2){
	f22 <- subset(F2, get(vari[v])< quantile(vari1[[v]], quantile1, na.rm=TRUE)) 
	f23 <- subset(F2, get(vari[v])<= quantile(vari1[[v]], 0.7, na.rm=TRUE) & get(vari[v])>= quantile(vari1[[v]], 0.3, na.rm=TRUE))
	}	#0.17 & 0.03 set slimate extreme bottom percentile
if(c==2 & v>2){
	f22 <- subset(F2, get(vari[v])> quantile(vari1[[v]], 1-quantile1, na.rm=TRUE)) 
	f23 <- subset(F2, get(vari[v])<= quantile(vari1[[v]], 0.7, na.rm=TRUE) & get(vari[v])>= quantile(vari1[[v]], 0.3, na.rm=TRUE))
	}	#0.83 & 0.965 set slimate extreme top percentile

F33 <- subset(f23, LATLON %in% (f22$LATLON)) # not neccesaary bc weighting will eliminate latlons not in f22!

f22$Csamp <- (ave(f22$DIAmean>0,f22$LATLON, FUN=function(x) sum(x, na.rm=TRUE)))/f22$nsampTOT #weighting factor of number of Cond data/total data in each grid cell
F33$Nsamp <- (ave(F33$DIAmean>0,F33$LATLON, FUN=function(x) sum(x, na.rm=TRUE)))/F33$nsampTOT #weighting factor of number of NULL data/total data in each grid cell

f24 <- subset(f22, !duplicated(LATLON), select=c(LATLON, Csamp))
F33 <- merge(F33, f24, all.x=TRUE, by="LATLON")
F33$CtoNratio <- 1 #F33$Csamp/F33$Nsamp  #needs to be for each disturb type tooo?

variCOND <- subset(f22, select=c("fire", "insect", "disease", "animal", "weather", "Drysummer", "Wetsummer", "Dryperiod", "Wetperiod", "unknowndamage"))
variNULL <- subset(F33, select=c("fire", "insect", "disease", "animal", "weather", "Drysummer", "Wetsummer", "Dryperiod", "Wetperiod", "unknowndamage", "CtoNratio"))

#find confidence intervals for null model disturbance freq. 
ratioFunc <- function(x, i){mean(x[i], na.rm=T)}  #mean(variCOND[[v2]], na.rm=T)/ #ratio not used instead just CI for NULL

variCOND1<-variCOND[[v2]]/f22$startsumTPA #*variNULL$CtoNratio

ratioboot2 <- boot(variCOND1, ratioFunc, R=10000) #variCOND[[v2]]
NULL.ci <- boot.ci(ratioboot2, conf=0.95, type="norm")
errorlow <- NULL.ci$normal[2]
errorhigh <- NULL.ci$normal[3]

#NULLTPA<- (std.error(F33$startsumTPA, na.rm=T))*1.96
#CONDTPA<- (std.error(f22$startsumTPA, na.rm=T))*1.96
NULLTPA<- (mean(F33$startsumTPA, na.rm=T))
CONDTPA<- (mean(f22$startsumTPA, na.rm=T))
TPAadj<- ((CONDTPA-NULLTPA)/NULLTPA)*100  # not needed if Cond to Null ratio (CtoNratio) used

#subset data based on condition to get actual ratio in Null climate and Conditional climates rather than Bootstrapped mean is approx. the same
#NullVarRatiod <- mean(variNULL[[v2]]*variNULL$CtoNratio)
#CondVarRatiod <- mean(variCOND[[v2]]) 

# average percent of plot disturbed
NullVarRatiod <- mean(variNULL[[v2]]/F33$startsumTPA)
CondVarRatiod <- mean(variCOND[[v2]]/f22$startsumTPA)
TPAadj<-0

rateRatiod <- (((CondVarRatiod/NullVarRatiod)*100)-100)-TPAadj
#rateRatiod.CIupper <- (((CondVarRatiod/errorhigh)*100)-100)-TPAadj
#rateRatiod.CIlower <- (((CondVarRatiod/errorlow)*100)-100)-TPAadj
rateRatiod.CIupper <- (((errorlow/NullVarRatiod)*100)-100)-TPAadj
rateRatiod.CIlower <- (((errorhigh/NullVarRatiod)*100)-100)-TPAadj

#NullVarRatiodex.unknown <- mean(variNULL[[v2]][variNULL$unknowndamage==0]*variNULL$CtoNratio[variNULL$unknowndamage==0]) #
#CondVarRatiodex.unknown <- mean(variCOND[[v2]][variCOND$unknowndamage==0])
NullVarRatiodex.unknown <- mean(variNULL[[v2]][variNULL$unknowndamage==0]/F33$startsumTPA[variNULL$unknowndamage==0])
CondVarRatiodex.unknown <- mean(variCOND[[v2]][variCOND$unknowndamage==0]/f22$startsumTPA[variCOND$unknowndamage==0]) 
TPAadj<-0

rateRatiodex.unknown <- (((CondVarRatiodex.unknown/NullVarRatiodex.unknown)*100)-100)-TPAadj

rateRatiod.CIupper2<- ifelse(rateRatiod.CIupper>=rateRatiod.CIlower, NA, rateRatiod.CIupper)
rateRatiod.CIlower2<- ifelse(rateRatiod.CIlower<=rateRatiod.CIupper, NA, rateRatiod.CIlower)

#sig <- ifelse((CondVarRatiod > errorhigh | CondVarRatiod < errorlow), "***", "")
sig <- ifelse(c(rateRatiod.CIupper2 > 0 | rateRatiod.CIlower2 < 0), "***", "")
sig <- ifelse(c(is.na(rateRatiod.CIupper2) | is.na(rateRatiod.CIlower2)), NA, sig)

Clim <- ifelse(vari[v]=="PREC5growchg" & c==2, "5-yr Wet period", NA)
Clim <- ifelse(vari[v]=="PREC5growchg" & c==1, "5-yr Dry period", Clim)
Clim <- ifelse(vari[v]=="TEMP5growchg" & c==1, "Cool period", Clim)
Clim <- ifelse(vari[v]=="TEMP5growchg" & c==2, "Warm period", Clim)
Clim <- ifelse(vari[v]=="SWpluvial5chg" & c==2, "Wet summer", Clim)
Clim <- ifelse(vari[v]=="SWdrought5chg" & c==1, "Dry summer", Clim)

n <- sum(variCOND[[v2]]>0) #==1)
n2 <- sum(f22$startsumTPA>0 ) #/6.07)
nnonAC <- sum(variNULL[[v2]]>0) #==1)
n2nonAC <- sum(F33$startsumTPA>0 ) #/6.07)

f34 <- c(vari2[v2],Clim, round(n,0), round(n2,0), round(nnonAC,0), round(n2nonAC), round(CondVarRatiod,6), round(NullVarRatiod,6), round(rateRatiod,4), round(rateRatiodex.unknown,2), round(rateRatiod.CIupper2,2), round(rateRatiod.CIlower2,2), sig) #, round(Vgrow2,2), round(Vgrow,2), error2)
F5b <- rbind(F5b, f34)
print(c(v2, v, c))
}
}
}

F6b <- as.data.table(F5b)
setnames(F6b, c("Disturbance type", "Climate condition", "n","COND n", "n.nonAC.trees", "nonAC.n.plots", "CondVarRatiod", "NullVarRatiod", "%rateRatiod","%rateratiodex.unknown", "CI.lower","CI.upper", "p < 0.05") )#, "Vgrow2","Vdead", "error2"))
print((F6b))
#write.csv((F6b))
 
 # Productivity was -62.5% lower in forests during a warm period 
 #weather disturbances were associated with -62% reduction in productivity during warm periods
 # Forests carbon change with disease with and without a wet period (witout latlon weighted average)
 #do confidence interval for %number of trees dead change and see if overlap with zero for significance 
 #subtract percent change expected because of differences in TPA and taken from same LATLON
 #should not need cond versesus null TPA correction if ratio to mortality within each !
#################################################################################################
##END MAIN FIGURES and climate conditions supplimentary figures


