
#EMPIRICAL SUCCESSION MAPPING Feb. 2015
#TRAVIS ANDREWS - TDA210@lehigh.edu

#### Figure 2 - by percentage difference

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



 data.url <- "https://raw.githubusercontent.com/wanderswest/ESM-FIA/master/ESM.data.csv"
 
    ESM.data <- getURL(data.url)                
    ESM.data <- read.csv(textConnection(ESM.data), header = TRUE, sep = ",", quote="\"", dec=".")
RM5mergedfullcut1991.11.17.2014  <- as.data.table(ESM.data)



##Load data from GITHUB
 # data.url <- "https://raw.githubusercontent.com/wanderswest/ESM-FIA/master/ESM.data.csv"
 
    # ESM.data <- getURL(data.url)                
    # ESM.data <- read.csv(textConnection(ESM.data), header = TRUE, sep = ",", quote="\"", dec=".")
# RM5mergedfullcut1991.11.17.2014  <- as.data.table(ESM.data)

#or load filtered data from local drive
#RM5mergedfullcut1991.11.17.2014  <- as.data.table(read.csv("/Users/travis/GitHub/ESM-FIA/ESM.data.csv", header = TRUE, sep = ",", quote="\"", dec="."))
#RM5mergedfullcut1991.11.17.2014  <- as.data.table(read.csv("/Users/travis/Desktop/ESMthin/Null2015PFTs.csv", header = TRUE, sep = ",", quote="\"", dec="."))

RM5mergedfullcut1991.11.17.2014 <-  subset(RM5mergedfullcut1991.11.17.2014, STDORGCD==0 ) # remove planted forests
S4swPREC <-  subset(RM5mergedfullcut1991.11.17.2014, STDORGCD==0 & cut==0 ) #exclude planted and harvested forests
F2 <- subset(S4swPREC)

print(c(sum(RM5mergedfullcut1991.11.17.2014$endsumTPA)/(1/0.166166884), "Trees resurveyed"))







#################################################################################################################
################ FIGURE 2 - Climate Condition DIVERGENCE ########################################################
####
##

#palette(adjustcolor(((rich.colors(60))), transform=diag(c(0.85,0.95,0.75,1)))) # alpha.f=0.15 #60

#Climate Colors:
pt <- "steelblue1" #pluvial trend
pe <- "royalblue"	#pluvial event
dt <- "goldenrod1" #drought trend
de <- "orange2"  	#drought event
nt <- "grey55"
ne <- "grey35"
wt <- "peachpuff"
ct <- "slategray1"
errorcol <- "grey75"
errorlwd <- 0.53
ag1 <- 25


#S4swPREC <- as.data.table(read.csv("/Volumes/m-z/tda210/USFS/RM5mergedfullcut1991.11.17.2014.csv", header = TRUE, sep = ",", quote="\"", dec="."))
F2 <- subset(S4swPREC)  #, SISP<300) 


undisturbed<-0
deciduous<-0
coniferous<-0
#dataset for undisturbed analysis Fig S2 <-exclude forests plots with identified disturbances 
if(undisturbed==1){F2 <- subset(F2,  disease==0 & weather==0 & fire==0 &  animal==0 & insect==0)  }
if(deciduous==1){F2 <- subset(F2, SISP>=300 )} #(EHstocksum+Evergreenstocksum+Hydricstocksum+ LHstocksum+NMHstocksum+SMHstocksum)>(LCstocksum +MCstocksum+ NPstocksum+ SPstocksum+ CDstocksum)
if(coniferous ==1){F2 <- subset(F2,  SISP<300)  }

#Convert variables to standard names and metric units
F2$DIAmean <- (F2$DIAendmean* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean <- (F2$DIAbeginmean *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum <- (F2$endsumTPA /0.404686)			#number of trees per hectare by plot at resurvey
F2$PREV_TPAsum <- (F2$startsumTPA /0.404686)	#number of trees per hectare by plot at initial survey (cm)
F2$PREV_STOCKINGmid <- F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey
#F2$TPAchg <- (F2$TPAsum-F2$PREV_TPAsum)/F2$REMPER
F2$SISPchg<-ifelse(F2$SISP>=300, 1, ifelse(F2$SISP<300, -1, NA) )
F2$DIAchgmean<- F2$DIAmean - F2$PREV_DIAmean 


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
F2$TEMP5growchg[is.na(F2$TEMP5growchg)==T] <- 0
F2$TEMP5growchg[F2$TEMP5growchg==0] <- rnorm(sum(F2$TEMP5growchg==0),0,sd=0.01) #make median 0 value a distribution for subseting
# note random distribution of median climate values will make results slightly different each run

#F2$PREC5grow <- F2$PREC5grow*92/3/10 #convert mm/day to cm per month
F2$PRECminLL <- ave(F2$PREC5, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$PRECgrowminLL <- ave(F2$PREC5grow, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$PREC5chg <- F2$PREC5-F2$PRECminLL
F2$PREC5growchg <- F2$PREC5grow-F2$PRECgrowminLL
F2$PREC5growchg[F2$PREC5growchg==0] <- rnorm(sum(F2$PREC5growchg==0),0,sd=0.01) #make median 0 value a distribution for subseting

#make median 0 value a distribution for subseting (but not random for consistent results!)
F2$SWdroughtLL <- ave(F2$SWdrought5, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$SWpluvialLL <- ave(F2$SWpluvial5, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$SWdrought5chg <- F2$SWdrought5-F2$SWdroughtLL
F2$SWpluvial5chg <- F2$SWpluvial5-F2$SWpluvialLL
F2$SWpluvial5chg[F2$SWpluvial5chg==0] <- rnorm(sum(F2$SWpluvial5chg ==0),0,sd=0.01)  #seq(0.01,-0.01, ((0.02/((sum(F2$SWpluvial5chg==0)-1)))))
F2$SWdrought5chg[F2$SWdrought5chg==0 ] <- rnorm(sum(F2$SWdrought5chg ==0),0,sd=0.01)  #seq(0.01, -0.01, ((0.002/((sum(F2$SWdrought5chg==0)/10)-1))))

F2$nsampTOT <- ave(F2$PREV_DIAmean>0,F2$LATLON, FUN=function(x) sum(x, na.rm=TRUE))


vari <- c("PREC5growchg", "TEMP5growchg", "SWpluvial5chg", "SWdrought5chg") #"PREC5growchg", "TEMP5growchg"
vari1 <- subset(F2, select=c(PREC5growchg, TEMP5growchg, SWpluvial5chg, SWdrought5chg))
vari2 <- c("5-yr JJA precip. anomaly (mm/day)", "5-yr JJA temp. anomaly (C)", "1-yr JJA soil moist. anomaly (z-score)", "1-yr JJA soil moist. anomaly (z-score)")

seq1<- c(0.2, 0.4, NA, 0.66, 0.86, 0.96, 1, NA)

if(coniferous==1 | deciduous==1){seq1<- c(0.2, 0.4, NA, 0.5, 0.75, 0.9, 1)} #broader quantiles for smaller dataset
#(100-(c(0.2, 0.4, NA, 0.5, 0.75, 0.9, 1)*100))/2

qt1a <- c(quantile(F2$PREC5growchg[F2$PREC5growchg<0],(1-rev(seq1))))  
 qt2a <- c(quantile(F2$PREC5growchg[F2$PREC5growchg>0], seq1)) 
 qt1b <- c(quantile(F2$TEMP5growchg[F2$TEMP5growchg<0],(1-rev(seq1)))) 
 qt2b <- c(quantile(F2$TEMP5growchg[F2$TEMP5growchg>0], seq1))  
 qt3 <- c(quantile(F2$SWdrought5chg[F2$SWdrought5chg<0],(1-rev(seq1)))) 
 qt4 <- c(quantile(F2$SWpluvial5chg[F2$SWpluvial5chg>0], seq1)) 
int <- 1


F5 <- NULL
for(v in 1:3)
{
	if(v==1){
		 qt1<- qt1a
		 qt2<- qt2a
	}
	if(v==2){
		 qt1<- qt1b
		 qt2<- qt2b
	}

	
F4<-foreach(c = 1:16, .combine=rbind) %dopar%
{
	if(c<=8 & v<3)
	{
		f22 <- subset(F2, get(vari[v])>= qt2[c] & get(vari[v])<=  qt2[c+int]) 		
		f23 <- subset(F2, get(vari[v])>= quantile(vari1[[v]], 0.3, na.rm=TRUE) & get(vari[v])<= quantile(vari1[[v]], 0.7, na.rm=TRUE), select=c(PLT_CN)) 
		ccond <- vari[v]
		col3=ifelse(qt2[c]>= 0.2, 3, 2)
		qta <- qt2[c]
		qtb <- qt2[c+int]
		}
	if(c>8 & v<3)
	{ 
		d=c-8
		f22 <- subset(F2, get(vari[v])<= qt1[d] & get(vari[v])>= qt1[d-int]) #		
		f23 <- subset(F2, get(vari[v])>= quantile(vari1[[v]], 0.3, na.rm=TRUE) & get(vari[v])<= quantile(vari1[[v]], 0.7, na.rm=TRUE), select=c(PLT_CN)) 
		ccond <- vari[v]
		col3=ifelse(qt1[d]<= -0.1, 1, 2)
		qtb <- qt1[d]
		qta <- qt1[d-int]
	}
		
	
	if(c<=8 & v==3)
	{
			f22 <- subset(F2, get(vari[v])>= qt4[c] & get(vari[v])<= qt4[c+int] & SWdrought5chg>= quantile(F2$SWdrought5chg, 0.15, na.rm=TRUE)) #exclude droughts
			f23 <- subset(F2, get(vari[v])<= quantile(vari1[[v]], 0.7, na.rm=TRUE) & get(vari[v])>= quantile(vari1[[v]], 0.3, na.rm=TRUE) & SWdrought5chg>= quantile(F2$SWdrought5chg, 0.15, na.rm=TRUE), select=c(PLT_CN)) 
		ccond <- "Wet summer"
		col3=ifelse(qt4[c]>=0.2, 3, 2)
		qta <- qt4[c]
		qtb <- qt4[c+int]
	}
	if(c>8 & v==3)
	{ 
		d=c-8
		f22 <- subset(F2, get(vari[v+1])< qt3[d] & get(vari[v+1])>= qt3[d-int] & SWpluvial5chg<= quantile(F2$SWpluvial5chg, 0.85, na.rm=TRUE)) #exclude wet summers  
		f23 <- subset(F2, get(vari[v+1])<= quantile(vari1[[v+1]], 0.7, na.rm=TRUE) & get(vari[v+1])>= quantile(vari1[[v+1]], 0.3, na.rm=TRUE) & SWpluvial5chg<= quantile(F2$SWpluvial5chg, 0.85, na.rm=TRUE), select=c(PLT_CN)) 
		ccond <- "Dry summer"
		col3=ifelse(qt3[d]<= -0.2, 1, 2)
		qtb <- qt3[d]
		qta <- qt3[d-int]
	}
	

f22$nsamp <- (ave(f22$PREV_DIAmean>0,f22$LATLON, FUN=function(x) sum(x, na.rm=TRUE)))/f22$nsampTOT #sample proportion of each grid cell 
#print(c(table(f22$nsamp<0.75), col3))

#f22$col3 <- col3
PLT_CNsamp1 <- f23$PLT_CN   

f24 <- subset(f22, !duplicated(LATLON), select=c(LATLON, nsamp))

F33 <- subset(F2, LATLON %in% (f22$LATLON)) # not neccesaary bc weighting will eliminate latlons not in f22!
F33 <- merge(F33, f24, all.x=T, by="LATLON")
F33 <- subset(F33, PLT_CN %in% (PLT_CNsamp1)) # not neccesaary bc weighting will eliminate latlons not in f22!

F33$aPREVSTOCKDIAbin <- floor(F33$PREV_STOCKINGmid/10)*10000+(floor(F33$PREV_DIAmean/2)*2)
F33$n1 <- 1
F33$aPREVSTOCKDIAbinlength <- ave(F33$n1, F33$aPREVSTOCKDIAbin, FUN=function(x) sum(x, na.rm=TRUE))
F33$nbin <- ave(F33$nsamp, F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE)) #average AC proportion in stock/dia bins

F33 <- subset(F33, aPREVSTOCKDIAbinlength> 10) #Null model needs at least 10 comparison plots at each diameter/stocking bin

##Carbon Model
c6 <- (((F33$carbon1sum-F33$PREVcarbon1sum)* 0.00112085116 )/F33$REMPER)
F33$avgCchg <- (ave(c6*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE)))/F33$nbin #multiply by AC proportion divided by mean proportion in each stocking/dia bin

##Carbon Model 5 inch and greater trees
c5 <- (((F33$carbon5sum-F33$PREVcarbon5sum)* 0.00112085116 )/F33$REMPER)
F33$avgC5chg <- (ave(c5*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE)))/F33$nbin

##radial growth  
rx5 <- (F33$DIAchgmean/F33$REMPER) ##
F33$avgx5 <- ave(rx5*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE))/F33$nbin

##radial growth
rx5grow <- (F33$DIAgrowmean/F33$REMPER)/2 ## convert from diameter to radial growth
F33$avgx5grow <- ave(rx5grow*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE))/F33$nbin

##Mortality diameter mean
mortrx5 <- F33$mortDIAmean  #
F33$avgmortx5 <- ave(mortrx5*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE))/F33$nbin

##Density Model
y6 <- (F33$TPAsum-F33$PREV_TPAsum)/F33$REMPER  #
F33$avgTPAchg <- ave(y6*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE))/F33$nbin

##SAPLINGS Density Model
sap6 <- ((F33$TPAsum05/0.404686)-(F33$PREVTPAsum05/0.404686))/F33$REMPER  # convert to hectares
F33$sapavgTPAchg <- ave(sap6*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE))/F33$nbin
##SAPLINGS mortality /recruitment
F33$sapavgMORT <- ave(F33$sapMortsum*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE))/F33$nbin #sapling mortality per plot
F33$sapavgRECR <- ave(F33$sapRecruitsum*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE))/F33$nbin #sapling plot recuritment

#SISP ratio
F33$sp1<- mean(F33$SISPchg, na.rm=T) #F33$nsamp
F33$spavgchg <- ave(F33$sp1,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE)) #/F33$nbin

#Merge null and conditional model data
F34 <- subset(F33, !duplicated(aPREVSTOCKDIAbin), select=c(aPREVSTOCKDIAbin, avgCchg, avgx5, avgx5grow, avgTPAchg, sapavgTPAchg, avgC5chg, avgmortx5, sapavgMORT, sapavgRECR))
f22$aPREVSTOCKDIAbin <- floor(f22$PREV_STOCKINGmid/10)*10000+(floor(f22$PREV_DIAmean/2)*2)
f22 <- merge(f22, F34, by="aPREVSTOCKDIAbin") #, all.x=T


#carbon 1 and 5 inch
c8 <- (((f22$carbon1sum-f22$PREVcarbon1sum)* 0.00112085116 )/f22$REMPER)
f22$xcarbon <- (c8-f22$avgCchg)
c55 <- (((f22$carbon5sum-f22$PREVcarbon5sum)* 0.00112085116 )/f22$REMPER)
f22$xcarbon5 <- ((c55)-f22$avgC5chg)


##mean diameter change  
x8 <- (f22$DIAchgmean/f22$REMPER)  # 
f22$x5 <-  (x8-f22$avgx5)

##radial growth  
x8grow <- (f22$DIAgrowmean/f22$REMPER)/2  # convert from diameter to radial growth
f22$x5grow <-  (x8grow-f22$avgx5grow)

##mortality size
mortx8 <- f22$mortDIAmean #
f22$mortx5 <-  (mortx8-f22$avgmortx5)

##Density Model
TPA <- (f22$TPAsum-f22$PREV_TPAsum)/f22$REMPER  #
f22$TPAchg <-  (TPA-f22$avgTPAchg)

##SAPLING Density Model
sapTPA <- ((f22$TPAsum05/0.404686)-(f22$PREVTPAsum05/0.404686))/f22$REMPER  #
f22$sapTPAchg <-  (sapTPA-f22$sapavgTPAchg)
f22$sapMORTchg <- (f22$sapMortsum/f22$sapavgMORT)*100 #sapling mortaltiy percent difference from null
f22$sapRECRchg <- (f22$sapRecruitsum/f22$sapavgRECR)*100 #sapling recruitment percent difference from null

##SISP ratio
# sp22 <- mean(f22$SISPchg-f22$sapavgTPAchg, na.rm=T)

#calculate mean condition for x-axis based on v loop
if(length(f22$xcarbon)>= 10)
{
if(v==2)
{
	x3 <- mean(f22$TEMP5growchg, na.rm=TRUE)
	minx3<-min(f22$TEMP5growchg, na.rm=TRUE)
	maxx3<-max(f22$TEMP5growchg, na.rm=TRUE)
	SDx3<-sd(f22$TEMP5growchg, na.rm=TRUE)
}
if(v==1)
{
	x3 <- mean(f22$PREC5growchg, na.rm=TRUE)
	minx3<-min(f22$PREC5growchg, na.rm=TRUE)
	maxx3<-max(f22$PREC5growchg, na.rm=TRUE)
	SDx3<-sd(f22$PREC5growchg, na.rm=TRUE)

}
if(v==3)
{
	x3 <- ifelse(c<=8,mean(f22$SWpluvial5chg, na.rm=TRUE), mean(f22$SWdrought5chg, na.rm=TRUE))
	minx3 <- ifelse(c<=8,min(f22$SWpluvial5chg, na.rm=TRUE), min(f22$SWdrought5chg, na.rm=TRUE))
	maxx3 <- ifelse(c<=8,max(f22$SWpluvial5chg, na.rm=TRUE), max(f22$SWdrought5chg, na.rm=TRUE))
	SDx3 <- ifelse(c<=8,sd(f22$SWpluvial5chg, na.rm=TRUE), sd(f22$SWdrought5chg, na.rm=TRUE))


}

#calculate y-axis variables
#USxcarbon<- mean(F2$carbon1sum-F2$PREVcarbon1sum)* 0.00112085116/mean(F2$REMPER)
y3 <- mean(f22$xcarbon, na.rm=TRUE) 
y3pctdif<-(mean(f22$xcarbon, na.rm=TRUE)/mean(f22$avgCchg, na.rm=TRUE))*100   #
#y3<- (y3pctdif/100)*USxcarbon
y3start <- mean(f22$avgCchg, na.rm=T) #null model zero line for plotting

y5 <- mean(f22$xcarbon5, na.rm=TRUE)  #
y5pctdif<-(mean(f22$xcarbon5, na.rm=TRUE)/mean(f22$avgC5chg, na.rm=TRUE))*100   #
y5start <- mean(f22$avgC5chg, na.rm=T) #null model zero line for plotting


#USTPA3 <- mean((F2$TPAsum-F2$PREV_TPAsum)/F2$REMPER)
TPA3pctdif<-(mean(f22$TPAchg, na.rm=TRUE)/mean(f22$avgTPAchg, na.rm=TRUE))*100   #
TPA3 <-mean(f22$TPAchg, na.rm=TRUE)
TPA3start <- mean(f22$avgTPAchg, na.rm=T) #null model zero line for plotting

#USsapTPA3 <- mean(((F2$TPAsum05/0.404686)-(F2$PREVTPAsum05/0.404686))/F2$REMPER)
sapTPA3pctdif <- (mean(f22$sapTPAchg, na.rm=TRUE)/ mean(f22$sapavgTPAchg, na.rm=T))*100#
sapTPA3 <- mean(f22$sapTPAchg, na.rm=TRUE) #
sapTPA3start <- mean(f22$sapavgTPAchg, na.rm=T)

#USx5 <- mean(F2$DIAchgmean/F2$REMPER, na.rm=T)/2  # convert from diameter to radial growth
x5pctdif <- (mean(f22$x5, na.rm=TRUE)/mean(f22$avgx5, na.rm=T))*100  #
x5<- mean(f22$x5, na.rm=TRUE)
x5start <- mean(f22$avgx5, na.rm=T)

#USx5grow <- mean(F2$DIAgrowmean/F2$REMPER, na.rm=T)/2  # convert from diameter to radial growth
x5growpctdif <- (mean(f22$x5grow, na.rm=TRUE)/mean(f22$avgx5grow, na.rm=T))*100  #
x5grow <- mean(f22$x5grow, na.rm=TRUE) #
x5growstart <- mean(f22$avgx5grow, na.rm=T)

mortx5 <- mean(f22$mortx5, na.rm=TRUE) #
mortx5start <-   mean(mortx8, na.rm=T)
sapMORT <- mean(f22$sapMORTchg, na.rm=T)
sapRECR <- mean(f22$sapRECRchg, na.rm=T)

#calc y-axis error
R2 <- 10 #repetitions set to 10 for faster looping
enderror <- boot(f22$xcarbon, meanFunc, R=R2)
enderrory3 <- sqrt(var(enderror$t))*1.96
enderror5 <- boot(f22$xcarbon5, meanFunc, R= R2)
enderrory5 <- sqrt(var(enderror5$t))*1.96
enderror2 <- boot(f22$TPAchg, meanFunc, R= R2)
enderrorTPA <- sqrt(var(enderror2$t))*1.96
enderrorsap <- boot(f22$sapTPAchg, meanFunc, R= R2)
enderrorsapTPA <- sqrt(var(enderrorsap$t))*1.96
enderror4 <- boot(f22$x5, meanFunc, R= R2)
enderrorx5 <- sqrt(var(enderror4$t))*1.96
enderrorx5grow4 <- boot(f22$x5grow, meanFunc, R= R2)
enderrorx5grow <- sqrt(var(enderrorx5grow4$t))*1.96
enderrorsapMORT <- std.error(f22$sapMORTchg, na.rm=T)
enderrorsapRECR <- std.error(f22$sapRECRchg, na.rm=T)

# enderrorsp2error <- boot((f22$SISPchg-f22$sapavgTPAchg), meanFunc, R= R2)
# enderrorsp22error <- sqrt(var(enderrorsp2error$t))*1.96

n <- length(f22$xcarbon) #conditional model n
nnull <- length(F33$avgx5) #null model n
#print(c(n, ccond, c, nnull))

f23 <- c(x3, y3, y5, TPA3, sapTPA3, x5, x5grow, enderrory3, enderrory5, enderrorTPA, enderrorsapTPA, enderrorx5, enderrorx5grow, n, v, col3, y3start, y5start, TPA3start, x5start, sapTPA3start, x5growstart, sapMORT, sapRECR, enderrorsapMORT, enderrorsapRECR, y3pctdif, y5pctdif, TPA3pctdif, sapTPA3pctdif ,x5pctdif, x5growpctdif, minx3, maxx3, SDx3)

}
}
F5 <- rbind(F5,F4)
}

F6 <- as.data.table(F5)
setnames(F6, c("x3","y3","y5", "TPA","sapTPA", "x5", "x5grow", "errory3","errory5", "errorTPA","errorsapTPA", "errorx5","errorx5grow", "n", "v", "col3", "y3start", "y5start", "TPA3start", "x5start", "sapTPA3start", "x5growstart", "sapMORTpct", "sapRECRpct", "enderrorsapMORT", "enderrorsapRECR", "y3pctdif", "y5pctdif", "TPA3pctdif", "sapTPA3pctdif", "x5pctdif", "x5growpctdif", "minx3", "maxx3", "SDx3"))
F6 <- subset(F6, n>100) 

y3start1 <- mean(F6$y3start)
y5start1 <- mean(F6$y5start)
TPA3start1 <-mean(F6$TPA3start) #mean(F2$TPAchg) #
sapTPA3start1 <-mean(F6$sapTPA3start)
x5start1 <- mean((F6$x5start))
x5growstart1 <- mean((F6$x5growstart)*5) #*10)

#convert difference btw ac and nonac to relative difference from average
F6$y3<- y3start1*(F6$y3pctdif/100)
F6$TPA<- TPA3start1*(F6$TPA3pctdif/100)
F6$sapTPA<- sapTPA3start1*(F6$sapTPA3pctdif/100)
F6$x5<- x5start1*(F6$x5pctdif/100)
F6$x5grow<- x5growstart1*(F6$x5growpctdif/100)

F6$errory3<-(((y3start1-F6$y3start)/F6$y3start)*F6$errory3)+F6$errory3
F6$errory5<-(((y5start1-F6$y5start)/F6$y5start)*F6$errory5)+F6$errory5
F6$errorTPA<-(((TPA3start1-F6$TPA3start)/F6$TPA3start)*F6$errorTPA)+F6$errorTPA
F6$errorsapTPA<-(((sapTPA3start1-F6$sapTPA3start)/F6$sapTPA3start)*F6$errorsapTPA)+F6$errorsapTPA
F6$errorx5<-(((x5start1-F6$x5start)/F6$x5start)*F6$errorx5)+F6$errorx5
F6$errorx5grow<-((((mean(F6$x5growstart)-F6$x5growstart)/F6$x5growstart)*F6$errorx5grow)+F6$errorx5grow)*5


#convert radial growth from cm to mm
#F6$x5grow <- F6$x5grow*5 #*10
#F6$errorx5grow <- F6$errorx5grow*5 #*10
#convert to percent difference
F6$sapMORTpct <-  (F6$sapMORTpct-100)
F6$sapRECRpct <-  (F6$sapRECRpct-100)
F6$sapCpct <-  ((F6$y3-F6$y5)/F6$y3)*100


#calc plot y-scales
sapscale  <-  30 #mean(F6$sapTPA3start)/mean(F6$TPA3start) 
Cymax  <-  max(F6$y3+F6$errory3)
Cymin  <-  min(F6$y3-F6$errory3)
#TPAymax <- max(F6$sapTPA/sapscale+F6$errorsapTPA/sapscale, na.rm=T)
TPAymax1 <- (max(F6$sapTPA+F6$errorsapTPA))/sapscale
TPAymax2 <- (max(F6$TPA+F6$errorTPA))
TPAymax <- max(TPAymax1, TPAymax2)
TPAymin2 <- (min(F6$TPA-F6$errorTPA))
TPAymin1 <- (min(F6$sapTPA-F6$errorsapTPA))/sapscale
TPAymin <- min(TPAymin1, TPAymin2)
DIAmax1 <- (max(F6$x5grow+F6$errorx5grow))
DIAmin1 <- (min(F6$x5grow-F6$errorx5grow))
DIAmax2 <- max(F6$x5+F6$errorx5)
DIAmin2 <- min(F6$x5-F6$errorx5)
DIAmin <- min(DIAmin2, DIAmin1)
DIAmax <- max(DIAmax2, DIAmax1)
sapmax <- max(F6$sapMORTpct+F6$enderrorsapMORT)
sapmin <- min(F6$sapMORTpct-F6$enderrorsapMORT)

#plot only sapling data that is significantly different with scaling than larger trees. 
#F6$sapTPA <- ifelse(((F6$sapTPA/sapscale)>(F6$TPA+F6$errorTPA))|((F6$sapTPA/sapscale)<(F6$TPA-F6$errorTPA)), F6$sapTPA, NA)
#F6$x5grow <- ifelse(((F6$x5grow-F6$errorx5grow)>(F6$x5))|((F6$x5grow)<(F6$x5-F6$errorx5)), F6$x5grow, NA)

#set up 5-year precip plots using ggplot
F7 <- subset(F6, v==1)
F7$sapTPA <- ifelse(F7$col3==2,NA, F7$sapTPA)
v=1
col4 <- c(dt,"grey35", pt)

#NPP plot
#limitssap <- aes(ymax=(y5/y5scale +errorsapTPA/sapscale)+TPA3start1, ymin=(sapTPA/sapscale-errorsapTPA/sapscale)+TPA3start1)
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
P1 <- ggplot(F7, aes(x=x3, y=(y3+y3start1), fill=as.factor(col3)))  +geom_hline(y=y3start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21)  + theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ theme(axis.text.x = element_blank())+ scale_y_continuous(name=expression(paste(" Productivity (Mg C ha"^"-1","yr" ^"-1",")")), limits=c(Cymin, Cymax)+y3start1, breaks=c(round(y3start1, 1)-0.5,  round(y3start1, 3),  round(y3start1, 1)+0.3,  round(y3start1, 1)+0.6), labels=c(round(y3start1, 1)-0.5, round(y3start1, 1), round(y3start1, 1)+0.3, round(y3start1, 1)+0.6)) + scale_x_continuous(element_blank(), breaks=c(-0.5, 0, 0.5)) + theme(legend.title=element_blank(), legend.key=element_blank()) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=element_blank(), values=c(col4))  + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) #+coord_cartesian(ylim = c(0.7, 2.65))

 
#Stem density plot
limits <- aes(ymax=(TPA+errorTPA)+TPA3start1, ymin=(TPA-errorTPA)+TPA3start1)
limitssap <- aes(ymax=(sapTPA/sapscale +errorsapTPA/sapscale)+TPA3start1, ymin=(sapTPA/sapscale-errorsapTPA/sapscale)+TPA3start1)
P2 <- ggplot(F7, aes(x=x3, y=(TPA+TPA3start1), fill=as.factor(col3)))  +geom_hline(y=TPA3start1, linetype=3) +geom_vline(linetype=3) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(-0.07, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ theme(axis.text.x = element_blank())+ scale_y_continuous(name=expression(paste(Delta, " Stem density (ha"^"-1","yr" ^"-1",")")), limits=c(TPAymin, TPAymax)+TPA3start1, breaks=c(2.0, TPA3start1, 8.0), labels=c(expression("2.0"[" (-35)"]), expression("4.6" [" (43)"]), expression("8.0" [" (145)"]))) + scale_x_continuous(element_blank(), breaks=c(-0.5, 0, 0.5)) +scale_fill_manual(values=c(col4), breaks=NULL) +geom_errorbar(limitssap, width=0.05, colour="grey70") + geom_point(aes( y=(sapTPA/sapscale)+ TPA3start1, col=as.factor(col3)), size=2, pch=8)+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.1),plot.margin=unit(c(0,0,0.65,0), "cm")) +scale_colour_manual(values=c(col9), name=NULL, breaks=c("2"), label=c("(Saplings)"))+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21) #+coord_cartesian(ylim = c(1.8, 9.85))


#Radial growth plot
limits <- aes(ymax=(x5+errorx5)+x5start1, ymin=(x5-errorx5)+x5start1)
limits2 <- aes(ymax=(x5grow+errorx5grow)+ x5start1, ymin=(x5grow-errorx5grow)+ x5start1)
P3 <- ggplot(F7, aes(x=x3, y=(x5+x5start1), fill=as.factor(col3)))  +geom_hline(y=x5start1, linetype=3) +geom_vline(linetype=3)  +theme(axis.line = element_line( size = 0.35)) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm"))  + theme(axis.ticks.margin = unit(-0.03, "cm")) + theme(axis.ticks.length = unit(0.1, "cm"))  +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(name=expression(paste(Delta, " Mean diameter (cm yr" ^"-1",")")), limits=c(DIAmin, DIAmax)+x5start1, breaks=c(0.07, x5start1, 0.17), labels=c(expression("0.07"[" (0.55)"]), expression("0.11" [" (0.63)"]), expression("0.17" [" (0.75)"]))) + scale_x_continuous(name=expression(paste("   Precipitation JJA 5-yr\nmean anomaly (mm day"^"-1",")")), breaks=c(-0.5, 0, 0.5)) +scale_fill_manual(values=c(col4), breaks=NULL)  + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.55,0.1), plot.margin=unit(c(-0.5,0,0.34, -0.04), "cm")) + theme(axis.title.x = element_text(vjust=-0.85)) + geom_errorbar(limits2, width=0.05, colour="grey70") + geom_point(aes(y=(x5grow+ x5start1), col=as.factor(col3)), size=2.5, pch=17)  +scale_colour_manual(values=c(col4), name=NULL, breaks=c("2"), label=expression(paste("(Radial growth mm yr" ^"-1",")"))) + geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21)

  
#set up 5-year TEMP plots using ggplot
F8 <- subset(F6, v==2) # & x3>0.2)
v=2
col8 <- c(ct,"grey35", wt)
#NPP plot
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
T1 <- ggplot(F8, aes(x=x3, y=(y3+y3start1), fill=as.factor(col3)))  +geom_hline(y=y3start1, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x =element_blank())+ scale_y_continuous(element_blank(), limits=c(Cymin, Cymax)+y3start1,  breaks=c(round(y3start1, 1)-0.5,  round(y3start1, 3),  round(y3start1, 1)+0.3,  round(y3start1, 1)+0.6)) + scale_x_continuous(element_blank(),breaks=c(-0.5, 0, 0.5, 1.0))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Cool 5-yr period", "Warm 5-yr period"), values=c(col8)) +theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) #+coord_cartesian(ylim = c(0.7, 2.65))

#stem density
limits <- aes(ymax=(TPA+errorTPA)+TPA3start1, ymin=(TPA-errorTPA)+TPA3start1)
T2 <- ggplot(F8,aes(x=x3, y=(TPA+TPA3start1), fill=as.factor(col3)))    +geom_hline(y=TPA3start1, linetype=3) +geom_vline(linetype=3) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm"))  + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_blank())+ theme(axis.text.x =(element_blank())) + scale_y_continuous(element_blank(), limits=c(TPAymin, TPAymax)+TPA3start1, breaks=c(2.0, TPA3start1, 8.0)) + scale_x_continuous(element_blank())  +scale_fill_manual(values=c(col8), breaks=NULL) +geom_errorbar(limitssap, width=0.05, colour="grey70") + geom_point(aes( y=(sapTPA/sapscale)+ TPA3start1, col=as.factor(col3)), size=2, pch=8) + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.1),plot.margin=unit(c(0,0,0.5,0), "cm")) +scale_colour_manual(values=c(col8), name=NULL, breaks=c("2"), label=c("(Saplings)"))+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21) #+coord_cartesian(ylim = c(1.8, 9.85)) 

#Radial growth plot
limits <- aes(ymax=(x5+errorx5)+x5start1, ymin=(x5-errorx5)+x5start1)
limits2 <- aes(ymax=(x5grow+errorx5grow)+ x5start1, ymin=(x5grow-errorx5grow)+ x5start1)
T3 <- ggplot(F8, aes(x=x3, y=(x5+x5start1), fill=as.factor(col3))) + geom_hline(y=x5start1, linetype=3) +geom_vline(linetype=3)  +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y =element_blank())+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(element_blank(), limits=c(DIAmin, DIAmax)+x5start1,breaks=c(0.02, x5start1, 0.2)) + scale_x_continuous(name=("Temperature JJA 5-yr\nmean anomaly (C)")) +scale_fill_manual(values=c(col8), breaks=NULL) + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.43,0.1), plot.margin=unit(c(-0.5,0,0, 0), "cm"))+ geom_errorbar(limits2, width=0.05, colour="grey70") + geom_point(aes( y=(x5grow+ x5start1), col=as.factor(col3)), size=2.5, pch=17)  +scale_colour_manual(values=c(col8), name=NULL, breaks=c("2"), label=expression(paste("(Radial growth mm yr" ^"-1",")"))) + geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21) 
 



#set up 1-year soil moist. plots using ggplot
F9 <- subset(F6, v==3)
v=3
col9 <- c(de,"grey35", pe)
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
#NPP plot
R1 <- ggplot(F9,aes(x=x3, y=(y3+y3start1), fill=as.factor(col3))) +geom_hline(y=y3start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x =(element_blank())) + scale_y_continuous(element_blank(), limits=c(Cymin, Cymax)+y3start1,  breaks=c(round(y3start1, 1)-0.5,  round(y3start1, 3),  round(y3start1, 1)+0.3,  round(y3start1, 1)+0.6)) + scale_x_continuous(element_blank()) + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c( "Dry summer", "Wet summer"), values=c(col9)) +theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) #+coord_cartesian(ylim = c(0.7, 2.65))

#Stem density plot
limits <- aes(ymax=(TPA+errorTPA)+TPA3start1, ymin=(TPA-errorTPA)+TPA3start1)
limitssap <- aes(ymax=(sapTPA/sapscale +errorsapTPA/sapscale)+TPA3start1, ymin=(sapTPA/sapscale-errorsapTPA/sapscale)+TPA3start1)
R2 <- ggplot(F9,aes(x=x3, y=(TPA+TPA3start1), fill=as.factor(col3)))   +geom_hline(y=TPA3start1, linetype=3) +geom_vline(linetype=3)+theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x =(element_blank()))+ scale_y_continuous(element_blank(), limits=c(TPAymin, TPAymax)+TPA3start1, breaks=c(2.0, TPA3start1, 8.0)) + scale_x_continuous(element_blank()) +scale_fill_manual(values=c(col9), breaks=NULL) +geom_errorbar(limitssap, width=0.10, colour="grey70") + geom_point(aes( y=(sapTPA/sapscale)+ TPA3start1, col=as.factor(col3)), size=2, pch=8) + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.1), plot.margin=unit(c(0,0,0.5,0), "cm")) +scale_colour_manual(values=c(col9), name=NULL, breaks=c("2"), label=c("(Saplings)"))+ geom_errorbar(limits, width= 0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21) #+coord_cartesian(ylim = c(1.8, 9.85))

#radial growth plot
limits <- aes(ymax=(x5+errorx5)+x5start1, ymin=(x5-errorx5)+x5start1)
limits2 <- aes(ymax=(x5grow+errorx5grow)+ x5start1, ymin=(x5grow-errorx5grow)+ x5start1)
R3 <- ggplot(F9, aes(x=x3, y=(x5+x5start1), fill=as.factor(col3)))   +geom_hline(y=x5start1, linetype=3) +geom_vline(linetype=3) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_blank())+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(element_blank(), limits=c(DIAmin, DIAmax)+x5start1, breaks=c(0.02, x5start1, 0.2)) + scale_x_continuous(name=("Soil moisture JJA 1-yr\nmin/max (z-score)"), breaks=c(-2.0, -1.0, 0.0, 1.0), labels=c("-2.0", "-1.0", "0.0", "1.0")) +scale_fill_manual(values=c(col9), breaks=NULL) + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.53,0.1), plot.margin=unit(c(-0.5,0,0, 0), "cm")) + geom_errorbar(limits2, width=0.10, colour="grey70")+ geom_point(aes( y=(x5grow+ x5start1), col=as.factor(col3)), size=2.5, pch=17)  +scale_colour_manual(values=c(col9), name=NULL, breaks=c("2"), label=expression(paste("(Radial growth mm yr" ^"-1",")"))) + geom_errorbar(limits, width=0.01, colour="grey50")  + geom_point(colour="grey25", size=2.5, pch=21)

#set up plot window and plot all together
dev.new(width=6.75, height=8) #8
grid.arrange(P1,T1,R1,P2,T2,R2,P3,T3,R3, ncol=3)
#grid.text(0.045, unit(0.575,"npc") - unit(1,"line"), rot=90, label="(30)", gp=gpar(cex=0.75))
print(F7)
#######END Fig 3.  Climate Correlations######################################################









####Fig. S1. undisturbed Carbon climate plots ########################
#NOTE: run above analysis excuding disturbance

#calc plot y-scales
Cymax  <-  max(F6$y5+F6$errory5)
Cymin  <-  min(F6$y5-F6$errory5)

F7 <- subset(F6, v==1)
F7$sapTPA <- ifelse(F7$col3==2,NA, F7$sapTPA)
v=1
y5start1 <-mean(F6$y5start)
TPA3start1 <-mean(F6$TPA3start) #mean(F2$TPAchg) #
sapTPA3start1 <-mean(F6$sapTPA3start)
x5start1 <- mean((F6$x5start)*10)
col4 <- c(dt,"grey35", pt)
#NPP plot
limitssap <- aes(ymax=(y5/y5scale +errorsapTPA/sapscale)+TPA3start1, ymin=(sapTPA/sapscale-errorsapTPA/sapscale)+TPA3start1)
limits <- aes(ymax=(y5+errory5)+y5start1, ymin=(y5-errory5+y5start1))
P1 <- ggplot(F7, aes(x=x3, y=(y5+y5start1), fill=as.factor(col3)))  +geom_hline(y=y5start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21)  + theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(name=expression(paste("Productivity (Mg C ha"^"-1","yr" ^"-1",")")), limits=c(Cymin, Cymax)+y5start1, breaks=c(1.1, round(y5start1,2), 1.8 )) + scale_x_continuous(name=("Precipitation JJA 5-yr\nmean anomaly (mm/day)"), breaks=c(-0.5, 0, 0.5))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Wet 5-yr period", "Dry 5-yr period"), values=c(col4)) 

#set up 5-year TEMP plots using ggplot
F8 <- subset(F6, v==2) # & x3>0.2)
v=2
col8 <- c(ct,"grey35", wt)
#NPP plot
limits <- aes(ymax=(y5+errory5)+y5start1, ymin=(y5-errory5+y5start1))
T1 <- ggplot(F8, aes(x=x3, y=(y5+y5start1), fill=as.factor(col3)))  +geom_hline(y=y5start1, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(element_blank(), limits=c(Cymin, Cymax)+y5start1, breaks=c(1.1, round(y5start1,2), 1.8 )) + scale_x_continuous(name=("Temperature JJA 5-yr\nmean anomaly (C)"),breaks=c(-0.5, 0, 0.5, 1.0))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Cool 5-yr period", "Warm 5-yr period"), values=c(col8))

#set up 1-year soil moist. plots using ggplot
F9 <- subset(F6, v==3)
v=3
col9 <- c(de,"grey35", pe)
limits <- aes(ymax=(y5+errory5)+y5start1, ymin=(y5-errory5+y5start1))
#NPP plot
R1 <- ggplot(F9,aes(x=x3, y=(y5+y5start1), fill=as.factor(col3))) +geom_hline(y=y5start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.02, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(element_blank(), limits=c(Cymin, Cymax)+y5start1, breaks=c(1.1, round(y5start1,2), 1.8 ))  + scale_x_continuous(name=("Soil moisture JJA 1-yr\nmin/max (z-score)"))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Wet summer", "Dry summer"), values=c(col9))

dev.new(width=7, height=3)
grid.arrange(P1,T1,R1, ncol=3)




####Fig. S1. Deciduous/Coniferous Carbon climate plots ########################
#NOTE: run above analysis excuding Coniferous or Deciduous plots (i.e. SISP>=300 is only deciduous plots)


#set up 5-year precip plots using ggplot
F7 <- subset(F6, v==1)
F7$sapTPA <- ifelse(F7$col3==2,NA, F7$sapTPA)
v=1
F7$col3[2]<-3
F7$col3[7]<-1
col4 <- c(dt,"grey35", pt)
#NPP plot
#limitssap <- aes(ymax=(y5/y5scale +errorsapTPA/sapscale)+TPA3start1, ymin=(sapTPA/sapscale-errorsapTPA/sapscale)+TPA3start1)
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
limitssap <- aes(ymax=(y3/y3scale +errorsapTPA/sapscale)+TPA3start1, ymin=(sapTPA/sapscale-errorsapTPA/sapscale)+TPA3start1)
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
P1 <- ggplot(F7, aes(x=x3, y=(y3+y3start1), fill=as.factor(col3)))  +geom_hline(y=y3start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21)  + theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(name=expression(paste("Productivity (Mg C ha"^"-1","yr" ^"-1",")")), limits=c(Cymin, Cymax)+y3start1, breaks=c( round(y3start1,2)-0.5, round(y3start1,2),  round(y3start1,2)+0.5 )) + scale_x_continuous(name=("Precipitation JJA 5-yr\nmean anomaly (mm/day)"), breaks=c(-0.5, 0, 0.5))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Wet 5-yr period", "Dry 5-yr period"), values=c(col4)) 

#set up 5-year TEMP plots using ggplot
F8 <- subset(F6, v==2) #
v=2
F8$col3[2]<-3
F8$col3[7]<-1
col8 <- c(ct,"grey35", wt)
#NPP plot
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
T1 <- ggplot(F8, aes(x=x3, y=(y3+y3start1), fill=as.factor(col3)))  +geom_hline(y=y3start1, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey30")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(element_blank(), limits=c(Cymin, Cymax)+y3start1, breaks=c(1.1, round(y3start1,2), 1.8 )) + scale_x_continuous(name=("Temperature JJA 5-yr\nmean anomaly (C)"),breaks=c(-0.5, 0, 0.5, 1.0))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.52,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Cool 5-yr period", "Warm 5-yr period"), values=c(col8))

#set up 1-year soil moist. plots using ggplot
F9 <- subset(F6, v==3)
v=3
F9$col3[2]<-3
F9$col3[7]<-1
col9 <- c(de,"grey35", pe)
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
#NPP plot
R1 <- ggplot(F9,aes(x=x3, y=(y3+y3start1), fill=as.factor(col3))) +geom_hline(y=y3start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.02, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(element_blank(), limits=c(Cymin, Cymax)+y3start1, breaks=c(1.1, round(y3start1,2), 1.8 ))  + scale_x_continuous(name=("Soil moisture JJA 1-yr\nmin/max (z-score)"))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Wet summer", "Dry summer"), values=c(col9))

dev.new(width=7, height=3)
grid.arrange(P1,T1,R1, ncol=3)

