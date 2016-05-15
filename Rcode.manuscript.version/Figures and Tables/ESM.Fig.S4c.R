#Nature CC resubmit Extended Data Figs

#Note other smaller spatial subsets do not have enough data and range of temperature anomalies to accurately identify response. 

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

#RM5mergedfullcut1991.11.17.2014  <- as.data.table(read.csv("/Users/travis/Desktop/Null2015full.csv", header = TRUE, sep = ",", quote="\"", dec="."))

#RM5mergedfullcut1991.11.17.2014 <-  subset(RM5mergedfullcut1991.11.17.2014) # remove planted forests
S4swPREC <-  subset(RM5mergedfullcut1991.11.17.2014,  STDORGCD==0 & cut==0)   #exclude planted and harvested forests
F2 <- subset(S4swPREC)


#################################################################################################################
################ FIGURE 2 - Climate Condition DIVERGENCE - CARBON and TEMPERATURE at Sub REGIONAL ########################################################
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

F2 <- subset(S4swPREC) #, temp5grow>0  ) 


#Convert variables to standard names and metric units
F2$DIAmean <- (F2$DIAendmean* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean <- (F2$DIAbeginmean *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum <- (F2$endsumTPA /0.404686)			#number of trees per hectare by plot at resurvey
F2$PREV_TPAsum <- (F2$startsumTPA /0.404686)	#number of trees per hectare by plot at initial survey (cm)
F2$PREV_STOCKINGmid <- F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey
F2$TPAchg <- (F2$TPAsum-F2$PREV_TPAsum)/F2$REMPER
F2$SISPchg<-ifelse(F2$SISP>=300, 1, ifelse(F2$SISP<300, -1, NA) )
#F2$PREVcarbon1sum<-F2$PREVcarbon5sum #F2$PREVNPP
#F2$carbon1sum<-F2$carbon5sum #F2$NPP


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


#F2$PREC5grow <- F2$PREC5grow*92/3/10 #convert mm/day to cm per month
F2$PRECminLL <- ave(F2$PREC5, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$PRECgrowminLL <- ave(F2$PREC5grow, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$PREC5chg <- F2$PREC5-F2$PRECminLL
F2$PREC5growchg <- F2$PREC5grow-F2$PRECgrowminLL
F2$PREC5growchg[F2$PREC5growchg==0] <- rnorm(sum(F2$PREC5growchg==0),0,sd=0.01) #make median 0 value a distribution for subseting

F2$SWdroughtLL <- ave(F2$SWdrought5, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$SWpluvialLL <- ave(F2$SWpluvial5, F2$LATLON, FUN=function(x) median(x, na.rm=TRUE))
F2$SWdrought5chg <- F2$SWdrought5-F2$SWdroughtLL
F2$SWpluvial5chg <- F2$SWpluvial5-F2$SWpluvialLL
F2$SWpluvial5chg[F2$SWpluvial5chg==0] <- rnorm(sum(F2$SWpluvial5chg==0),0,sd=0.01) #make median 0 value a distribution for subseting
F2$SWdrought5chg[F2$SWdrought5chg==0] <- rnorm(sum(F2$SWdrought5chg==0),0,sd=0.01) #make median 0 value a distribution for subseting
F2$nsampTOT <- ave(F2$DIAmean>0,F2$LATLON, FUN=function(x) sum(x, na.rm=TRUE))

F2$allmortality <- (F2$insect+F2$fire+F2$weather+F2$disease+F2$animal+F2$unknowndamage)
F2$alldisturb <- (F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)

#qt1 <- c(seq(0,1,1/10)/100,seq(1,59,2)/100)
#qt2 <- c(seq(41,99, 2)/100, seq(99,100,1/10)/100) #0.98, 0.99 added to ease moving average transistion to 1 at end


qt1 <- c(-2, -2, -2, -0.2, 0)
qt2 <- c(0, 0.2, 2, 2, 2) #c(seq(0, 1.25, 1.25/4)) #0.98, 0.99 added to ease moving average transistion to 1 at end
int <- 1


F5 <- NULL


for(l in 1:6)
{
	
for(v in 2:2)
{v=v
F4<-foreach(c = 1:10, .combine=rbind) %do%
{
	
	# #loop for each spatial subsection

	 if(l==1){ F2latlon <- subset(F2, LAT >= 45 & LON < -85) }
	 if(l==2){ F2latlon <- subset(F2, LAT >= 45 & LON >= -85)  }
	 if(l==3){ F2latlon <- subset(F2, LAT > 35 & LAT < 45 & LON <= -85) } 
	 if(l==4){ F2latlon <- subset(F2, LAT > 35 & LAT < 45 & LON > -85) } 
	 if(l==5){ F2latlon <- subset(F2, LAT <= 35 & LON <= -85) }
	 if(l==6){ F2latlon <- subset(F2, LAT <= 35 & LON > -85) }


vari <- c("PREC5growchg", "TEMP5growchg", "SWpluvial5chg", "SWdrought5chg") #"PREC5growchg", "TEMP5growchg"
vari1 <- subset(F2latlon, select=c(PREC5growchg, TEMP5growchg, SWpluvial5chg, SWdrought5chg))
vari2 <- c("5-yr JJA precip. anomaly (mm/day)", "5-yr JJA temp. anomaly (C)", "1-yr JJA soil moist. anomaly (z-score)", "1-yr JJA soil moist. anomaly (z-score)")
	

	if(c<=5 & v<3)
	{
		f22 <- subset(F2latlon, get(vari[v])>  qt2[c] & get(vari[v])<  qt2[c+int]) 		
		f23 <- subset(F2latlon, get(vari[v])>= quantile(vari1[[v]], 0.3, na.rm=TRUE) & get(vari[v])<= quantile(vari1[[v]], 0.7, na.rm=TRUE), select=c(PLT_CN)) 
		ccond <- vari[v]
		col3=ifelse(qt2[c]>= 0.2, 3, 2)
		qta <- qt2[c]
		qtb <- qt2[c+int]
		}
	if(c>5 & v<3)
	{ 
		d=c-5
		f22 <- subset(F2latlon, get(vari[v])< qt1[d] & get(vari[v])> qt1[d-int]) #		
		f23 <- subset(F2latlon, get(vari[v])>= quantile(vari1[[v]], 0.3, na.rm=TRUE) & get(vari[v])<= quantile(vari1[[v]], 0.7, na.rm=TRUE), select=c(PLT_CN)) 
		ccond <- vari[v]
		col3=ifelse(qt1[d]<= -0.2, 1, 2)
		qtb <- qt1[d]
		qta <- qt1[d-int]
	}

#set up spatial weighting of comparison model
f22$nsamp <- (ave(f22$DIAmean>0,f22$LATLON, FUN=function(x) sum(x, na.rm=TRUE)))/f22$nsampTOT #sample proportion of each grid cell 
#print(c(table(f22$nsamp<0.75), col3))

#f22$col3 <- col3
PLT_CNsamp1 <- f23$PLT_CN   

f24 <- subset(f22, !duplicated(LATLON), select=c(LATLON, nsamp))

F33 <- subset(F2latlon, LATLON %in% (f22$LATLON)) # not neccesaary bc weighting will eliminate latlons not in f22!
F33 <- merge(F33, f24, all.x=T, by="LATLON")
F33 <- subset(F33, PLT_CN %in% (PLT_CNsamp1)) # not neccesaary bc weighting will eliminate latlons not in f22!

F33$aPREVSTOCKDIAbin <- floor(F33$PREV_STOCKINGmid/10)*10000+(floor(F33$PREV_DIAmean/2)*2)
F33$n1 <- 1
F33$aPREVSTOCKDIAbinlength <- ave(F33$n1, F33$aPREVSTOCKDIAbin, FUN=function(x) sum(x, na.rm=TRUE))
F33$nbin <- ave(F33$nsamp, F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE))

F33 <- subset(F33, aPREVSTOCKDIAbinlength>= 1) #Null model needs at least 10 comparison plots at each diameter/stocking bin <- changing should not impact robust models


##Carbon Model
c6 <- (((F33$carbon1sum-F33$PREVcarbon1sum)* 0.00112085116 )/F33$REMPER)
F33$avgCchg <- (ave(c6*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE)))/F33$nbin

#Merge null and conditional model data
F34 <- subset(F33, !duplicated(aPREVSTOCKDIAbin), select=c(aPREVSTOCKDIAbin, avgCchg)) #, avgx5, avgTPAchg, sapavgTPAchg, avgC5chg, avgmortx5, sapavgMORT, sapavgRECR))
f22$aPREVSTOCKDIAbin <- floor(f22$PREV_STOCKINGmid/10)*10000+(floor(f22$PREV_DIAmean/2)*2)
f22 <- merge(f22, F34, by="aPREVSTOCKDIAbin") #, all.x=T

#carbon 1 inch
c8 <- (((f22$carbon1sum-f22$PREVcarbon1sum)* 0.00112085116 )/f22$REMPER)
f22$xcarbon <- ((c8)-f22$avgCchg)


#calculate mean condition for x-axis based on v loop
if(length(f22$xcarbon)>= 2)
{
if(v==1)
{
	x3 <- mean(f22$PREC5growchg, na.rm=TRUE)
}
if(v==2)
{
	x3 <- mean(f22$TEMP5growchg, na.rm=TRUE)
}


#calculate y-axis variables
y3 <- mean(f22$xcarbon, na.rm=TRUE) 
y3pctdif<-(mean(f22$xcarbon, na.rm=TRUE)/mean(f22$avgCchg, na.rm=TRUE))*100   #

y3start <- NA
if(col3==2){y3start <-  mean(c6, na.rm=T)} #mean change in forest npp during median climate conditions


#calc y-axis error
R2 <- 50 #repetitions set to 10 for faster looping
enderror <- boot(f22$xcarbon, meanFunc, R=R2)
enderrory3 <- sqrt(var(enderror$t))*1.96

n <- length(f22$xcarbon) #conditional model n
nnull <- length(F33$avgCchg) #null model n
y3null<- mean(f22$avgCchg, na.rm=T)

print(c(x3, n, nnull))

f23 <- c(x3, y3, enderrory3, n, v, l, col3, y3start, y3pctdif, y3null)

}
}
F5 <- rbind(F5,F4)
}
}


F6 <- as.data.table(F5)
setnames(F6, c("x3","y3", "errory3", "n", "v", "l", "col3", "y3start", "y3pctdif", "y3null"))
F6 <- subset(F6, n>=2) 


#set up plots using ggplot
F81 <- subset(F6, v==2 & l==1)
v=1
y3start1 <-mean(F81$y3start, na.rm=T)
col8 <- c(ct,"grey35", wt)
Cymax1  <- max(F81$y3+ F81$errory3)
Cymin1  <- min(F81$y3-F81$errory3)
Cmid1 <- round(y3start1, 2)
Cmin1<- round(Cmid1-Cmid1*0.2, 2)
Cmax1<- round(Cmid1+Cmid1*0.2, 2)

#NPP plot
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
L1 <- ggplot(F81, aes(x=x3, y=(y3+y3start1), fill=as.factor(col3)))  +geom_hline(y=y3start1, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + scale_y_continuous(name=(expression(atop(paste("North of 45"^"o", " N"), paste("Productivity (Mg C ha"^"-1","yr" ^"-1",")")))), limits=c(Cymin1-0.2, Cymax1)+y3start1, breaks=c(Cmin1, Cmid1, Cmax1)) + scale_x_continuous(element_blank(), limits=c(-1.0, 1.0),breaks=c(-0.5, 0, 0.5, 1.0))+ scale_fill_manual(values=c(col8))+ ggtitle(expression(paste("-95"^"o"," W to -85"^"o"," W" ))) + theme(plot.title=element_text( size=12))


F82 <- subset(F6, v==2 & l==2)
v=1
y3start2 <-mean(F82$y3start, na.rm=T)
col8 <- c(ct,"grey35", wt)
Cymax2  <- max(F82$y3+ F82$errory3)
Cymin2  <- min(F82$y3-F82$errory3)
Cmid2 <- round(y3start2, 2)
Cmin2<- round(Cmid2-Cmid2*0.1, 2)
Cmax2<- round(Cmid2+Cmid2*0.1, 2)

#NPP plot
limits <- aes(ymax=(y3+errory3)+y3start2, ymin=(y3-errory3+y3start2))
L2 <- ggplot(F82, aes(x=x3, y=(y3+y3start2), fill=as.factor(col3)))  +geom_hline(y=y3start2, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + scale_y_continuous(element_blank(), limits=c(Cymin2-0.2, Cymax2+0.15)+y3start2, breaks=c(Cmin2, Cmid2, Cmax2)) + scale_x_continuous(element_blank(), limits=c(-1.0, 1.0),breaks=c(-0.5, 0, 0.5, 1.0))  + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.75,0.9),plot.margin=unit(c(0,0,0,0), "cm"))  +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Cool 5-yr period", "Warm 5-yr period"), values=c(col8)) + ggtitle(expression(paste("East of -85"^"o"," W" ))) + theme(plot.title=element_text( size=12))


F83 <- subset(F6, v==2 & l==3)
v=1
y3start3 <-mean(F83$y3start, na.rm=T)
col8 <- c(ct,"grey35", wt)
Cymax3  <- max(F83$y3+ F83$errory3)
Cymin3  <- min(F83$y3-F83$errory3)
Cmid3 <- round(y3start3, 2)
Cmin3 <- round(Cmid3-Cmid3*0.2, 2)
Cmax3 <- round(Cmid3+Cmid3*0.2, 2)

#NPP plot
limits <- aes(ymax=(y3+errory3)+y3start3, ymin=(y3-errory3+y3start3))
L3 <- ggplot(F83, aes(x=x3, y=(y3+y3start3), fill=as.factor(col3)))  +geom_hline(y=y3start3, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + scale_y_continuous(name= (expression(atop(paste("35"^"o", " N to 45"^"o", " N"), paste("Productivity (Mg C ha"^"-1","yr" ^"-1",")")))), limits=c(Cymin3, Cymax3)+y3start3, breaks=c(Cmin3, Cmid3, Cmax3)) + scale_x_continuous(element_blank(), limits=c(-1.0, 1.0),breaks=c(-0.5, 0, 0.5, 1.0))+scale_fill_manual(values=c(col8)) 

F84 <- subset(F6, v==2 & l==4)
v=1
y3start4 <-mean(F84$y3start, na.rm=T)
col8 <- c(ct,"grey35", wt)
Cymax4  <- max(F84$y3+ F84$errory3)
Cymin4  <- min(F84$y3-F84$errory3)
Cmid4 <- round(y3start4, 2)
Cmin4 <- round(Cmid4-Cmid4*0.2, 2)
Cmax4 <- round(Cmid4+Cmid4*0.2, 2)

#NPP plot
limits <- aes(ymax=(y3+errory3)+y3start4, ymin=(y3-errory3+y3start4))
L4 <- ggplot(F84, aes(x=x3, y=(y3+y3start4), fill=as.factor(col3)))  +geom_hline(y=y3start4, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + scale_y_continuous(element_blank(), limits=c(Cymin4, Cymax4)+y3start4, breaks=c(Cmin4, Cmid4, Cmax4)) + scale_x_continuous(element_blank(), limits=c(-1.0, 1.0),breaks=c(-0.5, 0, 0.5, 1.0)) +scale_fill_manual(values=c(col8)) 

F85 <- subset(F6, v==2 & l==5)
v=1
y3start5 <-mean(F85$y3start, na.rm=T)
col8 <- c(ct,"grey35", wt)
Cymax5  <- max(F85$y3+ F85$errory3)
Cymin5  <- min(F85$y3-F85$errory3)
Cmid5 <- round(y3start5, 2)
Cmin5 <- round(Cmid5-Cmid5*0.2, 2)
Cmax5 <- round(Cmid5+Cmid5*0.2, 2)

#NPP plot
limits <- aes(ymax=(y3+errory3)+y3start5, ymin=(y3-errory3+y3start5))
L5 <- ggplot(F85, aes(x=x3, y=(y3+y3start5), fill=as.factor(col3)))  +geom_hline(y=y3start5, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + scale_y_continuous(name=(expression(atop(paste("South of 35"^"o", " N"), paste("Productivity (Mg C ha"^"-1","yr" ^"-1",")")))), limits=c(Cymin5, Cymax5)+y3start5, breaks=c(Cmin5, Cmid5, Cmax5)) + scale_x_continuous(name=("Temperature JJA 5-yr\nmean anomaly (C)"), limits=c(-1.0, 1.0),breaks=c(-0.5, 0, 0.5, 1.0))+scale_fill_manual(values=c(col8)) 

F86 <- subset(F6, v==2 & l==6)
v=1
y3start6 <-mean(F86$y3start, na.rm=T)
col8 <- c(ct,"grey35", wt)
Cymax6  <- max(F86$y3+ F86$errory3)
Cymin6  <- min(F86$y3-F86$errory3)
Cmid6 <- round(y3start6, 2)
Cmin6 <- round(Cmid6-Cmid6*0.2, 2)
Cmax6 <- round(Cmid6+Cmid6*0.2, 2)

#NPP plot
limits <- aes(ymax=(y3+errory3)+y3start6, ymin=(y3-errory3+y3start6))
L6 <- ggplot(F86, aes(x=x3, y=(y3+y3start6), fill=as.factor(col3)))  +geom_hline(y=y3start6, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + scale_y_continuous(element_blank(), limits=c(Cymin6, Cymax6)+y3start6, breaks=c(Cmin6, Cmid6, Cmax6)) + scale_x_continuous(name=("Temperature JJA 5-yr\nmean anomaly (C)"), limits=c(-1.0, 1.0),breaks=c(-0.5, 0, 0.5, 1.0))+scale_fill_manual( values=c(col8)) 


#set up plot window and plot all together
dev.new(width=6, height=7)
grid.arrange(L1, L2, L3, L4, L5, L6, ncol=2)

