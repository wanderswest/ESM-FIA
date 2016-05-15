#EMPIRICAL SUCCESSION MAPPING Feb. 2015
# Vegetation carbon forward projections
#TRAVIS ANDREWS - TDA210@lehigh.edu

####CAUTION for() LOOPS AHEAD####

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
 # data.url <- "https://raw.githubusercontent.com/wanderswest/ESM-FIA/master/ESM.data.csv"
 
    # ESM.data <- getURL(data.url)                
    # ESM.data <- read.csv(textConnection(ESM.data), header = TRUE, sep = ",", quote="\"", dec=".")
# RM5mergedfullcut1991.11.17.2014  <- as.data.table(ESM.data)

#or load filtered data from local drive
#RM5mergedfullcut1991.11.17.2014  <- as.data.table(read.csv("/Users/travis/GitHub/ESM-FIA/ESM.data.csv", header = TRUE, sep = ",", quote="\"", dec="."))
RM5mergedfullcut1991.11.17.2014  <- as.data.table(read.csv("/Users/travis/Desktop/ESMthin/Null2015PFTs.csv", header = TRUE, sep = ",", quote="\"", dec="."))

RM5mergedfullcut1991.11.17.2014 <-  subset(RM5mergedfullcut1991.11.17.2014, STDORGCD==0 ) # remove planted forests
S4swPREC <-  subset(RM5mergedfullcut1991.11.17.2014, STDORGCD==0 & cut==0 ) #exclude planted and harvested forests
F2 <- subset(S4swPREC)

print(c(sum(RM5mergedfullcut1991.11.17.2014$endsumTPA)/(1/0.166166884), "Trees resurveyed"))


##########
#######################   SUPP. FIGURE -CARBON GROWTH MODEL  #######################################################
######
####
###
#subset resurveyed plots to those with at least one 5" tree in order to eliminate plots only actually surveyed once.
F2<-subset(S4swPREC) #

#carbon
F2$DIAmean<-(F2$DIAendmean* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean<-(F2$DIAbeginmean *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum<-F2$carbon1sum* 0.00112085116  #	#convert living carbon to Mg/ha
F2$PREV_TPAsum<-F2$PREVcarbon1sum* 0.00112085116 ##convert living carbon to Mg/ha
F2$PREV_STOCKINGmid<-F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey

#calculate change in Stem density and mean tree diameter for each plot
F2$TPAchg<-F2$TPAsum-F2$PREV_TPAsum
F2$DIAchg<-F2$DIAmean-F2$PREV_DIAmean
#######################   Extended data FIGURE 4 - EMPIRICAL FOREST CARBON MODEL  #######################################################



#set up plot
dev.new(width=8, height=6)
plot(c((F2$PREV_DIAmean),(F2$DIAmean)),c((F2$PREV_TPAsum),(F2$TPAsum)), type ="n", col="grey94", 
	main=c("Empirical Model of 5-year forest carbon change, DBH>=12.7 cm", "FIA survey 1997-2012, Eastern USA") 
	,ylim=c(60,165), xlim=c(23, 40), mgp=c(2,1,0)
	#,ylim=c(58,90), xlim=c(23, 29), mgp=c(2,1,0)
	,xlab=NA, ylab=NA, axes=FALSE) #1250

axis(side=1, tck=-0.01, labels=NA, lwd=0.75)
axis(side=2, tck=-0.01, labels=NA, lwd=0.75)
axis(side=1, lwd=0, line= -0.7)
axis(side=2, lwd=0, line= -0.7)

mtext(side=1, "Mean tree diameter (cm)", line=1.2)
mtext(side=2, expression(paste("Carbon, living trees (Mg ha"^"-1",")")), line=1.2)
box()

#plot all plots as individual background vectors colored by inital stocking 
arrows((F2$PREV_DIAmean),(F2$PREV_TPAsum),(F2$DIAmean),(F2$TPAsum), col= "grey95", length=0.065, angle=22,lwd=0.9)   #col= "grey85"

#Set up plot with relative plot locations and legend - ideal for large window
x1<-14-0.5
x2<-14+0.5
y1<-(c(seq(32,26.05,-.55))*2.5-52)*6.1
y2<-(c(seq(32.2,26.25,-0.085))*2.5-52)*6.1
y3<-(c(seq(31.8,26,-0.73))*2.5-52)*6.1
x1b<-x1-7 #x1b<-x1-7.5
y1b<-y1[1]-264 #y1b<-y1[1]-215
cex1<-0.8
t3<-32 #transition mean diameter between independent vectors and moving average
lw2<-2 #arrow width

stockcols2<-"grey75" #rev((seq(1,120,15)/4)+18)
segments(x1,y2,x2, y2, c(seq(45,18, -0.384)),lwd=2, lend=2) #c(seq(66,38,-0.3))  col="grey85"
points(30,100, pch=16, cex=500, col="#FFFFFF75" ) #FFFFFF75

text(x1-0.25, y1[1]+26,"5-year Temperate Forest Carbon Change", pos=4, cex=1.05, font=3)

#segments(33.5,0,33.5,600, col="white", lwd=0.75)
segments(c(x1,x1,x1),c( y2[1], y2[36],  y2[71]),x2+0.15,c(y2[1], y2[36],  y2[71]), col="black",lwd=0.9, lend=2)
segments(c(x1),c(  y2[10]),x2+0.15,c(y2[10]), col="black",lwd=0.9, lend=2, lty="dashed")
arrows(x1+0.1,y3,x2-0.1,y3, col="black", length=0.085, angle=25,lwd=lw2+1)  
arrows(x1+0.1,y3,x2-0.1,y3, col=c(stockcols2), length=0.08, angle=28,lwd=lw2)  
text(x1-.5, y1[1]+10,"Initial relative density", cex=0.9, pos=4) #initial plot stocking and mean change vector increments
#text(x1+0.83, c(y1, y1[11]),c("120","100","90", "80","70", "60","50","40","30", "20", "10", "1"), cex=0.9, pos=4)
text(x2, y1[1]-3,"Over", cex=cex1, pos=4) #initial plot stocking and mean change vector increments
text(x2, y1[5],"Full", cex=cex1, pos=4) #initial plot stocking and mean change vector increments
text(x2, y1[10],"Under", cex=cex1, pos=4) #initial plot stocking and mean change vector increments

#Carbon thinning line 
x=seq(0,50,0.5) #seq(22,35.5,0.5)
y=(6*x-34)  #Best fit for early succession
lines((y)~(x), col="#00000090", lty="dotted", lwd=5) # "#00000085"
#xb=seq(14.7,20.2,0.5)
#yb=(2.2*xb-30)  #Best fit for early succession
#lines((yb)~(xb), col="#00000090", lty="dotted", lwd=5) # "#00000085"
#arrows(x[5], y[18], x[8], y[8]+0.2, col="#00000090", length=0.085, angle=25,lwd=1)  
#text(x[5], y[18]-0.5,"y = x - 6", cex=cex1, pos=3) #initial plot stocking and mean change vector increments


#text(19.2, 1200, "Self thinning -3/2 power function \n           (Reineke 1933)" , pos=4, srt=-48, col="#00000098", cex=cex1)

#set up plot with data info - works with par$pin 8.5,4.7
#text(x1-0.5, y1[11]-5,"Eastern US Forest Interiors\n n =            resurveyed plots\n trees >= 12.7cm DBH\n FIA 1998-2012", pos=4, cex=0.85) # col="grey30", #min. tree diameter = 12.7cm\n
#text(x1+0.5, y1[11]-3.5, (length(F2$PREV_DIAmean)), pos=4,  cex=0.85)


####Organize and subset data to plot with FOR loops 
#looped by phase space location to have sufficient underlying data and to exclude excessive vectors 


####Organize and subset data to plot with FOR loops 
#looped by phase space location to have sufficient underlying data and to exclude excessive vectors
DuplicatedVectors <- c(5031, 5034.6, 5036.4, 15036.4, 40040, 65040, 85038.2, 65038.2, 85036.4, 65036.4 ) 
slopetable <- NULL
for(o in 1:2){
for(l in 1:6){

if(l==1){	  ##Reforestation
STOCKval <- c(4,15)
QMDvalues <- c(seq(21.7,30.7, 3), 36, 45) #c(18,21,24,29, 45)  #
s1 <- 5
s2 <- 5
s3 <- 5
s4 <- 5
increment=1}

if(l==2){	  ##Early succession		
STOCKval <- c(seq(0, 70, 121/7), 82, 101, 121) #6
QMDvalues <- c(seq(14.5,21.7,1.8)) #c(seq(12,18.4,3.2))  #
s1 <- 5 #30
s2 <- 5 #10
s3 <- 5 #10
s4 <- 5 #5
increment=1}

if(l==3){	  ##middle
STOCKval <- c(15, 32, seq(41, 85, 65/5),98,121) #121,121,121
QMDvalues <- c(seq(21.7,t3, 10.3/6))  #
s1 <- 30
s2 <- 10
s3 <- 20
s4 <- 2
increment=1}

if(l==4){	  ##moving average end
STOCKval <- c(5, 15, NA, 30, 38, NA, NA,20, 40, 55, 55, 60, 75, NA, NA, NA, 65, 85, NA, 98, 100) #121, 121) #121,121,121
QMDvalues <- c(seq(t3-1,45,1.8), 50)  #
s1 <- 40 #100
s2 <- 30 #50
s3 <- 25 #40
s4 <- 30 #5
increment=3}

if(l==5){	  ##Early early succession
STOCKval <- c(0, 5, 15, seq(30, 121, 91/6)) #6
QMDvalues <- c(12.7, 14.5) 
s1 <- 3 #30
s2 <- 3 #10
s3 <- 3#10
s4 <- 3#5
increment=1}

if(l==6){	  ##moving average end top
STOCKval <- c(100, NA, NA, 125) 
QMDvalues <- c(seq(t3-1,46,1.8), 50)  #
s1 <- 5 #100
s2 <- 5 #50
s3 <- 5# 3 40
s4 <- 5 #5
increment=3}

#Subset by stocking range then mean diameter range

#Replot with bold vectors

g <- length(QMDvalues-1)
r <- length(STOCKval-1)
 for(s in 1:r)
 { S11r <- subset(F2, PREV_STOCKINGmid>= STOCKval[s] & PREV_STOCKINGmid <STOCKval[s+increment]) 	
	for(h in 1:g)
	{ 
	S11q <- subset(S11r, PREV_DIAmean>=(QMDvalues[h]) & PREV_DIAmean<(QMDvalues[h+increment]))
#print(length(S11q$PREV_DIAmean))
	if(length(S11q$PREV_DIAmean)>s1 & s==1 | length(S11q$PREV_DIAmean)>s2 & s==2 | length(S11q$PREV_DIAmean)>s3  & s==3 | length(S11q$PREV_DIAmean)>s4  & s>=4)
		{	
		stkqmd <- STOCKval[s]*1000+QMDvalues[h]
		if(stkqmd %in% DuplicatedVectors){next} #remove overlapping vectors			
	arrows((mean(S11q$PREV_DIAmean, na.rm=TRUE)),(mean(S11q$PREV_TPAsum, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE)), col= "grey90", length=0.08, angle=28,lwd=3)  # VectCols[s] col=s*1.8+10,  s*1.2+7 # s*2+10 #stockcol
print(c(stkqmd, l, s, h, mean(S11q$TPAchg, na.rm=T), sum(S11q$TPAchg>0, na.rm=T)))

	}}}
 } #l loop
} #o loop
box()


############  END CARBON PLOT SET UP
#	
#	

############### Carbon Empirical GROWTH Model - Climate condtions and disturbance PROJECTIONS!! ########################################

meanFunc<-function(x,i){mean(x[i], na.rm=TRUE)}

#subset resurveyed plots to those with at least one 5" tree in order to eliminate plots only actually surveyed once.
F2<-subset(RM5mergedfullcut1991.11.17.2014, PREVcarbon1sum>0) # & STDORGCD==0 ) 

#carbonC
F2$DIAmean<-(F2$DIAendmean* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean<-(F2$DIAbeginmean *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum<-(F2$endsumTPA /0.404686)			#number of trees per hectare by plot at resurvey
F2$PREV_TPAsum<-(F2$startsumTPA /0.404686)	#number of trees per hectare by plot at initial survey (cm)
F2$Csum<-F2$carbon1sum* 0.00112085116  #	#convert living carbon to Mg/ha
F2$PREV_Csum<-F2$PREVcarbon1sum* 0.00112085116 ##convert living carbon to Mg/ha
F2$PREV_STOCKINGmid<-F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey
F2$STOCKINGmid<-F2$STOCKING5mid


#calculate change in Stem density and mean tree diameter for each plot
F2$Cchg<-F2$Csum-F2$PREV_Csum
F2$DIAchg<-F2$DIAmean-F2$PREV_DIAmean
F2$STOCKchg<-F2$STOCKINGmid-F2$PREV_STOCKINGmid
#set grid cell size
a<-b<-2  #2=0.5degrees, 0.5=2degrees
	F2$iLAT<-((ceiling(F2$LAT*a))/a)-0.5 #+0.25  #round down to integer then add .25 ; switched /* removed ).25 addition
	F2$iLON<-((floor(F2$LON*b))/b)+0.5 #-0.25  #round down to integer then subtract .25; switched /*  ).25 addition
	F2$LATLONYR<-(F2$iLAT*10000000000000) +(F2$iLON*(-100000000)) +F2$INVYR
	F2$LATLON<-as.factor(F2$iLAT*100000000+ F2$iLON*-1000) #LATLON identifyier

#Find climate variable mean and range for each grid cell
F2$TEMPminLL<-ave(F2$temp5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$TEMPgrowminLL<-ave(F2$temp5grow, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$TEMP5chg<-F2$temp5-F2$TEMPminLL
F2$TEMP5growchg<-F2$temp5grow-F2$TEMPgrowminLL

F2$PRECminLL<-ave(F2$PREC5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$PRECgrowminLL<-ave(F2$PREC5grow, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$PREC5chg<-F2$PREC5-F2$PRECminLL
F2$PREC5growchg<-F2$PREC5grow-F2$PRECgrowminLL

F2$SWdroughtLL<-ave(F2$SWdrought5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$SWpluvialLL<-ave(F2$SWpluvial5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$SWdrought5chg<-F2$SWdrought5-F2$SWdroughtLL
F2$SWpluvial5chg<-F2$SWpluvial5-F2$SWpluvialLL

F2$NOdisturb<-ifelse((F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)==0,1,0) 
F2$alldisturb<-ifelse((F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)>0 & (F2$unknowndamage+F2$vegetation)==0, 1,0) 
F2$allmort<-(F2$insect+F2$fire+F2$weather+F2$disease+F2$animal+F2$unknowndamage+F2$vegetation)

Feqmerge<-F2

#Comparison model
F2<-subset(F2, select=c(PLT_CN, PREV_STOCKINGmid, PREV_DIAmean, PREV_Csum, DIAmean, Csum, STOCKINGmid,  REMPER, DIAchg, Cchg, STOCKchg, TEMP5growchg, PREC5growchg, SWdrought5chg, SWpluvial5chg, cutting, STDORGCD, cut, cutDIAmean, TPAsum, NOdisturb, alldisturb, allmort, insect, fire, weather, disease, animal, unknowndamage, vegetation, EHstocksum, PREVNMHstocksum, PREVSMHstocksum, PREVLHstocksum, PREVLCstocksum)) #, PREV_SPCDdiversity, SPCDdiversity, PREV_SPGRPCD, SPGRPCD))# 
c2<-40 # max diameter bin to collapse outliers
#previous whole dataset with change
#F2$PREV_DIAmean[F2$PREV_DIAmean>40]<-40
#F2$PREVDIAbin<-ceiling(F2$PREV_DIAmean*2)/2
F2$PREVDIAbin<-ifelse(F2$PREV_DIAmean<=32, ceiling(F2$PREV_DIAmean*2)/2, ceiling(F2$PREV_DIAmean/2)*2)
F2$PREVDIAbin[F2$PREVDIAbin> c2]<-c2  #collapse outliers

F2$PREVSTOCKbin<-ifelse(F2$PREV_DIAmean<=32, ceiling(F2$PREV_STOCKINGmid/10)*10,ceiling(F2$PREV_STOCKINGmid/20)*20) 
F2$PREVSTOCKbin[F2$PREVSTOCKbin>100]<-110  #collapse outliers

F2$PREVDIASTOCKbin<-F2$PREVDIAbin*10000+F2$PREVSTOCKbin
F2$mergeDIASTOCKbin<-F2$PREVDIAbin*10000+F2$PREVSTOCKbin
#F2$PREVDIASTOCKbinN<-ave(F2$PREVDIASTOCKbin>0, F2$PREVDIASTOCKbin, FUN=function(x) sum(x, na.rm=TRUE))
#F2$rando<-ave(F2$PREVDIASTOCKbin>0, F2$PREVDIASTOCKbin, FUN=function(x) sample(1:sum(x, na.rm=TRUE),sum(x, na.rm=TRUE),replace=FALSE))

#large bins
F2$PREVDIAbin2<-ceiling(F2$PREVDIAbin/2)*2
F2$PREVDIAbin2[F2$PREVDIAbin2> c2]<-c2 
F2$PREVSTOCKbin2<-ceiling(F2$PREVSTOCKbin/25)*25
F2$PREVSTOCKbin2[F2$PREVSTOCKbin2>100]<-100  #collapse outliers
F2$PREVDIASTOCKbin2<-F2$PREVDIAbin2*10000+F2$PREVSTOCKbin2
F2$mergeDIASTOCKbin2<-F2$PREVDIAbin2*10000+F2$PREVSTOCKbin2

#F4<-F2
F5<-subset(F2, cut==0 & STDORGCD==0 )
F5harvest<-subset(F2, STDORGCD==0) #, cutDIAmean>8.7 | is.na(cutDIAmean)==T )



time<-300 #years forward modeling
#Climate Colors:
pt<-"steelblue1" #pluvial trend
pe<-"royalblue"	#pluvial event
dt<-"goldenrod1" #drought trend
de<-"orange2"  	#drought event
nt<-"green" #no trend
av<-"grey75" #average
wt<-"peachpuff" #warming
ct<-"slategray1" #cooling
errorcol<-"grey75"
ud<-"greenyellow" #undisturbed
ds<-"forestgreen" #"#B78700" #50% less disturbance
lg<-"red" #logging
lg2<-"#B78700" #logging2
d175<-"plum"
d200<- "palevioletred4"

col6<-c("ct",  "de",  "wt",  "dt", "pe", "pt",  "av", "ud", "ds", "lg", "lg2", "n12", "n13") #, "n14", "n15") 
labels<-c("5-yr cool period",  "dry summer",  "5-yr warm period",  "5-yr dry period", "wet summer", "5-yr wet period",  "average", "undisturbed", "50% less disturbances", "average w/logging", "50% more disturbances", "75% more disturbances", "100% more disturbances" ) 
col5<-c(ct,de,wt,dt, pe, pt, av, ud, ds,lg, lg2,"plum", "palevioletred4" )

qq2<-c(seq(0.3,0.6,0.2),rep(0.6, 500))
qq3<-c(seq(0.9,0.4,-0.2), rep(0.4, 500))
Fbind2<-NULL
Fbindall2<-NULL

###
##for loop start
for(c in 7:7){
for(i in 1:(time/5)){

if(1==2)  #randomly sample for unstable solution testing/ model uncertainty 
{ 
F5$samp<-sample(1:sum(F5$PLT_CN>0, na.rm=T),sum(F5$PLT_CN>0, na.rm=T),replace=FALSE)
F2<-subset(F5, samp< (sum(F5harvest$PLT_CN>0)/2))
}

if(1==1)	
{
#use random 99% of underlying data to introduce slight randomness
F5$samp<-sample(1:sum(F5$PLT_CN>0, na.rm=T),sum(F5$PLT_CN>0, na.rm=T),replace=FALSE)
F5b<-subset(F5, samp<  (sum(F5$PLT_CN>0)*0.995 )) #7200
#F5c<-subset(F5, samp< (sum(F5$PLT_CN>0)*0.99))
F5harvest$samp<-sample(1:sum(F5harvest$PLT_CN>0, na.rm=T),sum(F5harvest$PLT_CN>0, na.rm=T),replace=FALSE)
F5harvestb<-subset(F5harvest, samp< (sum(F5harvest$PLT_CN>0)*0.99))

#climate conditions max. intensity with stable eq points
#more extreme conditions will not reach equilibrium, less extreme will produce less diferentiation from average. 
if(c==2){F2<-subset(F5b, SWdrought5chg<=quantile(F5$SWdrought5chg, qq3[i]+0.1, na.rm=T) & SWpluvial5chg<=quantile(F5$SWpluvial5chg, 0.85, na.rm=T))}
if(c==1){F2<-subset(F5b, TEMP5growchg<=quantile(F5$TEMP5growchg, qq3[i]-0.15, na.rm=T))}
if(c==3){F2<-subset(F5b, TEMP5growchg>=quantile(F5$TEMP5growchg, qq2[i]+0.1,  na.rm=T))} 
if(c==4){F2<-subset(F5b, PREC5growchg<=quantile(F5$PREC5growchg, qq3[i]+0.05, na.rm=T))}
if(c==5){F2<-subset(F5b, SWpluvial5chg>=quantile(F5$SWpluvial5chg, qq2[i]+0.15, na.rm=T) & SWdrought5chg>=quantile(F5$SWdrought5chg, 0.15, na.rm=T))}
if(c==6){F2<-subset(F5b, PREC5growchg>=quantile(F5$PREC5growchg, qq2[i]+0.15, na.rm=T)) }
if(c==7){F2<-subset(F5b)}
if(c==8){F2<-subset(F5b, NOdisturb==1)}
if(c==9){
	F5b$samp<-sample(1:sum(F5b$PLT_CN>0, na.rm=T),sum(F5b$PLT_CN>0, na.rm=T),replace=FALSE)
	F5b$samp<- ifelse(((F5b$insect+F5b$fire+F5b$animal+F5b$disease+F5b$weather) > 0), F5b$samp, 0)
	F2<-subset(F5b, (samp< (sum((F5b$insect+F5b$fire+F5b$animal+F5b$disease+F5b$weather)>0)*0.5))) #50% less disturbance
	}	
if(c==10){F2<-subset(F5harvestb)}
if(c==11){ #52% more disturbance
	F5b$samp<- ifelse(((F5b$insect+F5b$fire+F5b$animal+F5b$disease+F5b$weather) > 0), 0, F5b$samp)
	F2<-subset(F5b, (samp< (sum(F5$PLT_CN>0)*0.5)))
	} 

if(c==12){ #207% more disturbance
	F5b$samp<- ifelse(((F5b$insect+F5b$fire+F5b$animal+F5b$disease+F5b$weather) > 0), 0, F5b$samp)
	F2<-subset(F5b, (samp< (sum(F5$PLT_CN>0)*0.25)))
	}  
if(c==13){F2<-subset(F5b, (insect+fire+animal+disease+weather)>0) } #322% more disturbance (all forests disturbed)


#difference is warming increases disease and weather damage of early successional trees (lower density, smaller diameter), whereas weather and disease generally cause mortality of larger trees and at higher density (late succession) as occurs in the 5-year wet period!

}
print(c(c,i))
if(c==1 | c==3){meancond<-mean(F2$TEMP5growchg, na.rm=TRUE)}
if(c==4 | c==6){meancond <-mean(F2$PREC5growchg, na.rm=TRUE)}
if(c==5){meancond<-mean(F2$SWpluvial5chg, na.rm=TRUE)}
if(c==2){meancond<-mean(F2$SWdrought5chg, na.rm=TRUE)}
if(c>6){meancond<-NA}

F2$PREVDIASTOCKbinN<-ave(F2$PREVDIASTOCKbin>0, F2$PREVDIASTOCKbin, FUN=function(x) sum(x, na.rm=TRUE))
F2$rando<-ave(F2$PREVDIASTOCKbinN>0, F2$PREVDIASTOCKbin, FUN=function(x) sample(1:sum(x, na.rm=TRUE),sum(x, na.rm=TRUE),replace=FALSE ))
#f1<-subset(F2, !duplicated(mergeDIASTOCKbin), select=c(mergeDIASTOCKbin, PREVDIASTOCKbinN))
F2$PREVDIASTOCKbinN2<-ave(F2$PREVDIASTOCKbin2>0, F2$PREVDIASTOCKbin2, FUN=function(x) sum(x, na.rm=TRUE))
F2$rando2<-ave(F2$PREVDIASTOCKbinN2>0, F2$PREVDIASTOCKbin2, FUN=function(x) sample(1:sum(x, na.rm=TRUE),sum(x, na.rm=TRUE),replace=FALSE))

f1<-subset(F2, !duplicated(mergeDIASTOCKbin), select=c(mergeDIASTOCKbin, PREVDIASTOCKbinN))
f1c<-subset(F2, !duplicated(mergeDIASTOCKbin2), select=c(mergeDIASTOCKbin2, PREVDIASTOCKbinN2))

if(i==1){F4<-F5}
if(i>1){F4<-F3}
#print(qq2[i])
F3<-NULL
F3<-subset(F4) #, select= -c(PREVDIASTOCKbin, PREVDIASTOCKbinN))
##new vectors to bin and merge with old

F3$DIAbin<-ifelse(F3$DIAmean<=32, ceiling(F3$DIAmean*2)/2, ceiling(F3$DIAmean/2)*2)
F3$DIAbin[F3$DIAbin> c2]<-c2  #collapse outliers
F3$STOCKbin<-ifelse(F3$DIAmean<=32, ceiling(F3$STOCKINGmid/10)*10, ceiling(F3$STOCKINGmid/20)*20)  
F3$STOCKbin[F3$STOCKbin>100]<-110  #collapse outliers
F3$mergeDIASTOCKbin<-F3$DIAbin*10000+F3$STOCKbin
F3$mergeDIASTOCKbin[is.na(F3$mergeDIASTOCKbin)==T]<-130010 #restart plots that end with no tree
F3$DIASTOCKbin<-F3$DIAbin*10000+F3$STOCKbin
F3$DIASTOCKbin[is.na(F3$DIASTOCKbin)==T]<-130010
F3$DIASTOCKbinN<-ave(F3$DIASTOCKbin>0, F3$DIASTOCKbin, FUN=function(x) sum(x, na.rm=TRUE))

F3$DIAbin2<-ceiling(F3$DIAbin/2)*2
F3$DIAbin2[F3$DIAbin2> c2]<-c2
F3$STOCKbin2<-ceiling(F3$STOCKbin/25)*25
F3$STOCKbin2[F3$STOCKbin2>100]<-100  #collapse outliers
F3$DIASTOCKbin2<-F3$DIAbin2*10000+F3$STOCKbin2
F3$mergeDIASTOCKbin2<-F3$DIAbin2*10000+F3$STOCKbin2

F3$CsumSTART<-F3$Csum
F3$DIAmeanSTART<-F3$DIAmean
F3$STOCKINGmidSTART<-F3$STOCKINGmid

F3b<-merge(F3, f1, by=c("mergeDIASTOCKbin"))

#assign random number to each bin
#if rando number based on number of observations in matching dataset
F3b$rando<-ave(F3b$PREVDIASTOCKbinN, F3b$DIASTOCKbin, FUN=function(x) sample(1:mean(x, na.rm=TRUE),(sum(x, na.rm=TRUE)/mean(x, na.rm=TRUE)),replace=TRUE))

F3b<-subset(F3b, select=c(CsumSTART, DIAmeanSTART, STOCKINGmidSTART, mergeDIASTOCKbin, rando))
F3b<-merge(F3b, F2, by=c("mergeDIASTOCKbin", "rando"))
F3d<-subset(F3b, select=c(CsumSTART, DIAmeanSTART, STOCKINGmidSTART, STOCKINGmid, DIAmean, Csum, PLT_CN , PREV_DIAmean, PREV_Csum, Cchg, DIAchg, REMPER, TPAsum)) #, EHstocksum, NMHstocksum, SMHstocksum, LHstocksum))


##For bins without any plots matching initial and final surveys - use data from larger bins

F3c<-merge(F3, f1, all.x=T, by=c("mergeDIASTOCKbin"))

if(sum(is.na(F3c$PREVDIASTOCKbinN)==T)>0)
{
f3c<-subset(F3c, is.na(F3c$PREVDIASTOCKbinN)==T) #, select= -c( PREVDIASTOCKbinN2))
F3c<-merge(f3c, f1c, by=c("mergeDIASTOCKbin2"))
F3c$rando2<-ave(F3c$PREVDIASTOCKbinN2, F3c$DIASTOCKbin, FUN=function(x) sample(1:mean(x, na.rm=TRUE),(sum(x, na.rm=TRUE)/mean(x, na.rm=TRUE)),replace=TRUE))

F3c<-subset(F3c, select=c(CsumSTART, DIAmeanSTART, STOCKINGmidSTART, mergeDIASTOCKbin2, rando2))
F3c<-merge(F3c, F2, by=c("mergeDIASTOCKbin2", "rando2"))
setnames(F3c, c("mergeDIASTOCKbin2", "rando2"), c("mergeDIASTOCKbin", "rando") )

F3c<-subset(F3c, select=c(CsumSTART, DIAmeanSTART, STOCKINGmidSTART, STOCKINGmid, DIAmean, Csum, PLT_CN, PREV_DIAmean, PREV_Csum, Cchg, DIAchg, REMPER, TPAsum))# , EHstocksum, NMHstocksum, SMHstocksum, LHstocksum))
F3b<-subset(F3b, select=c(CsumSTART, DIAmeanSTART, STOCKINGmidSTART, STOCKINGmid, DIAmean, Csum, PLT_CN, PREV_DIAmean, PREV_Csum, Cchg, DIAchg, REMPER, TPAsum))# , EHstocksum, NMHstocksum, SMHstocksum, LHstocksum))
F3d<-rbind(F3b, F3c)
}
F3<-F3d

F<-subset(F3, !duplicated(PLT_CN))
plots<-sum(F$DIAmean>0, na.rm=T)

x1<-mean(F3$DIAmeanSTART, na.rm=TRUE)
x2<-mean(F3$DIAmean, na.rm=TRUE)
y1<-mean(F3$CsumSTART, na.rm=TRUE)
y2<-mean(F3$Csum, na.rm=TRUE)
TPA2<-mean(F3$TPAsum, na.rm=TRUE)
Caccum<-(mean(F3$Cchg, na.rm=TRUE))/(mean(F3$REMPER,na.rm=T))
#EH<- mean(F3$EHstocksum, na.rm=T)
#NMH<- mean(F3$NMHstocksum, na.rm=T)
#SMH<- mean(F3$SMHstocksum, na.rm=T)
#LH<- mean(F3$LHstocksum, na.rm=T)

#points(x2,y2, col="grey15", pch=(c-1), cex=0.3)
#points(x2,y2, col=col5[c], pch=c, cex=0.35)
#lines(c(x1, x2), c(y1, y2) ,col=("grey15"), lwd=0.25) #2 

#arrows(x1,y1, x2, y2 ,col=("grey15"), length=0.02, angle=25,lwd=0.85) #2 
#arrows(x1,y1, x2, y2 ,col=col5[c], length=0.018, angle=25,lwd=0.7) 
arrows(x1,y1, x2, y2 ,col=("grey15"), length=0.06, angle=25,lwd=1.5) #2 
arrows(x1,y1, x2, y2 ,col=col5[c], length=0.06, angle=25,lwd=1.35) 
arrows(x1,y1, x2, y2 ,col=col5[c], length=0.06, angle=25,lwd=1.35) 


if( i==10 | i==20 ) 
{
	arrows(x1,y1, x2, y2 ,col=("salmon"), length=0.06, angle=25,lwd=1.2)	 #denote 100 years of growth with red vectors
		}

# if(i==time/5 & c<12)
# {
	# points(x2,y2, col="white", pch=(c-1), lwd=1.45, cex=1.35)
 	# points(x2,y2, col="black", pch=(c-1), lwd=1.75, cex=1.25)
	# points(x2,y2, col=col5[c], pch=(c-1), lwd=1.25, cex=1.25)
	 # }
	 
# if(i==time/5 & c>=12)
# {
	# points(x2,y2, col="white", pch=(c), lwd=1.45, cex=1.35)
 	# points(x2,y2, col="black", pch=(c), lwd=1.75, cex=1.25)
	# points(x2,y2, col=col5[c], pch=(c), lwd=1.25, cex=1.25)
	 # } 	 
	 
	 
#gather arrow data
Fbindall<-c(x1,y1, x2, y2, TPA2, c, Caccum, i) #, EH, NMH, SMH, LH)
Fbindall2<-rbind(Fbindall2, Fbindall)

#turn error bars on last increment # dont plot
if(i==time/5){ 
Cchg<-F$Cchg
DIAchg<-F$DIAchg
#Plot bootstraped error associated with mean diameter and stem density change 
#not actually used anywhere...
endCerror<-boot(Cchg, meanFunc, R=1000)	
endCerror<-sqrt(var(endCerror$t))*1.96
endQMDerror<-boot(DIAchg, meanFunc, R=1000)				
endQMDerror<-sqrt(var(endQMDerror$t))*1.96

Fbind<-c(x2, y2, TPA2, endCerror, endQMDerror, col6[c], Caccum, meancond, col5[c])
Fbind2<-rbind(Fbind2, Fbind)
print(Fbind2)


}
}
#i
}
#box()





#plot vectors on carbon mapping
Fbindall2<-as.data.table(Fbindall2)
setnames(Fbindall2, c("x1", "y1","x2", "y2", "TPA2", "c", "Caccum", "i"))#, "EH", "NMH", "SMH", "LH"))


Fbind3<-as.data.table(Fbind2)
setnames(Fbind3,c("x2", "y2", "TPA2", "endCerror", "endQMDerror", "cond", "Caccum", "meancond", "col5"))

#write.csv(Fbindall2, file= "/Users/travis/Desktop/futureCFbindall2.csv")
#write.csv(Fbind3, file= "/Users/travis/Desktop/futureCFbind3.csv")
#Fbindall2 <- as.data.table(read.csv("/Users/travis/Desktop/futureCFbindall2.csv"))


##
######
# replot #########
labs3<- c("ct",  "de",  "wt",  "dt",  "pe",  "pt",  "av", "ud",  "d0.5x",  "lg",  "d1.5x", "d2x", "d3x")
labels3<-c("(ud) Undisturbed forests", "(d0.5x) 50% less disturbances", expression(paste("(dt) -0.24 mm day"^"-1","precipitation" )), expression(paste("(ct) -0.35"^"o","C cooling" )), "(de) Recurring dry summers", "(av) Average of all forests", expression(paste("(wt) +0.33"^"o","C warming" )), expression(paste("(pt) +0.37 mm day"^"-1","precipitation" )), "(pe) Recurring wet summers", "(d1.5x) 50% more disturbances", "(d2x) 100% more disturbances", "(d3x) 200% more disturbances", "(lg) All forests w/ logging", "  5-year future projection"   )

legend(22.5, 170, labels3, pch= c(rep(21,13), NA), col=c(ud, ds, dt, ct, de, av, wt, pt, pe, lg2, d175, d200, lg, NA  ), cex=0.85, bty="o", pt.cex=0.95, pt.lwd=2, pt.bg="black", ncol=2, bg="#FFFFFF90", box.col="#FFFFFF90")
arrows(28.85, 137.5, 29.35, 137.5,col="black", length=0.06, angle=25,lwd=1.5 )
arrows(28.85, 137.5, 29.35, 137.5,col="grey25", length=0.06, angle=25,lwd=1.25 )
text(Fbindall2$x2[1]-1.75, Fbindall2$y2[1]-10, "Eastern US\naverage", cex=0.85, col="black", pos=4)

y3<- Fbindall2$y2[Fbindall2$i==60]
x3<- Fbindall2$x2[Fbindall2$i==60]
lm3<-lm(y3 ~ x3)
summary(lm3)
abline(lm3, col="#00000070", lwd=3)


c1<-(1:13)
foreach(r=1:780) %do%{
#fb2<-subset(Fbindall2, c==c1[r])
fb2<- Fbindall2[r]
arrows(fb2$x1, fb2$y1, fb2$x2, fb2$y2,col="black", length=0.06, angle=25,lwd=1.5 )
arrows(fb2$x1, fb2$y1, fb2$x2, fb2$y2,col=col5[fb2$c], length=0.06, angle=25,lwd=1.33 )
arrows(fb2$x1, fb2$y1, fb2$x2, fb2$y2,col=col5[fb2$c], length=0.06, angle=25,lwd=1.33 )
#Fbindall3<-subset(fb2, i==10 | i==20 )
if(fb2$i==10 | fb2$i==20){arrows(Fbindall3$x1, Fbindall3$y1, Fbindall3$x2, Fbindall3$y2,col="salmon", length=0.06, angle=25,lwd=1.35 )}
}

points(x3, y3, pch=19, col="black", cex=0.75)
xoff<-c(rep(0.6,10), -0.62, 0.6,0.6)
xoff2<-c(rep(0.6,10), 0, 0.6,0.6)
yoff<- c(0, 0, -3, 5, -2, -2, 0, 0,0,0,-6,0,0)
yoff2<- c(0, 0, -3, 5, -2, -2, 0, 0,0,0,-4.5,0,0)
segments(x3, y3, x3+xoff2, y3+yoff2)
text(x3+(xoff-0.1), y3+yoff, labs3, pos=4, cex=0.85)
box()

##############END CARBON FORWARD PROJECTION
##

##############END CARBON FORWARD PROJECTION
##
#
#






#
############### Carbon Empirical GROWTH Model - UNCERTAINTY ANALYSIS ##################################################################################
meanFunc<-function(x,i){mean(x[i], na.rm=TRUE)}

#subset resurveyed plots to those with at least one 5" tree in order to eliminate plots only actually surveyed once.
F2<-subset(RM5mergedfullcut1991.11.17.2014, PREVcarbon1sum>0) # & STDORGCD==0 ) 

#carbonC
F2$DIAmean<-(F2$DIAendmean* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean<-(F2$DIAbeginmean *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum<-(F2$endsumTPA /0.404686)			#number of trees per hectare by plot at resurvey
F2$PREV_TPAsum<-(F2$startsumTPA /0.404686)	#number of trees per hectare by plot at initial survey (cm)
F2$Csum<-F2$carbon1sum* 0.00112085116  #	#convert living carbon to Mg/ha
F2$PREV_Csum<-F2$PREVcarbon1sum* 0.00112085116 ##convert living carbon to Mg/ha
F2$PREV_STOCKINGmid<-F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey
F2$STOCKINGmid<-F2$STOCKING5mid


#calculate change in Stem density and mean tree diameter for each plot
F2$Cchg<-F2$Csum-F2$PREV_Csum
F2$DIAchg<-F2$DIAmean-F2$PREV_DIAmean
F2$STOCKchg<-F2$STOCKINGmid-F2$PREV_STOCKINGmid
#set grid cell size
a<-b<-2  #2=0.5degrees, 0.5=2degrees
	F2$iLAT<-((ceiling(F2$LAT*a))/a)-0.5 #+0.25  #round down to integer then add .25 ; switched /* removed ).25 addition
	F2$iLON<-((floor(F2$LON*b))/b)+0.5 #-0.25  #round down to integer then subtract .25; switched /*  ).25 addition
	F2$LATLONYR<-(F2$iLAT*10000000000000) +(F2$iLON*(-100000000)) +F2$INVYR
	F2$LATLON<-as.factor(F2$iLAT*100000000+ F2$iLON*-1000) #LATLON identifyier

#Find climate variable mean and range for each grid cell
F2$TEMPminLL<-ave(F2$temp5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$TEMPgrowminLL<-ave(F2$temp5grow, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$TEMP5chg<-F2$temp5-F2$TEMPminLL
F2$TEMP5growchg<-F2$temp5grow-F2$TEMPgrowminLL

F2$PRECminLL<-ave(F2$PREC5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$PRECgrowminLL<-ave(F2$PREC5grow, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$PREC5chg<-F2$PREC5-F2$PRECminLL
F2$PREC5growchg<-F2$PREC5grow-F2$PRECgrowminLL

F2$SWdroughtLL<-ave(F2$SWdrought5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$SWpluvialLL<-ave(F2$SWpluvial5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$SWdrought5chg<-F2$SWdrought5-F2$SWdroughtLL
F2$SWpluvial5chg<-F2$SWpluvial5-F2$SWpluvialLL

F2$NOdisturb<-ifelse((F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)==0,1,0) 
F2$alldisturb<-ifelse((F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)>0 & (F2$unknowndamage+F2$vegetation)==0, 1,0) 
F2$allmort<-(F2$insect+F2$fire+F2$weather+F2$disease+F2$animal+F2$unknowndamage+F2$vegetation)

Feqmerge<-F2
 
#Comparison model
F2<-subset(F2, select=c(PLT_CN, PREV_STOCKINGmid, PREV_DIAmean, PREV_Csum, DIAmean, Csum, STOCKINGmid,  REMPER, DIAchg, Cchg, STOCKchg, TEMP5growchg, PREC5growchg, SWdrought5chg, SWpluvial5chg, cutting, STDORGCD, cut, cutDIAmean, TPAsum, NOdisturb, alldisturb, allmort, insect, fire, animal, disease, weather)) #, PREV_SPCDdiversity, SPCDdiversity, PREV_SPGRPCD, SPGRPCD))# 
c2<-40 # max diameter bin to collapse outliers
#previous whole dataset with change

F2$PREVDIAbin<-ifelse(F2$PREV_DIAmean<=32, ceiling(F2$PREV_DIAmean*2)/2, ceiling(F2$PREV_DIAmean/2)*2)
F2$PREVDIAbin[F2$PREVDIAbin> c2]<-c2  #collapse outliers

F2$PREVSTOCKbin<-ifelse(F2$PREV_DIAmean<=32, ceiling(F2$PREV_STOCKINGmid/10)*10,ceiling(F2$PREV_STOCKINGmid/20)*20) 
F2$PREVSTOCKbin[F2$PREVSTOCKbin>100]<-110  #collapse outliers

F2$PREVDIASTOCKbin<-F2$PREVDIAbin*10000+F2$PREVSTOCKbin
F2$mergeDIASTOCKbin<-F2$PREVDIAbin*10000+F2$PREVSTOCKbin
#F2$PREVDIASTOCKbinN<-ave(F2$PREVDIASTOCKbin>0, F2$PREVDIASTOCKbin, FUN=function(x) sum(x, na.rm=TRUE))
#F2$rando<-ave(F2$PREVDIASTOCKbin>0, F2$PREVDIASTOCKbin, FUN=function(x) sample(1:sum(x, na.rm=TRUE),sum(x, na.rm=TRUE),replace=FALSE))

#large bins
F2$PREVDIAbin2<-ceiling(F2$PREVDIAbin/2)*2
F2$PREVDIAbin2[F2$PREVDIAbin2> c2]<-c2 
F2$PREVSTOCKbin2<-ceiling(F2$PREVSTOCKbin/25)*25
F2$PREVSTOCKbin2[F2$PREVSTOCKbin2>100]<-100  #collapse outliers
F2$PREVDIASTOCKbin2<-F2$PREVDIAbin2*10000+F2$PREVSTOCKbin2
F2$mergeDIASTOCKbin2<-F2$PREVDIAbin2*10000+F2$PREVSTOCKbin2

#F4<-F2
F5<-subset(F2, cut==0 & STDORGCD==0 )
F5harvest<-subset(F2, STDORGCD==0) #


time<-300 #years forward modeling

#Climate Colors:
pt<-"steelblue1" #pluvial trend
pe<-"royalblue"	#pluvial event
dt<-"goldenrod1" #drought trend
de<-"orange2"  	#drought event
nt<-"green" #no trend
av<-"grey75" #average
wt<-"peachpuff" #warming
ct<-"slategray1" #cooling
errorcol<-"grey75"
ud<-"greenyellow" #undisturbed
ds<-"brown" #"#B78700"
lg<-"red" #logging
lg2<-"#B78700" #logging2

#col6<-c("ct",  "de",  "wt",  "dt", "pe", "pt",  "av", "ud", "ds", "lg", "lg2")  #use col6 from above and labels and col5
#labels<-c("5-yr cool period",  "dry summer",  "5-yr warm period",  "5-yr dry period", "wet summer", "5-yr wet period",  "average", "undisturbed", "disturbed", "average w/logging", "all mortality" ) 
#col5<-c(ct,de,wt,dt, pe, pt, av, ud, ds,lg, lg2)


qq2<-c(seq(0.3,0.6,0.2),rep(0.6, 500))
qq3<-c(seq(0.9,0.4,-0.2), rep(0.4, 500))
Fbindeq<-NULL
Fbindeq2<-NULL
Fbindeq3<-NULL

###
##for loop start

#ERROReqbind<-foreach(r=1:10, .combine=rbind) %dopar%{  #iterate model scenario 10 times and rbind results
for(r in 1:10) {  #iterate model scenario 10 times and rbind results #CAUTION THIS TAKES A LOONNG TIME  ~20 minutes - consider foreach parrallel proccessing

for(c in 1:13){   #run each scenario
if(r==1 & c==1){print("Run countdown:")}
	print(c(10-r, 11-c))

for(i in 1:(time/5))  #run each model scenario to steady-state and collect last 10 iterations
{

#same climate extremes as above projections, now just repeated 10 times each to assess model uncertainty 
if(c==2){F2<-subset(F5, SWdrought5chg<=quantile(F5$SWdrought5chg, qq3[i]+0.1, na.rm=T) & SWpluvial5chg<=quantile(F5$SWpluvial5chg, 0.85, na.rm=T))}
if(c==1){F2<-subset(F5, TEMP5growchg<=quantile(F5$TEMP5growchg, qq3[i]-0.15, na.rm=T))}
if(c==3){F2<-subset(F5, TEMP5growchg>=quantile(F5$TEMP5growchg, qq2[i]+0.1,  na.rm=T))} 
if(c==4){F2<-subset(F5, PREC5growchg<=quantile(F5$PREC5growchg, qq3[i]+0.05, na.rm=T))}
if(c==5){F2<-subset(F5, SWpluvial5chg>=quantile(F5$SWpluvial5chg, qq2[i]+0.15, na.rm=T) & SWdrought5chg>=quantile(F5$SWdrought5chg, 0.15, na.rm=T))}
if(c==6){F2<-subset(F5, PREC5growchg>=quantile(F5$PREC5growchg, qq2[i]+0.15, na.rm=T)) }
if(c==7){F2<-subset(F5)}
if(c==8){F2<-subset(F5, NOdisturb==1)}
#if(c==9){F2<-subset(F5, alldisturb>0)}
#if(c==10){F2<-subset(F5harvest)}
#if(c==11){F2<-subset(F5, allmort>0)}	

if(c==9){
	F5$samp<-sample(1:sum(F5$PLT_CN>0, na.rm=T),sum(F5$PLT_CN>0, na.rm=T),replace=FALSE)
	F5$samp<- ifelse(((F5$insect+F5$fire+F5$animal+F5$disease+F5$weather) > 0), F5$samp, 0)
	F2<-subset(F5, (samp< (sum((F5$insect+F5$fire+F5$animal+F5$disease+F5$weather)>0)*0.5))) #50% less disturbance
	}	
if(c==10){F2<-subset(F5harvest)}
if(c==11){
	F5$samp<-sample(1:sum(F5$PLT_CN>0, na.rm=T),sum(F5$PLT_CN>0, na.rm=T),replace=FALSE)
	F5$samp<- ifelse(((F5$insect+F5$fire+F5$animal+F5$disease+F5$weather) > 0), 0, F5$samp)
	F2<-subset(F5, (samp< (sum(F5$PLT_CN>0)*0.5)))
	} 

if(c==12){
	F5$samp<-sample(1:sum(F5$PLT_CN>0, na.rm=T),sum(F5$PLT_CN>0, na.rm=T),replace=FALSE)
	F5$samp<- ifelse(((F5$insect+F5$fire+F5$animal+F5$disease+F5$weather) > 0), 0, F5$samp)
	F2<-subset(F5, (samp< (sum(F5$PLT_CN>0)*0.25)))
	}  
if(c==13){F2<-subset(F5, (insect+fire+animal+disease+weather)>0) }




F2$samp<-sample(1:sum(F2$PLT_CN>0, na.rm=T),sum(F2$PLT_CN>0, na.rm=T),replace=FALSE)
F2<-subset(F2, samp<=(7582*0.95))  #sample 95% of smallest condion

F2$PREVDIASTOCKbinN<-ave(F2$PREVDIASTOCKbin>0, F2$PREVDIASTOCKbin, FUN=function(x) sum(x, na.rm=TRUE))
F2$rando<-ave(F2$PREVDIASTOCKbinN>0, F2$PREVDIASTOCKbin, FUN=function(x) sample(1:sum(x, na.rm=TRUE),sum(x, na.rm=TRUE),replace=FALSE ))
#f1<-subset(F2, !duplicated(mergeDIASTOCKbin), select=c(mergeDIASTOCKbin, PREVDIASTOCKbinN))
F2$PREVDIASTOCKbinN2<-ave(F2$PREVDIASTOCKbin2>0, F2$PREVDIASTOCKbin2, FUN=function(x) sum(x, na.rm=TRUE))
F2$rando2<-ave(F2$PREVDIASTOCKbinN2>0, F2$PREVDIASTOCKbin2, FUN=function(x) sample(1:sum(x, na.rm=TRUE),sum(x, na.rm=TRUE),replace=FALSE))

f1<-subset(F2, !duplicated(mergeDIASTOCKbin), select=c(mergeDIASTOCKbin, PREVDIASTOCKbinN))
f1c<-subset(F2, !duplicated(mergeDIASTOCKbin2), select=c(mergeDIASTOCKbin2, PREVDIASTOCKbinN2))

if(i==1){F4<-F5}
if(i>1){F4<-F3}
#print(qq2[i])
F3<-NULL
F3<-F4 #, select= -c(PREVDIASTOCKbin, PREVDIASTOCKbinN))

##new vectors to bin and merge with old
F3$DIAbin<-ifelse(F3$DIAmean<=32, ceiling(F3$DIAmean*2)/2, ceiling(F3$DIAmean/2)*2)
F3$DIAbin[F3$DIAbin> c2]<-c2  #collapse outliers
F3$STOCKbin<-ifelse(F3$DIAmean<=32, ceiling(F3$STOCKINGmid/10)*10, ceiling(F3$STOCKINGmid/20)*20)  
F3$STOCKbin[F3$STOCKbin>100]<-110  #collapse outliers
F3$mergeDIASTOCKbin<-F3$DIAbin*10000+F3$STOCKbin
F3$mergeDIASTOCKbin[is.na(F3$mergeDIASTOCKbin)==T]<-130010 #restart plots that end with no tree
F3$DIASTOCKbin<-F3$DIAbin*10000+F3$STOCKbin
F3$DIASTOCKbin[is.na(F3$DIASTOCKbin)==T]<-130010
F3$DIASTOCKbinN<-ave(F3$DIASTOCKbin>0, F3$DIASTOCKbin, FUN=function(x) sum(x, na.rm=TRUE))

F3$DIAbin2<-ceiling(F3$DIAbin/2)*2
F3$DIAbin2[F3$DIAbin2> c2]<-c2
F3$STOCKbin2<-ceiling(F3$STOCKbin/25)*25
F3$STOCKbin2[F3$STOCKbin2>100]<-100  #collapse outliers
F3$DIASTOCKbin2<-F3$DIAbin2*10000+F3$STOCKbin2
F3$mergeDIASTOCKbin2<-F3$DIAbin2*10000+F3$STOCKbin2

F3$CsumSTART<-F3$Csum
F3$DIAmeanSTART<-F3$DIAmean
F3$STOCKINGmidSTART<-F3$STOCKINGmid

F3b<-merge(F3, f1, by=c("mergeDIASTOCKbin"))

#assign random number to each bin
#if rando number based on number of observations in matching dataset
F3b$rando<-ave(F3b$PREVDIASTOCKbinN, F3b$DIASTOCKbin, FUN=function(x) sample(1:mean(x, na.rm=TRUE),(sum(x, na.rm=TRUE)/mean(x, na.rm=TRUE)),replace=TRUE))

F3b<-subset(F3b, select=c(CsumSTART, DIAmeanSTART, STOCKINGmidSTART, mergeDIASTOCKbin, rando))
F3b<-merge(F3b, F2, by=c("mergeDIASTOCKbin", "rando"))
F3d<-subset(F3b, select=c(CsumSTART, DIAmeanSTART, STOCKINGmidSTART, STOCKINGmid, DIAmean, Csum, PLT_CN , PREV_DIAmean, PREV_Csum, Cchg, DIAchg, REMPER, TPAsum))

##For bins without any plots matching initial and final surveys - use data from larger bins
F3<-F3d
x2<-mean(F3$DIAmean, na.rm=TRUE)
y2<-mean(F3$Csum, na.rm=TRUE)
TPA2<-mean(F3$TPAsum, na.rm=TRUE)
STOCKING2<-mean(F3$STOCKINGmid, na.rm=TRUE)
#gather eq data
if(i>=45)
	{
	Fbindeq<-c( x2, y2, TPA2, STOCKING2, c, i)
	Fbindeq2<-rbind(Fbindeq2, Fbindeq)
	}
	
} # i loop
} # c loop
} # r loop



Fbindeq<-as.data.table(Fbindeq2)
setnames(Fbindeq,c("x2", "y2", "TPA2", "STOCKING2", "cond", "i"))
aggregate(y2~cond,data= Fbindeq,  sd)
aggregate(x2~cond,data= Fbindeq,  sd)
aggregate(TPA2~cond,data= Fbindeq,  sd)
aggregate(STOCKING2~cond,data= Fbindeq,  sd)
Fbindeq$cSD<-(ave(Fbindeq$y2, Fbindeq$cond, FUN=function(x) sd(x, na.rm=TRUE)))
Fbindeq$xSD<-(ave(Fbindeq$x2, Fbindeq$cond, FUN=function(x) sd(x, na.rm=TRUE)))
Fbindeq$tpaSD<-(ave(Fbindeq$TPA2, Fbindeq$cond, FUN=function(x) sd(x, na.rm=TRUE)))

Feq2<-subset(Fbindeq, !duplicated(cond), select=c(cSD))
cond<-col6 #as.data.table(c("ct",  "de",  "wt",  "dt", "pe", "pt",  "av", "ud", "ds", "lg", "lg2" ) )
Feq3<-cbind(Feq2, cond)
setnames(Feq3,c("ySD", "cond"))
write.csv(Feq3, file= "/Users/travis/Desktop/futureCFeq3.csv")

#######END PREJECTION ERROR ANALYSIS######################################
#
#
#




#
############## CLimate stats -  using data from projections above! ################
#########
###

#ratio of mean conditions with Null weighted spatially by Condition
meanFunc2<-function(x,i){(mean(x$variCOND[i], na.rm=TRUE)/mean((x$variNULL[i]), na.rm=TRUE))} 
meanFunc3<-function(x,i){(mean(x$deadcarbon5sumTREE[i], na.rm=TRUE)/mean((x$mortality2[i]), na.rm=TRUE))} 
 


F2<-subset(RM5mergedfullcut1991.11.17.2014, PREVcarbon1sum>0 & PREVSTOCKING5mid>= 60  ) 
#Convert variables to standard names and metric units
F2$DIAmean<-(F2$DIA5meanalive* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean<-(F2$PREVDIA5meanalive *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum<-(F2$endsumTPA /0.404686)			#number of trees per hectare by plot at resurvey
F2$PREV_TPAsum<-(F2$startsumTPA /0.404686)	#number of trees per hectare by plot at initial survey (cm)
F2$PREV_STOCKINGmid<-F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey
F2$cutDIAdif<- (F2$DIA5meanalive-F2$cutDIAmean)
F2$mortDIAdif<- (F2$DIA5meanalive-F2$mortDIAmean)

#F2$mortDIAdif[is.na(F2$mortDIAdif)==T]<-0  #when no mortality, no change from average
 
#set grid cell size
a<-b<-4  #2=0.5degrees, 0.5=2degrees
	F2$iLAT<-((ceiling(F2$LAT*a))/a)-0.5 #+0.25  #round down to integer then add .25 ; switched /* removed ).25 addition
	F2$iLON<-((floor(F2$LON*b))/b)+0.5 #-0.25  #round down to integer then subtract .25; switched /*  ).25 addition
	F2$LATLONYR<-(F2$iLAT*10000000000000) +(F2$iLON*(-100000000)) +F2$INVYR
	F2$LATLON<-as.factor(F2$iLAT*100000000+ F2$iLON*-1000) #LATLON identifyier

#Find climate variable mean and range for each grid cell
F2$TEMPminLL<-ave(F2$temp5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$TEMPgrowminLL<-ave(F2$temp5grow, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$TEMP5chg<-F2$temp5-F2$TEMPminLL
F2$TEMP5growchg<-F2$temp5grow-F2$TEMPgrowminLL

F2$PRECminLL<-ave(F2$PREC5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$PRECgrowminLL<-ave(F2$PREC5grow, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$PREC5chg<-F2$PREC5-F2$PRECminLL
F2$PREC5growchg<-F2$PREC5grow-F2$PRECgrowminLL

F2$SWdroughtLL<-ave(F2$SWdrought5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$SWpluvialLL<-ave(F2$SWpluvial5, F2$LATLON, FUN=function(x) mean(x, na.rm=TRUE))
F2$SWdrought5chg<-F2$SWdrought5-F2$SWdroughtLL
F2$SWpluvial5chg<-F2$SWpluvial5-F2$SWpluvialLL

F2$average<-1
F2$nsampTOT<-ave(F2$DIAmean>0,F2$LATLON, FUN=function(x) sum(x, na.rm=TRUE))
F2$NOdisturb<-ifelse((F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)==0,1,NA) 
F2$alldisturb<-ifelse((F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)>0 & (F2$unknowndamage+F2$vegetation)==0, 1, NA)
F2$allmort<-ifelse((F2$insect+F2$fire+F2$weather+F2$disease+F2$animal+F2$unknowndamage+F2$vegetation)>0, 1, NA)
F2$disturbsum<- (F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)

vari<-c( "TEMP5growchg", "SWdrought5chg",  "TEMP5growchg",  "PREC5growchg", "SWpluvial5chg", "PREC5growchg", "average", "NOdisturb", "alldisturb", "allmort", "cutting") 

cond<- col6 # from above #c("ct",  "de",  "wt",  "dt", "pe", "pt",  "av", "ud", "ds", "lg2", "lg" ) 
vari1<-subset(F2, select=c(noquote(vari))) 
#vari2<-c("fire")

F5<-subset(F2, cut==0 & STDORGCD==0) #& mortDIAmean>0)
F5harvest<-subset(F2, STDORGCD==0) # & cutDIAmean>0) #, cutDIAmean>8.7 | is.na(cutDIAmean)==T )

qq3[i]<-0.4
qq2[i]<-0.6
F5b<-NULL
for(c in c(1:13))
{
	
if(c==1){
	f22<-subset(F5, TEMP5growchg<=quantile(F5$TEMP5growchg, qq3[i]-0.15, na.rm=T))
	F33<-subset(F5, LATLON %in% (f22$LATLON) & TEMP5growchg<=quantile(F5$TEMP5growchg, 0.7, na.rm=T) & TEMP5growchg>=quantile(F5$TEMP5growchg, 0.3, na.rm=T))
	} # Do all null conditions individually to simplify
		
#same adjustment as above for largest effect with stable eq points
if(c==2){
	f22<-subset(F5, SWdrought5chg<=quantile(F5$SWdrought5chg, qq3[i]+0.1, na.rm=T) & SWpluvial5chg<=quantile(F5$SWpluvial5chg, 0.85, na.rm=T))
	F33<-subset(F5, LATLON %in% (f22$LATLON) & SWdrought5chg<=quantile(F5$SWdrought5chg, 0.7, na.rm=T)  & SWdrought5chg>=quantile(F5$SWdrought5chg, 0.3, na.rm=T) & SWpluvial5chg<=quantile(F5$SWpluvial5chg, 0.85, na.rm=T))  # SWdrought5chg above zero is meaningless 
	}
	
if(c==3){
	f22<-subset(F5, TEMP5growchg>=quantile(F5$TEMP5growchg, qq2[i]+0.1,  na.rm=T))
	F33<-subset(F5, LATLON %in% (f22$LATLON) & TEMP5growchg<=quantile(F5$TEMP5growchg, 0.7,  na.rm=T) & TEMP5growchg>=quantile(F5$TEMP5growchg, 0.3,  na.rm=T))
	} 
if(c==4){
	f22<-subset(F5, PREC5growchg<=quantile(F5$PREC5growchg, qq3[i]+0.05, na.rm=T))
	F33<-subset(F5, LATLON %in% (f22$LATLON) & PREC5growchg<=quantile(F5$PREC5growchg, 0.7, na.rm=T) & PREC5growchg>=quantile(F5$PREC5growchg, 0.3, na.rm=T))
	}
if(c==5){
	f22<-subset(F5, SWpluvial5chg>=quantile(F5$SWpluvial5chg, qq2[i]+0.15, na.rm=T) & SWdrought5chg>=quantile(F5$SWdrought5chg, 0.15, na.rm=T))
	F33<-subset(F5, LATLON %in% (f22$LATLON) & SWpluvial5chg<=quantile(F5$SWpluvial5chg, 0.7, na.rm=T)  & SWpluvial5chg>quantile(F5$SWpluvial5chg, 0.3, na.rm=T) & SWdrought5chg>=quantile(F5$SWdrought5chg, 0.15, na.rm=T)) # SWpluvial5chg below zero is meaningless
	}

if(c==6){
	f22<- subset(F5, PREC5growchg>=quantile(F5$PREC5growchg, qq2[i]+0.15, na.rm=T)) 
	F33<- subset(F5, LATLON %in% (f22$LATLON) & PREC5growchg<=quantile(F5$PREC5growchg, 0.7, na.rm=T) & PREC5growchg>=quantile(F5$PREC5growchg, 0.3, na.rm=T)) 
	}
if(c==7){
	f22	<- subset(F5)
	F33<-subset(F5, LATLON %in% (f22$LATLON))
	}
if(c==8){
	f22<-subset(F5, NOdisturb==1)
	F33<-subset(F5, LATLON %in% (f22$LATLON))
	}
if(c==9){
	F5$samp<-sample(1:sum(F5$PLT_CN>0, na.rm=T),sum(F5$PLT_CN>0, na.rm=T),replace=FALSE)
	F5$samp<- ifelse(((F5$insect+F5$fire+F5$animal+F5$disease+F5$weather) > 0), F5$samp, 0)
	f22<-subset(F5, (samp< (sum((F5$insect+F5$fire+F5$animal+F5$disease+F5$weather)>0)*0.5))) #50% less disturbance
	F33<-subset(F5, LATLON %in% (f22$LATLON) )
	} 
 if(c==10){
 	f22<-subset(F5harvest)
 	F33<-subset(F5, LATLON %in% (f22$LATLON))
 	}
if(c==11){
	F5$samp<-sample(1:sum(F5$PLT_CN>0, na.rm=T),sum(F5$PLT_CN>0, na.rm=T),replace=FALSE)
	F5$samp<- ifelse(((F5$insect+F5$fire+F5$animal+F5$disease+F5$weather) > 0), 0, F5$samp)
	f22<-subset(F5, (samp< (sum(F5$PLT_CN>0)*0.5)))
	F33<-subset(F5, LATLON %in% (f22$LATLON))
	} 
if(c==12){
	F5$samp<-sample(1:sum(F5$PLT_CN>0, na.rm=T),sum(F5$PLT_CN>0, na.rm=T),replace=FALSE)
	F5$samp<- ifelse(((F5$insect+F5$fire+F5$animal+F5$disease+F5$weather) > 0), 0, F5$samp)
	f22<-subset(F5, (samp< (sum(F5$PLT_CN>0)*0.25)))
	F33<-subset(F5, LATLON %in% (f22$LATLON) )
	}
if(c==13){
	f22<-subset(F5, (insect+fire+animal+disease+weather)>0) 
	F33<-subset(F5, LATLON %in% (f22$LATLON))
	}


if(c!=10)
{
Vgrow<- mean(f22$mortDIAmean, na.rm=T)*2.54 #
Vgrow2<- mean(F33$mortDIAmean, na.rm=T)*2.54 #divide by CtoNratio incase it does not = 1
error2<-std.error(f22$mortDIAmean, na.rm=T)*2.54
#Vgrow3<- mean(f22$mortDIAdif, na.rm=T)*10 
}

if(c==10)
{
Vgrow<-((mean(f22$cutDIAmean*f22$cut, na.rm=T)+mean(f22$mortDIAmean*f22$mortality2, na.rm=T))/(mean(f22$mortality2[f22$mortality2>0])+mean(f22$cut[f22$cut>0])))*2.54 #
Vgrow2<-mean(F33$mortDIAmean, na.rm=T)*2.54 
error2<-std.error(f22$cutDIAmean*2.54, na.rm=T)
}

Vgrow3<-Vgrow-Vgrow2 #
n<-round(sum(f22$mortality2/6.01)/sum(f22$startsumTPA/6.01), 1)*100
pch1<- ifelse(c<12, c-1, c)

f34<-c(cond[c], n, round(Vgrow2,2), round(Vgrow,2), round(Vgrow3, 2), error2, pch1) #vari2[v2],
F5b<-rbind(F5b, f34)
print(c(c))
} 

F6b<-as.data.table(F5b)
setnames(F6b, c( "cond", "n","Vgrow2","Vdead", "difference", "error2", "pch1")) #"Disturbance type",
print(na.omit(F6b))
#write.csv(na.omit(F6b))

#Fbind3 <- as.data.table(read.csv("/Users/travis/Desktop/futureCFbind3.csv"))
Fbind4<-as.data.table(Fbind3)
setnames(Fbind4,c("x2", "y2mean", "TPA2", "endCerror", "endQMDerror", "cond", "Caccum", "meancond", "col5"))

FF<-as.data.table(merge(Fbind4, F6b, by="cond")) #merge forward carbon projections with mean dia mortality data

FF<-as.data.table(merge(FF, Feq3, by="cond")) #Using model error data (Feq3) derived above
#write.csv(FF, file= "/Users/travis/Desktop/futureCFF.csv")
#FF <- as.data.table(read.csv("/Users/travis/Desktop/futureCFF.csv"))


FF2<-subset(FF) #, cond!="ct") 
y<-as.numeric(as.character(FF2$y2mean))
x<-as.numeric(as.character(FF2$difference))+as.numeric(as.character(FF2$Vdead[1]))
x[is.na(x)]<-as.numeric(as.character(FF2$Vdead[1]))
xerror<-(as.numeric(as.character(FF2$error2)))*1.96
yerror<-FF2$ySD*1.96
cond<-as.character(FF2$cond)
pch1<- as.numeric(FF2$pch1)


#y<- y[x<1000]
#x<- x[x<1000]


#labels<-c( "all forests", expression(paste("-0.35"^"o","C cooling" )), "recurring dry summers", "50% less disturbances",  expression(paste("-0.24 mm day"^"-1","precip." )), "all forests w/ logging", "50% more disturbances",  "75% more disturbances", "100% more disturbances", "recurring wet summers", expression(paste("+0.37 mm day"^"-1","precip." )), "undisturbed forests", expression(paste("+0.33"^"o","C warming" )))


#palette(colorRampPalette(c("black","black","black","black", "grey60", "grey60", "grey60", "grey60"))( 100)) ## (n)
dev.new(width=4.25, height=3.7) #width=6, height=5)
plot(x,y, type="n",xlab=NA, ylab=NA, axes=FALSE) #, ylim=c(70,165), xlim=c(21.5, 24.4)) #main=c("Forest carbon equilibrium\n compared to tree morality size average over 5-years")

 
axis(side=1, tck=-0.01, labels=NA, lwd=0.75, at=c( 21.5, 22,22.5, 23,23.5, 24))
axis(side=2, tck=-0.01, labels=NA, lwd=0.75)
axis(side=1, lwd=0, line= -0.7, at=c( 21.5, 22,22.5, 23,23.5, 24))
axis(side=2, lwd=0, line= -0.7)

mtext(side=1, "Mean tree mortality diameter (cm)", line=1.2)
mtext(side=2, expression(paste("C steady-state (Mg ha"^"-1",")")), line=1.2) 

xoff3<-c(0.1, -0.18, -0.06, 0.18, 0.18, -0.18, 0.18, 0.18, 0.02, -0.06, 0.18, 0.18, 0.18) # rep(0.18,7))
yoff3<- c(10, -5, -9, 5, 6, 2, 2, 5,9, -9, 5,5,8)
segments(x, y, x+xoff3, y+yoff3) #, col="#00000090")
points(x,y, pch= 21, cex=1.2, bg= "black", col=col5[c(7,1,2,9,4,10,11,12,13,5,6,8,3)], lwd=2.5) #c(av, ct, de, ds, dt, lg, lg2, pe, pt, ud, wt))
segments(x-xerror, y, x+xerror, y, lwd=0.8, col="grey55")
segments(x, y-yerror, x, y+yerror, lwd=0.8, col="grey55")
abline(lm(y~x), col="grey5", lwd=0.95)
#points(x,y, pch=19)
points(x,y, pch= 19, cex=0.7, col= "black") #c(av, ct, de, ds, dt, lg, lg2, pe, pt, ud, wt))
labs3<- c("av", "ct", "de", "d0.5x", "dt", "lg", "d1.5x", "d2x", "d3x", "pe", "pt", "ud", "wt")
text(x+c(0.15,0, -0.1, 0,0,0,0,0,0.1,-0.1,0,0,0), y+c(1, -6, 0, 7, 8, 3, 3, 8, 0, 0, 6, 6, 8), labs3, pos=c(3,2,1,4,4,2,4,4,3,1,4,4,4), col="grey25", cex=0.98, offset=1)
summary(lm(y~x))










# # lab1<- c(12, 4, 2, 5, 3, 1)
# cols51<- c(8,9, 1, 4, 2, 7)
# legend(22.7, 170, c(labels[lab1]), pch=c(pch1[lab1]), text.col="white", col="black", cex=0.75, bty="n", pt.cex=0.95, pt.lwd=1.2)
# legend(22.7, 170, c(labels[lab1]), pch=c(pch1[lab1]), col=c(col5[cols51]), cex=0.75, bty="n", pt.cex=0.9)

# lab1<- c(13, 11, 10, 7, 8, 9, 6)
# cols51<- c(3, 6, 5, 11, 12, 13, 10)
# legend(21.3, 110, c(labels[lab1]), pch=c(pch1[lab1]), text.col="white", col="black", cex=0.75, bty="n", pt.cex=0.95, pt.lwd=1.2)
# legend(21.3, 110, c(labels[lab1]), pch=c(pch1[lab1]), col=c(col5[cols51]), cex=0.75, bty="n", pt.cex=0.9)
#box()
#####################################################################################
#THE END


# pe smaller tree mortallity shows that on average trees are larger???









