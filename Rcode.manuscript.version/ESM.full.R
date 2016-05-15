#EMPIRICAL SUCCESSION MAPPING Feb. 2015
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
 data.url <- "https://raw.githubusercontent.com/wanderswest/ESM-FIA/master/ESM.data.csv"
 
    ESM.data <- getURL(data.url)                
    ESM.data <- read.csv(textConnection(ESM.data), header = TRUE, sep = ",", quote="\"", dec=".")
RM5mergedfullcut1991.11.17.2014  <- as.data.table(ESM.data)


#or load filtered data from local drive
#RM5mergedfullcut1991.11.17.2014  <- as.data.table(read.csv("/Users/travis/GitHub/ESM-FIA/ESM.data.csv", header = TRUE, sep = ",", quote="\"", dec="."))

RM5mergedfullcut1991.11.17.2014 <-  subset(RM5mergedfullcut1991.11.17.2014, STDORGCD==0 ) # remove planted forests
S4swPREC <-  subset(RM5mergedfullcut1991.11.17.2014, STDORGCD==0 & cut==0 ) #exclude planted and harvested forests
F2 <- subset(S4swPREC)


#Intro staatistics
print(c(sum(RM5mergedfullcut1991.11.17.2014$startsumTPA)/(1/0.166166884), "Trees initially surveyed"))
print(c(sum(RM5mergedfullcut1991.11.17.2014$endsumTPA)/(1/0.166166884), "Trees resurveyed"))
print(c(sum(RM5mergedfullcut1991.11.17.2014$PLT_CN>0), "Plots resurveyed"))
print(c(mean(F2$REMPER), "mean resurvey period")) #most appropriate dataset for analyses 
print(c(std.error(F2$REMPER), "SE resurvey period"))
print(c(sum(RM5mergedfullcut1991.11.17.2014$cut>0), "Plots with harvesting"))
print(c(sum(RM5mergedfullcut1991.11.17.2014$STDORGCD>0, na.rm=T), "Planted plots"))

print(c(mean(F2$carbon1sum, na.rm=T)*0.00112085116, "mean forest carbon Mg/ha")) #0.00112085116 = convert ABG carbon to Mg/ha
print(c(sd(F2$carbon1sum, na.rm=T)*0.00112085116, "Standard devation forest carbon Mg/ha"))
LegacyCarbon <- ((F2$carbon1sum-F2$PREVcarbon1sum)/F2$REMPER)*0.00112085116	#0.00112085116 = convert ABG carbon to Mg/ha
print(c(mean(LegacyCarbon), "Legacy carbon accumulation Mg/ha/year")) 	
print(c(std.error(LegacyCarbon), "SE Legacy carbon accumulation Mg/ha/year")) 
NPPchg <- ((F2$NPP-F2$PREVNPP)/F2$REMPER)*0.00112085116	#0.00112085116 = convert ABG carbon to Mg/ha
print(c(mean(NPPchg), "NPP Mg/ha/year")) 	
print(c(std.error(NPPchg), "SE NPP Mg/ha/year")) 

###SUCCESSION STATS
#carbon and density variable name and unit conversions
F2$DIAmean <- (F2$DIAendmean* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean <- (F2$DIAbeginmean *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$Csum <- F2$carbon1sum* 0.00112085116  #	#convert living carbon to Mg/ha
F2$PREV_Csum <- F2$PREVcarbon1sum* 0.00112085116 ##convert living carbon to Mg/ha
F2$TPAsum <- (F2$endsumTPA /0.404686)			#number of trees per hectare by plot at resurvey
F2$PREV_TPAsum <- (F2$startsumTPA /0.404686)	#number of trees per hectare by plot at initial survey 
F2$PREV_STOCKINGmid <- F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey
#calculate change in Stem density and mean tree diameter for each plot
F2$TPAchg <- (F2$TPAsum-F2$PREV_TPAsum)/F2$REMPER
print(c(mean(F2$TPAchg, na.rm=T), " Overall TPAchg stems/ha/year")) 	
print(c(std.error(F2$TPAchg, na.rm=T), "SE TPA")) 

#maximum early successional ingrowth and carbon mean change vectors (FIG. 1)
STOCKval <- c(1, 9, 16, 25, 33, 42, 50, 59, 69, 79, 89,  101) 
QMDvalues <- c(12.7, 14.1, seq(15.5, 32, 1.5)) #,23.4, 24.5 ,26.8, 29, 32, 34.5, 36.5, 41 ) #10
s=4
h=1
S11r <- subset(F2, PREV_STOCKINGmid>= STOCKval[s] & PREV_STOCKINGmid <STOCKval[s+1]) 	
S11q <- subset(S11r, PREV_DIAmean>=(QMDvalues[h]) & PREV_DIAmean<(QMDvalues[h+1]))
print(c(mean(S11q$TPAchg, na.rm=T), "early succession TPAchg stems/ha/year")) 	
print(c(std.error(S11q$TPAchg, na.rm=T), "SE early succession")) 

#calculate NPP for each plot
F2$NPPchg <- ((F2$NPP-F2$PREVNPP)/F2$REMPER)*0.00112085116	#0.00112085116 = convert ABG carbon to Mg/ha
STOCKval <- c(0, 5, 15, seq(30, 121, 91/6)) #6
QMDvalues <- c(12.7, 14.5) 
s=4
h=1
S11r <- subset(F2, PREV_STOCKINGmid>= STOCKval[s] & PREV_STOCKINGmid <STOCKval[s+1]) 	
S11q <- subset(S11r, PREV_DIAmean>=(QMDvalues[h]) & PREV_DIAmean<(QMDvalues[h+1]))
print(c(round(mean(S11q$NPPchg, na.rm=T),1), "early succession NPP Mg/ha/year")) 	
print(c(std.error(S11q$NPPchg, na.rm=T), "SE early succession NPP")) 


#reforestation ingrowth
F2e <- subset(F2, PREV_DIAmean>18.5 & PREV_STOCKINGmid< 20) 
print(c(mean(F2e$NPP, na.rm=T), "early succession NPP Mg/ha/year")) 	
print(c(std.error(F2e$NPP, na.rm=T), "SE early succession NPP")) 
print(c(mean(F2e$TPAchg, na.rm=T), "early succession TPAchg stems/ha/year")) 	
print(c(std.error(F2e$TPAchg, na.rm=T), "SE early succession NPP")) 

print(c(sum(F2e$ingrowth2>0)/sum(F2e$PLT_CN>0)*100, "percent ingrowth in low density plots"))
#print(c(sum(F2e$mortality2>0 & F2e$mortDIAmean>F2e$DIAendmean, na.rm=T)/sum(F2e$PLT_CN>0)*100, "percent mortality large trees in low density plots"))

#### Largest mean diameter forest stats
F2e <- subset(F2,PREV_DIAmean>35 ) #& PREV_STOCKINGmid>60)
print(c(sum(F2e$PLT_CN>0)/sum(F2$PLT_CN>0)*100, "largest mean diameter plots percent of overall dataset"))


####DISTURBANCE STATS
F2$allmortalitynoAD <- (F2$unknowndamage+F2$vegetation) 
F2$alldisturb <- (F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)
F2$allmortality <- (F2$insect+F2$fire+F2$weather+F2$disease+F2$animal+F2$unknowndamage+F2$vegetation)

print(c(sum(F2$DSTRBCD1>0)/sum(F2$PLT_CN>0)*100, "Percent dataset large disturbance in 5 years" ) )
print(c((1-(sum(F2$allmortality>0)/sum(F2$PLT_CN>0)))*100, "Percent dataset without mortality in 5 years" ) )

print(c(sum(F2$alldisturb>0)/sum(F2$PLT_CN>0)*100, "Percent dataset disturbed in 5 years" ) )
print(c(sum(F2$unknowndamage>0 & F2$alldisturb==0 & F2$vegetation==0)/sum(F2$PLT_CN>0)*100, "Percent dataset only unknown disturbance in 5 years" ) )
print(c(sum(F2$unknowndamage)/(sum(F2$alldisturb)+ sum(F2$allmortalitynoAD))*100, "Percent mortality attributed to unknown disturbance in 5 years" ) )
print(c((sum(F2$alldisturb)/sum(F2$allmortality))*100, "Percent mortality attributed to unknown disturbance in 5 years" )) 

F2D <- subset(F2, PREV_STOCKINGmid>60)	#Disturbed forests
print(c(sum(F2D$alldisturb>0)/sum(F2D$PLT_CN>0)*100, "Percent self thinning forests disturbed in 5 years" ) )



############################################################################################################################
#######################   FIGURE 1 - EMPIRICAL FOREST DENSITY MODEL  #######################################################
####

dev.new(width=8, height=6) #set up plot
S4swPREC <-  subset(RM5mergedfullcut1991.11.17.2014, STDORGCD==0 & cut==0 ) #exclude planted and logged forests
F2 <- subset(S4swPREC) 

#Convert variables to standard names and metric units
F2$DIAmean <- (F2$DIAendmean* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean <- (F2$DIAbeginmean *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum <- (F2$endsumTPA /0.404686)			#number of trees per hectare by plot at resurvey
F2$PREV_TPAsum <- (F2$startsumTPA /0.404686)	#number of trees per hectare by plot at initial survey (cm)
F2$PREV_STOCKINGmid <- F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey

#calculate change in Stem density and mean tree diameter for each plot
F2$TPAchg <- F2$TPAsum-F2$PREV_TPAsum
F2$DIAchg <- F2$DIAmean-F2$PREV_DIAmean

#setup color ramp
 palette(adjustcolor(((rich.colors(60))), transform=diag(c(0.95,1,0.85,1)))) 
stockcols <- rev(palette())
#Climate Colors:
pt <- "steelblue1" #pluvial trend
pe <- "royalblue"	#pluvial event
dt <- "goldenrod1" #drought trend
de <- "orange2"  	#drought event
nt <- "grey55"
ne <- "grey35"
ag1 <- 25


#select model type and bootstraped error on/off
ModelType <- 3 #1=Even increment w/ moving average, 2=Even increment only, 3=Textbook 
ErrorOn <- 1  #1=Yes, 0=No

#create plot
plot(c((F2$PREV_DIAmean),(F2$DIAmean)),c((F2$PREV_TPAsum),(F2$TPAsum)), type ="n", col="grey94", 
	main=c("Empirical Model of 5-year forest density change, DBH>=12.7 cm", "FIA survey 1998-2012, Eastern USA") 
	,ylim=c(50,950), xlim=c(14, 43), mgp=c(2,1,0)
	,xlab=NA, ylab=NA, axes=FALSE) #1250

axis(side=1, tck=-0.01, labels=NA, lwd=0.75)
axis(side=2, tck=-0.01, labels=NA, lwd=0.75)
axis(side=1, lwd=0, line= -0.7)
axis(side=2, lwd=0, line= -0.7)
mtext(side=1, "Mean tree diameter (cm)", line=1.2)
mtext(side=2, expression(paste("Stem density (trees ha"^"-1",")")), line=1.2)
box()

#plot all plots as individual background vectors colored by inital stocking 
arrows((F2$PREV_DIAmean),(F2$PREV_TPAsum),(F2$DIAmean),(F2$TPAsum), col= (F2$PREV_STOCKINGmid/4)+18, length=0.065, angle=22,lwd=0.9)  

#Set up plot with relative plot locations and legend - for dev.new(width=8, height=6) 
x1 <- 39*1.05
x2 <- 39.8*1.05
y1 <- (c(seq(1200,820,-38))-50)*0.73#-22
y2 <- (c(seq(1220,800,-5))-75)*.75 #-22
x1b <- (x1-5.5)*0.965 #x1b <- x1-7.5
y1b <- (y1[1]+110)*.69 #-264 #y1b <- y1[1]-215
cex1 <- 0.8
cexA <- 3 #succession large labels
t3 <- 32 #33.5  #transition mean diameter between independent vectors and moving average
lw2 <- 2 #arrow width
slopetable <- NULL
colS <- "#00000090"

#legend
segments(x1-0.2,y2,x2+0.2, y2, col=c(stockcols[seq(16,46,0.35)]),lwd=2, lend=2)
points(30,750, pch=16, cex=200, col="#FFFFFF75" ) #FFFFFF50
points(x1b-0.5,y1b+165, pch=15, cex=25, col="#FFFFFF90" ) #FFFFFF75
#segments(t3,0,t3,600, col="white", lwd=0.75)
segments(c(x1-1.5,x1-1.5,x1-1.5),c( y2[1]+5, y2[39], y2[85]),x2+0.3,c(y2[1]+5, y2[39], y2[85]), col="black",lwd=0.9, lend=2)
segments((x1-1.5),( y2[10]),x2+0.3,(y2[10]), col="black",lwd=0.9, lend=2, lty="dashed")
arrows(x1b-2.5, y1b, x1b+1.5, y1b, code=3, length=0.08)
arrows(x1b, y1b+200,x1b, y1b-75, code=3, length=0.08)
arrows(x1,y1,x2,y1, col="black", length=0.085, angle=25,lwd=lw2+0.75)  
arrows(x1,y1,x2,y1, col=c(stockcols[seq(17,44,2.4)]), length=0.08, angle=28,lwd=lw2)  
#text(x1+2.25, y1[1]+55,"Initial relative density", cex=0.9, pos=2) #initial plot stocking and mean change vector increments
text(x1, y1[1]-5,"Over", cex=cex1, pos=2) #initial plot stocking and mean change vector increments
text(x1, y1[4]+10,"Full", cex=cex1, pos=2) #initial plot stocking and mean change vector increments
text(x1, y1[9]+10,"Under", cex=cex1, pos=2) #initial plot stocking and mean change vector increments

F2l <- subset(F2, PREVSTOCKING5mid>60  ) #
fit4<-rq(log10(TPAsum) ~ log10(DIAmean), tau=0.90, data=F2l)
lml<-summary(fit4)
x<-seq(0,50, 0.5)
y=((x^ lml$coefficients[2])* 10^lml$coefficients[1]) # #325000  #Best fit for early succession without saplings
lines((y)~(x), col="#00000090", lty="dotted", lwd=6.5) # "#00000090"

#Future average:	 BASED ON results from future equilibrium model below
cond <- "3"
x <- 34.27
y <- 280.8
xer <- 0.372*2
yer <- 5.096*2
points(x,y, pch=21, col=NA, bg="grey35", cex=0.75)
segments(x-xer,y, x+xer, y, col="#00000075", lwd=1.5)
segments(x,y-yer, x, y+yer, col="#00000075", lwd=1.5)

#manually remove overlapping vectors
DuplicatedVectors <- c(31, 31.5,33.5, 36, 38, 8026.2, 25012.7, 50014, 9027.5, 1024.5, 1026, 1027.5, 16032, 16012.7, 1029, 1032, 1030.5, 1612.7, 9030.5, 9032, 15032, 9029, 50014.1, 10031.5, 10036, 9024.5, 42030.5, 16030.5, 9531.5, 9534.5, 9536, 12033.5, 12036, 33029, 40, 9540, 12040, 80040, 80036, 80038) #16029  9633.5, 9636, 9532

##ESM mean vector calculation and plotting with for loops
L <- ifelse(ModelType==3, 4,1)
for(o in ErrorOn:2){
for(l in 1:L){

####Organize and subset data to plot with FOR loops 
if(ModelType==1){
###Null model even increments, moving average vectors
STOCKval <- c(0,0,seq(0, 141, 8))
QMDvalues <- seq(0,50, 0.5)
s1 <- s2 <- s3 <- s4 <- 20
increment <- 3
}

if(ModelType==2){
###Null model even increment vectors
STOCKval <- seq(0, 121, 10)
QMDvalues <- seq(0,50, 2)
s1 <- s2 <- s3 <- s4 <- 20
increment <- 1
}

if(ModelType==3){
#ESM Density with non-overlapping vectors
#looped by phase space location to have sufficient underlying data and to exclude excessive vectors

if(l==1){	#independent dense
STOCKval <- c(1, 9, 16, 25, 33, 42, 50, 59, 69, 79, 89,  101) 
QMDvalues <- c(12.7, 14.1, seq(15.5, t3, 1.5)) 
s1 <- s2 <- s3 <- s4 <- 3
increment <- 1}

if(l==4){	  ##moving average end
STOCKval <- c(0, 9.5, 12, 25, NA, NA, 12, 20, 35, 50, 65, 80, 100, 120, 140 ) 
#QMDvalues <- c(seq(t3, 42, 1.8), 42)
QMDvalues <- c(t3-0.5, 33.5, 34.5, 36, 40, 42, NA, NA, 38, 40, 45, 47)
#QMDvalues <- c(t3, 33.5, 35, 38, 42, 45)
s1 <- 3
s2 <- 3
s3 <- 3
s4 <- 3
increment <- 2}

if(l==3){	  #overstocked
STOCKval <- c(100, 122)
QMDvalues <- c(12.5, 15, 18.5, 20, 21.5, 23, 24.5, 26.8, 28.5,  33)
s1 <- s2 <- s3 <- s4 <- 3
increment <- 1}

if(l==2){	    #Bottom right###
STOCKval <- c(0, 8, 16)  
QMDvalues <- c(26.2, 28, 31, 34) #, 38, 44)
s1 <- s2 <- s3 <- s4 <- 4 
increment <- 1}
} #for Model type 3
	
#Subset by stocking range then mean diameter range
if(o==1)
{
g <- length(QMDvalues-1)
r <- length(STOCKval-1)
 for(s in 1:r)
 { S11r <- subset(F2, PREV_STOCKINGmid>= STOCKval[s] & PREV_STOCKINGmid <STOCKval[s+increment]) 	
	for(h in 1:g)
	{ 
	S11q <- subset(S11r, PREV_DIAmean>=(QMDvalues[h]) & PREV_DIAmean<(QMDvalues[h+increment]))
	if(length(S11q$PREV_DIAmean)>s1 & s==1 | length(S11q$PREV_DIAmean)>s2 & s==2 | length(S11q$PREV_DIAmean)>s3  & s==3 | length(S11q$PREV_DIAmean)>s4  & s>=4)
		{		
		stockcol <- (mean(S11q$PREV_STOCKINGmid,na.rm=TRUE)/4)+18
		stkqmd <- STOCKval[s]*1000+QMDvalues[h]
if(stkqmd %in% DuplicatedVectors){next} #remove overlapping vectors
if(mean(S11q$TPAsum, na.rm=TRUE)>1000){next} #exclude outliers

#Plot mean change vector
arrows((mean(S11q$PREV_DIAmean, na.rm=TRUE)),(mean(S11q$PREV_TPAsum, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE)), col=stockcol, length=0.08, angle=28,lwd=lw2+1.25)  

#Plot bootstraped error associated with mean diameter and stem density change 
endTPAerror <- boot(S11q$TPAchg, meanFunc, R=1000)	
endTPAerror <- sqrt(var(endTPAerror$t))*1.96
			lines(c((mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE))),
				 c((mean(S11q$TPAsum, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE))+endTPAerror),  
					col="grey60", type="l", lwd=0.7)
			lines(c((mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE))),
				 c((mean(S11q$TPAsum, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE))-endTPAerror),  
					col="grey60", type="l", lwd=0.7)
endQMDerror <- boot(S11q$DIAchg, meanFunc, R=1000)				
endQMDerror <- sqrt(var(endQMDerror$t))*1.96
			lines(c((mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)+endQMDerror)),
				 c((mean(S11q$TPAsum, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE))),  
					col="grey60", type="l", lwd=0.7)
			lines(c((mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)-endQMDerror)),
				 c((mean(S11q$TPAsum, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE))),  
					col="grey60", type="l", lwd=0.7)
		}
	}
  }
}



#Replot with bold vectors
if(o==2)
{
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
#print(stkqmd)
if(stkqmd %in% DuplicatedVectors){next} #remove overlapping vectors
if(mean(S11q$TPAsum, na.rm=TRUE)>1000){next} #exclude outliers

stockcol <- (mean(S11q$PREV_STOCKINGmid,na.rm=TRUE)/4)+18 

#summery change inset - summary of mean vectors
 	arrows(x1b,y1b,(x1b+mean(S11q$DIAchg, na.rm=TRUE)),(y1b+mean(S11q$TPAchg, na.rm=TRUE)), col="grey15", length=0.055, angle=25,lwd=lw2+0.4) 
	arrows(x1b,y1b,(x1b+mean(S11q$DIAchg, na.rm=TRUE)),(y1b+mean(S11q$TPAchg, na.rm=TRUE)), col=stockcol, length=0.05, angle=28,lwd=lw2)
	slope <- mean(S11q$TPAchg, na.rm=TRUE)/mean(S11q$DIAchg, na.rm=TRUE)
	slopetable <- rbind(slopetable, slope)

arrows((mean(S11q$PREV_DIAmean, na.rm=TRUE)),(mean(S11q$PREV_TPAsum, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE)), col= "grey15", length=0.08, angle=25,lwd=lw2+1)  
arrows((mean(S11q$PREV_DIAmean, na.rm=TRUE)),(mean(S11q$PREV_TPAsum, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE)), col= stockcol, length=0.08, angle=28,lwd=lw2)	

#text(mean(S11q$DIAmean, na.rm=TRUE),mean(S11q$TPAsum, na.rm=TRUE), label= stkqmd)
		}
  	}
  }
}
}  #l loop
} #o loop

summary(slopetable)
#segments(x1-6,y1[1]+115,x1+2, y1[1]+100, col="#FFFFFF95",lwd=23, lend=2)
boxed.labels(x1-5, y1[1]+110, "5-year Temperate Forest Density Change",bg="white", border=F, cex=1.05, col="black", font=3, ypad=0.7)

#####Plot summary INSET text
arrows(x1b-2.5, y1b, x1b+1.5, y1b, code=3, length=0.08, col="#00000050")
arrows(x1b, y1b+200,x1b, y1b-75, code=3, length=0.08, col="#00000050")
text(x1b+1.7, y1b-45,"Self\nthinning", cex=cex1) #(letters) could correspond to % of each in each quadrant (A%|B%|C%)
text(x1b+2.2, y1b,"Gap filling", pos=3, cex=cex1)
par(lheight=.75)
text(x1b-2.5, y1b+50, "Low density\noverstory\nrecruitment", pos=3, cex=cex1) 
text(x1b+0.5, y1b+140,"Early\nsuccession", pos=4, cex=cex1)
text(x1b-0.1, y1b+180,"+200 trees", pos=3, cex=cex1-0.1)
text(x1b-0.1, y1b-55,"-75 trees", pos=1, cex=cex1-0.1)
text(x1b+1.1, y1b,"+1.5 cm", pos=4, cex=cex1-0.1)
text(x1b-2.2, y1b,"-2.5 cm", pos=2, cex=cex1-0.1)
#text(x1b-1.5, y1b-60,"Equilibrium", pos=3, cex=cex1, col="black")
#filledellipse(rx1=0.75, ry1=20, mid=c(x1b-0.3,y1b), angle=1.5, lcol="black", col=NULL, lty=3, lwd=0.5)
#text(x1b, y1b+265,(("Successional trajectory summary")), pos=3, cex=cex1)
#boxed.labels(x1b-0.65, y1[1]+55,"Successional summary",bg="#FFFFFF50", border=F, cex=0.9)
#boxed.labels(x1b-5.25, y1[1]+35,"Successional\nsummary",bg="#FFFFFF00", border=F, cex=0.9, pos=4)
boxed.labels(x1b-3, y1[1]+52,"Inset A",bg="#FFFFFF50", border=F, cex=0.9, xpad=1.3, ypad=1.3)
lines(c(x1b-4.75, x1b-4.75), c(y1b+10, y1b+260), lwd=0.5)
#lines(c(x1b-4, x1b-2.05), c(y1b-55, y1b-55), lwd=0.5)
text(x1+0.83, c(y1+10, y1[11]-15),c("120","100","90", "80","70", "60","50","40","30", "20", "10", "1"), cex=0.85, pos=4)
text(x1+2.7, y1[1]+52,"Initial relative density", cex=0.9, pos=2) #initial plot stocking and mean change vector increments
box()
	
	

############### END DENSITY MODEL #######################################################################################







#########################################################################################################################
###############FIGURE 1B - ESM CARBON MODEL #############################################################################

F2 <- subset(S4swPREC)
#Re-use previous plot variables - Change Y-axis to vegetation carbon
F2$DIAmean <- (F2$DIAendmean* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean <- (F2$DIAbeginmean *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum <- F2$carbon1sum* 0.00112085116  #	#convert living carbon to Mg/ha
F2$PREV_TPAsum <- F2$PREVcarbon1sum* 0.00112085116 ##convert living carbon to Mg/ha
F2$PREV_STOCKINGmid <- F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey

#calculate change in CARBON and mean tree diameter for each plot
F2$TPAchg <- F2$TPAsum-F2$PREV_TPAsum
F2$DIAchg <- F2$DIAmean-F2$PREV_DIAmean
######################################################################
#set up plot
dev.new(width=8, height=4.5)

#create plot
plot(c((F2$PREV_DIAmean),(F2$DIAmean)),c((F2$PREV_TPAsum),(F2$TPAsum)), type ="n", col="grey94", 
	main=c("Empirical Model of 5-year forest carbon change, DBH>=12.7 cm", "FIA survey 1997-2012, Eastern USA") 
	,ylim=c(6.5,200), xlim=c(14, 43), mgp=c(2,1,0)
	,xlab=NA, ylab=NA, axes=FALSE) #1250

axis(side=1, tck=-0.01, labels=NA, lwd=0.75)
axis(side=2, tck=-0.01, labels=NA, lwd=0.75)
axis(side=1, lwd=0, line= -0.7)
axis(side=2, lwd=0, line= -0.7)

mtext(side=1, "Mean tree diameter (cm)", line=1.2)
mtext(side=2, expression(paste("Carbon, living trees (Mg ha"^"-1",")")), line=1.2)

box()

#plot all plots as individual background vectors colored by inital stocking 
arrows((F2$PREV_DIAmean),(F2$PREV_TPAsum),(F2$DIAmean),(F2$TPAsum), col= (F2$PREV_STOCKINGmid/4)+18, length=0.065, angle=22,lwd=0.9)
  
#Set up plot with relative plot locations and legend - ideal for dev.new(width=8, height=4.5)
x1 <- 14-0.5
x2 <- 14+0.5
y1 <- (c(seq(32,26.05,-.55))*2.5-52)*6.1
y2 <- (c(seq(32.2,26.25,-0.085))*2.5-52)*6.1
y3 <- (c(seq(31.8,26,-0.73))*2.5-52)*6.1
x1b <- x1-7 #x1b <- x1-7.5
y1b <- y1[1]-264 #y1b <- y1[1]-215
cex1 <- 0.8
t3 <- 32 #transition mean diameter between independent vectors and moving average
lw2 <- 2 #arrow width

#LEGEND
stockcols2 <- rev((seq(1,120,15)/4)+18)
segments(x1,y2,x2, y2, c(seq(45,18, -0.384)),lwd=2, lend=2) #c(seq(66,38,-0.3))  col="grey85"
points(30,100, pch=16, cex=200, col="#FFFFFF75" ) #FFFFFF75
text(x1-0.25, y1[1]+26,"5-year Temperate Forest Carbon Change", pos=4, cex=1.05, font=3)
segments(c(x1,x1,x1),c( y2[1], y2[36],  y2[71]),x2+1.25 ,c(y2[1], y2[36],  y2[71]), col="black",lwd=0.9, lend=2)
segments(c(x1),c(y2[10]),x2+1.25,c(y2[10]), col="black",lwd=0.9, lend=2, lty="dashed")
arrows(x1+0.1,y3,x2-0.1,y3, col="black", length=0.085, angle=25,lwd=lw2+1)  
arrows(x1+0.1,y3,x2-0.1,y3, col=c(stockcols2), length=0.08, angle=28,lwd=lw2)  
text(x1-.5, y1[1]+10,"Initial relative density", cex=0.9, pos=4) #initial plot stocking and mean change vector increments
text(x2, y1[1]-5,"Over", cex=cex1, pos=4) #initial plot stocking and mean change vector increments
text(x2, y1[4],"Full", cex=cex1, pos=4) #initial plot stocking and mean change vector increments
text(x2, y1[9],"Under", cex=cex1, pos=4) #initial plot stocking and mean change vector increments

#Carbon thinning line - fit to tree diameters >5 inches with relative densities between 90-110
F2l <- subset(F2, PREVSTOCKING5mid>=90 & PREVSTOCKING5mid<=110 & DIAmean<40)
lml <- summary(lm((F2l$TPAsum)~(F2l$DIAmean) ))
print(c(lml$coefficients[2], "self thinning line best fit scaling exponent"))
x=seq(0,50,0.5) #seq(22,35.5,0.5)
y=(lml$coefficients[2]*x+lml$coefficients[1])  #Best fit 
lines((y)~(x), col="#00000090", lty="dotted", lwd=5) # "#00000085"


#Future average: BASED ON results from future equilibrium model below
cond <- "3"	
x <- 34.27
y <- 116.5
xer <- 0.372*2
yer <- 3.54*2
points(x,y, pch=21, col=NA, bg="grey35", cex=0.75)
segments(x-xer,y, x+xer, y, col="#00000075", lwd=1.5)
segments(x,y-yer, x, y+yer, col="#00000075", lwd=1.5)


####Organize and subset data to plot with FOR loops 
#looped by phase space location to have sufficient underlying data and to exclude excessive vectors
DuplicatedVectors <- c(5031, 5034.6, 5036.4, 15036.4) 
slopetable <- NULL
for(o in 1:2){
for(l in 1:5){

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
STOCKval <- c(5, 15, NA, 30, 38, NA, NA,20, 40, 55, 55, 60, 75, NA, NA,NA, 65, 85, 100, 85, 98, 121) #121, 121) #121,121,121
QMDvalues <- c(seq(t3-1,45,1.8), 50)  #
s1 <- 40 #100
s2 <- 30 #50
s3 <- 3#40
s4 <- 5
increment=3}

if(l==5){	  ##Early early succession
STOCKval <- c(0, 5, 15, seq(30, 121, 91/6)) #6
QMDvalues <- c(12.7, 14.5) 
s1 <- 3 #30
s2 <- 3 #10
s3 <- 3#10
s4 <- 3#5
increment=1}

#Subset by stocking range then mean diameter range
if(o==1){
g <- length(QMDvalues-1)
r <- length(STOCKval-1)
 for(s in 1:r)
 { S11r <- subset(F2, PREV_STOCKINGmid>= STOCKval[s] & PREV_STOCKINGmid <STOCKval[s+increment]) 	
	for(h in 1:g)
	{ 
	S11q <- subset(S11r, PREV_DIAmean>=(QMDvalues[h]) & PREV_DIAmean<(QMDvalues[h+increment]))
	if(length(S11q$PREV_DIAmean)>s1 & s==1 | length(S11q$PREV_DIAmean)>s2 & s==2 | length(S11q$PREV_DIAmean)>s3  & s==3 | length(S11q$PREV_DIAmean)>s4  & s>=4)
		{		
		stkqmd <- STOCKval[s]*1000+QMDvalues[h]
		if(stkqmd %in% DuplicatedVectors){next} #remove overlapping vectors
		stockcol <- (mean(S11q$PREV_STOCKINGmid,na.rm=TRUE)/4)+18
#Plot mean change vector
		arrows((mean(S11q$PREV_DIAmean, na.rm=TRUE)),(mean(S11q$PREV_TPAsum, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE)),col=(stockcol), length=0.085, angle=25,lwd=lw2+1.25) 

#Plot bootstraped error associated with mean diameter and stem density change 
endTPAerror <- boot(S11q$TPAchg, meanFunc, R=1000)	
endTPAerror <- sqrt(var(endTPAerror$t))*1.96
			lines(c((mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE))),
				 c((mean(S11q$TPAsum, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE))+endTPAerror),  
					col="grey50", type="l", lwd=0.75)
			lines(c((mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE))),
				 c((mean(S11q$TPAsum, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE))-endTPAerror),  
					col="grey50", type="l", lwd=0.75)
endQMDerror <- boot(S11q$DIAchg, meanFunc, R=1000)				
endQMDerror <- sqrt(var(endQMDerror$t))*1.96
			lines(c((mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)+endQMDerror)),
				 c((mean(S11q$TPAsum, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE))),  
					col="grey50", type="l", lwd=0.75)
			lines(c((mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)-endQMDerror)),
				 c((mean(S11q$TPAsum, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE))),  
					col="grey50", type="l", lwd=0.75)
	}}}}

#Replot with bold vectors
if(o==2){
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
	arrows((mean(S11q$PREV_DIAmean, na.rm=TRUE)),(mean(S11q$PREV_TPAsum, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE)),col=("grey15"), length=0.085, angle=25,lwd=lw2+0.75) 
	arrows((mean(S11q$PREV_DIAmean, na.rm=TRUE)),(mean(S11q$PREV_TPAsum, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE)), col= (S11q$PREV_STOCKINGmid/4)+18, length=0.08, angle=28,lwd=lw2)  # VectCols[s] col=s*1.8+10,  s*1.2+7 # s*2+10 #stockcol

print(c(stkqmd, l, s, h, mean(S11q$TPAchg, na.rm=T), sum(S11q$TPAchg>0, na.rm=T)))

	}}}}
 } #l loop
} #o loop
box()


####################### FIGURE 1 - END  #################################################################
####
#







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
F2 <- subset(S4swPREC, temp5grow>0 ) #& SISP<300) #remove any NA's from temperature data

#dataset for undisturbed analysis Fig S2
if(1==2)
{
F2$vegetation[F2$vegetation!=0]<-F2$vegetation[F2$vegetation>0]+rnorm(sum(F2$vegetation>0),0, 0.1)
vg<-quantile(F2$vegetation[F2$vegetation!=0], 0.99, na.rm=T)
F2$unknowndamage[F2$unknowndamage!=0]<-F2$unknowndamage[F2$unknowndamage>0]+rnorm(sum(F2$unknowndamage>0),0, 0.1)
un<-quantile(F2$unknowndamage[F2$unknowndamage!=0], 0.98, na.rm=T)
F2 <- subset(F2,  unknowndamage<=un & vegetation<=vg & disease==0 & animal==0 & insect==0 & weather==0 & fire==0) 
}

#Convert variables to standard names and metric units
F2$DIAmean <- (F2$DIAendmean* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean <- (F2$DIAbeginmean *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum <- (F2$endsumTPA /0.404686)			#number of trees per hectare by plot at resurvey
F2$PREV_TPAsum <- (F2$startsumTPA /0.404686)	#number of trees per hectare by plot at initial survey (cm)
F2$PREV_STOCKINGmid <- F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey
F2$TPAchg <- (F2$TPAsum-F2$PREV_TPAsum)/F2$REMPER
F2$SISPchg<-ifelse(F2$SISP>=300, 1, ifelse(F2$SISP<300, -1, NA) )
F2$PREVcarbon1sum<-F2$PREVNPP
F2$carbon1sum<-F2$NPP


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

vari <- c("PREC5growchg", "TEMP5growchg", "SWpluvial5chg", "SWdrought5chg") #"PREC5growchg", "TEMP5growchg"
vari1 <- subset(F2, select=c(PREC5growchg, TEMP5growchg, SWpluvial5chg, SWdrought5chg))
vari2 <- c("5-yr JJA precip. anomaly (mm/day)", "5-yr JJA temp. anomaly (C)", "1-yr JJA soil moist. anomaly (z-score)", "1-yr JJA soil moist. anomaly (z-score)")

qt1 <- c(seq(0,1,1/9)/100,seq(1,59,2)/100)
qt2 <- c(seq(41,99, 2)/100, seq(99,100,1/9)/100) #0.98, 0.99 added to ease moving average transistion to 1 at end
int <- 10
F5 <- NULL
for(v in 1:3)
{v=v
F4<-foreach(c = 1:60, .combine=rbind) %dopar%
{
	if(c<=30 & v<3)
	{
		f22 <- subset(F2, get(vari[v])> quantile(vari1[[v]], qt2[c], na.rm=TRUE) & get(vari[v])< quantile(vari1[[v]], qt2[c+int], na.rm=TRUE)) 		
		f23 <- subset(F2, get(vari[v])>= quantile(vari1[[v]], 0.3, na.rm=TRUE) & get(vari[v])<= quantile(vari1[[v]], 0.7, na.rm=TRUE), select=c(PLT_CN)) 
		ccond <- vari[v]
		col3=ifelse(qt2[c+(int/2)]>=0.85, 3, 2)
		qta <- qt2[c]
		qtb <- qt2[c+int]
		}
	if(c>30 & v<3)
	{ 
		d=c-20
		f22 <- subset(F2, get(vari[v])< quantile(vari1[[v]], qt1[d], na.rm=TRUE) & get(vari[v])> quantile(vari1[[v]], qt1[d-int], na.rm=TRUE)) #		
		f23 <- subset(F2, get(vari[v])>= quantile(vari1[[v]], 0.3, na.rm=TRUE) & get(vari[v])<= quantile(vari1[[v]], 0.7, na.rm=TRUE), select=c(PLT_CN)) 
		ccond <- vari[v]
		col3=ifelse(qt1[d-(int/2)]<=0.15, 1, 2)
		qtb <- qt1[d]
		qta <- qt1[d-int]
	}
		
	
	if(c<=30 & v==3)
	{
			f22 <- subset(F2, get(vari[v])> quantile(vari1[[v]], qt2[c], na.rm=TRUE) & get(vari[v])< quantile(vari1[[v]], qt2[c+int], na.rm=TRUE) & SWdrought5chg>= quantile(F2$SWdrought5chg, 0.15, na.rm=TRUE)) #exclude droughts
			f23 <- subset(F2, get(vari[v])<= quantile(vari1[[v]], 0.7, na.rm=TRUE) & get(vari[v])>= quantile(vari1[[v]], 0.3, na.rm=TRUE) & SWdrought5chg>= quantile(F2$SWdrought5chg, 0.15, na.rm=TRUE), select=c(PLT_CN)) 
		ccond <- "Wet summer"
		col3=ifelse(qt2[c+(int/2)]>=0.85, 3, 2)
		qta <- qt2[c]
		qtb <- qt2[c+int]
	}
	if(c>30 & v==3)
	{ 
		d=c-20
		f22 <- subset(F2, get(vari[v+1])< quantile(vari1[[v+1]], qt1[d], na.rm=TRUE) & get(vari[v+1])> quantile(vari1[[v+1]], qt1[d-int], na.rm=TRUE) & SWpluvial5chg<= quantile(F2$SWpluvial5chg, 0.85, na.rm=TRUE)) #exclude wet summers  
		f23 <- subset(F2, get(vari[v+1])<= quantile(vari1[[v+1]], 0.7, na.rm=TRUE) & get(vari[v+1])>= quantile(vari1[[v+1]], 0.3, na.rm=TRUE) & SWpluvial5chg<= quantile(F2$SWpluvial5chg, 0.85, na.rm=TRUE), select=c(PLT_CN)) 
		ccond <- "Dry summer"
		col3=ifelse(qt1[d-(int/2)]<=0.15, 1, 2)
		qtb <- qt1[d]
		qta <- qt1[d-int]
	}
	


f22$nsamp <- (ave(f22$DIAmean>0,f22$LATLON, FUN=function(x) sum(x, na.rm=TRUE)))/f22$nsampTOT #sample proportion of each grid cell 
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

F33 <- subset(F33, aPREVSTOCKDIAbinlength> 25) #Null model needs at least 25 comparison plots at each diameter/stocking bin

##Carbon Model
c6 <- (((F33$carbon1sum-F33$PREVcarbon1sum)* 0.00112085116 )/F33$REMPER)
F33$avgCchg <- (ave(c6*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE)))/F33$nbin #multiply by AC proportion divided by mean proportion in each stocking/dia bin

##Carbon Model 5 inch and greater trees
c5 <- (((F33$carbon5sum-F33$PREVcarbon5sum)* 0.00112085116 )/F33$REMPER)
F33$avgC5chg <- (ave(c5*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE)))/F33$nbin

##radial growth 5,10,15 inch min DIA 
rx5 <- (F33$DIAgrowmean/F33$REMPER)/2 ## convert from diameter to radial growth
F33$avgx5 <- ave(rx5*F33$nsamp,F33$aPREVSTOCKDIAbin, FUN=function(x) mean(x, na.rm=TRUE))/F33$nbin

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
F34 <- subset(F33, !duplicated(aPREVSTOCKDIAbin), select=c(aPREVSTOCKDIAbin, avgCchg, avgx5, avgTPAchg, sapavgTPAchg, avgC5chg, avgmortx5, sapavgMORT, sapavgRECR))
f22$aPREVSTOCKDIAbin <- floor(f22$PREV_STOCKINGmid/10)*10000+(floor(f22$PREV_DIAmean/2)*2)
f22 <- merge(f22, F34, by="aPREVSTOCKDIAbin") #, all.x=T

#carbon 1 and 5 inch
c8 <- (((f22$carbon1sum-f22$PREVcarbon1sum)* 0.00112085116 )/f22$REMPER)
f22$xcarbon <- ((c8)-f22$avgCchg)

c55 <- (((f22$carbon5sum-f22$PREVcarbon5sum)* 0.00112085116 )/f22$REMPER)
f22$xcarbon5 <- ((c55)-f22$avgC5chg)

##radial growth  
x8 <- (f22$DIAgrowmean/f22$REMPER)/2  # convert from diameter to radial growth
f22$x5 <-  (x8-f22$avgx5)

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
sp22 <- mean(f22$SISPchg-f22$sapavgTPAchg, na.rm=T)

#calculate mean condition for x-axis based on v loop
if(length(f22$xcarbon)>= 10)
{
if(v==2)
{
	x3 <- mean(f22$TEMP5growchg, na.rm=TRUE)
}
if(v==1)
{
	x3 <- mean(f22$PREC5growchg, na.rm=TRUE)
}
if(v==3)
{
x3 <- ifelse(c<=30,mean(f22$SWpluvial5chg, na.rm=TRUE), mean(f22$SWdrought5chg, na.rm=TRUE))
}

#calculate y-axis variables
y3 <- mean(f22$xcarbon, na.rm=TRUE) 
y3pctdif<-(mean(f22$xcarbon, na.rm=TRUE)/mean(f22$avgCchg, na.rm=TRUE))*100   #
y3start <-  mean(c6, na.rm=T) #null model zero line for plotting
y5 <- mean(f22$xcarbon5, na.rm=TRUE)  #
y5start <- mean(c5, na.rm=T) #null model zero line for plotting
TPA3 <- mean(f22$TPAchg, na.rm=TRUE) #
sapTPA3 <- mean(f22$sapTPAchg, na.rm=TRUE) #
TPA3start <-  mean(y6, na.rm=T)
sapTPA3start <- mean(sap6, na.rm=T)
x5 <- mean(f22$x5, na.rm=TRUE) #
x5start <-  mean(rx5, na.rm=T)
mortx5 <- mean(f22$mortx5, na.rm=TRUE) #
mortx5start <-   mean(mortx8, na.rm=T)

sapMORT <- mean(f22$sapMORTchg, na.rm=T)
sapRECR <- mean(f22$sapRECRchg, na.rm=T)

#calc y-axis error
R2 <- 100 #repetitions set to 10 for faster looping
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
enderrormort4 <- boot(f22$x5, meanFunc, R= R2)
enderrormortx5 <- sqrt(var(enderrormort4$t))*1.96
enderrorsapMORT <- std.error(f22$sapMORTchg, na.rm=T)
enderrorsapRECR <- std.error(f22$sapRECRchg, na.rm=T)
enderrorsp2error <- boot((f22$SISPchg-f22$sapavgTPAchg), meanFunc, R= R2)
enderrorsp22error <- sqrt(var(enderrorsp2error$t))*1.96

n <- length(f22$xcarbon) #conditional model n
nnull <- length(F33$avgx5) #null model n
print(c(n, ccond, c, nnull))

f23 <- c(x3, y3, y5, TPA3, sapTPA3, x5, mortx5, enderrory3, enderrory5, enderrorTPA, enderrorsapTPA, enderrorx5, enderrormortx5, n, v, col3, y3start, y5start, TPA3start, x5start, sapTPA3start, mortx5start, sapMORT, sapRECR, enderrorsapMORT, enderrorsapRECR, sp22, enderrorsp22error, y3pctdif)

}
}
F5 <- rbind(F5,F4)
}

F6 <- as.data.table(F5)
setnames(F6, c("x3","y3","y5", "TPA","sapTPA", "x5", "DIAmortality", "errory3","errory5", "errorTPA","errorsapTPA", "errorx5","DIAmorterror", "n", "v", "col3", "y3start", "y5start", "TPA3start", "x5start", "sapTPA3start", "DIAmortstart", "sapMORTpct", "sapRECRpct", "enderrorsapMORT", "enderrorsapRECR", "SISPratio", "SISPerror", "y3pctdif"))
F6 <- subset(F6, n>350) 

#convert radial growth from cm to mm
F6$x5 <- F6$x5*10 
F6$errorx5 <- F6$errorx5*10
#convert to percent difference
F6$sapMORTpct <-  (F6$sapMORTpct-100)
F6$sapRECRpct <-  (F6$sapRECRpct-100)
F6$sapCpct <-  ((F6$y3-F6$y5)/F6$y3)*100


#calc plot y-scales
sapscale  <-  30 #mean(F6$sapTPA3start)/mean(F6$TPA3start) 
Cymax  <-  max(F6$y3+F6$errory3)
Cymin  <-  min(F6$y3-F6$errory3)
TPAymax <- max(F6$sapTPA/sapscale+F6$errorsapTPA/sapscale, na.rm=T)
#TPAymax <- (max(F6$sapTPA+F6$errorsapTPA))/sapscale
TPAymin <- min(F6$TPA-F6$errorTPA)
#TPAymin <- (min(F6$sapTPA-F6$errorsapTPA))/sapscale
DIAmax <- max(F6$x5+F6$errorx5)
DIAmin <- min(F6$x5-F6$errorx5)
sapmax <- max(F6$sapMORTpct+F6$enderrorsapMORT)
sapmin <- min(F6$sapMORTpct-F6$enderrorsapMORT)

#plot only sapling data that is significantly different with scaling than larger trees. 
F6$sapTPA <- ifelse(((F6$sapTPA/sapscale)>(F6$TPA+F6$errorTPA))|((F6$sapTPA/sapscale)<(F6$TPA-F6$errorTPA)), F6$sapTPA, NA)

#set up 5-year precip plots using ggplot
F7 <- subset(F6, v==1)
F7$sapTPA <- ifelse(F7$col3==2,NA, F7$sapTPA)
v=1
y3start1 <-mean(F6$y3start)
TPA3start1 <-mean(F6$TPA3start) #mean(F2$TPAchg) #
sapTPA3start1 <-mean(F6$sapTPA3start)
x5start1 <- mean((F6$x5start)*10)
col4 <- c(dt,"grey35", pt)

#NPP plot
#limitssap <- aes(ymax=(y5/y5scale +errorsapTPA/sapscale)+TPA3start1, ymin=(sapTPA/sapscale-errorsapTPA/sapscale)+TPA3start1)
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
P1 <- ggplot(F7, aes(x=x3, y=(y3+y3start1), fill=as.factor(col3)))  +geom_hline(y=y3start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21)  + theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ theme(axis.text.x = element_blank())+ scale_y_continuous(name=expression(paste(Delta, " NPP (Mg C ha"^"-1","yr" ^"-1",")")), limits=c(Cymin, Cymax)+y3start1, breaks=c(1.1, 1.4, round(y3start1, 2), 2.0, 2.3), labels=c(1.1, 1.4, round(y3start1, 1), 2.0, 2.3)) + scale_x_continuous(element_blank(), breaks=c(-0.5, 0, 0.5))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c( "Dry 5-yr period", "Wet 5-yr period"), values=c(col4)) 
 
#Stem density plot
limits <- aes(ymax=(TPA+errorTPA)+TPA3start1, ymin=(TPA-errorTPA)+TPA3start1)
limitssap <- aes(ymax=(sapTPA/sapscale +errorsapTPA/sapscale)+TPA3start1, ymin=(sapTPA/sapscale-errorsapTPA/sapscale)+TPA3start1)
P2 <- ggplot(F7, aes(x=x3, y=(TPA+TPA3start1), fill=as.factor(col3)))  +geom_hline(y=TPA3start1, linetype=3) +geom_vline(linetype=3) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(-0.07, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ theme(axis.text.x = element_blank())+ scale_y_continuous(name=expression(paste(Delta, " Stem density (ha"^"-1","yr" ^"-1",")")), limits=c(TPAymin, TPAymax)+TPA3start1, breaks=c(2.0, 4.4866, 7.0), labels=c(expression("2.0"[" (-40)"]), expression("4.5" [" (35)"]), expression("7.0" [" (110)"]))) + scale_x_continuous(element_blank(), breaks=c(-0.5, 0, 0.5)) +scale_fill_manual(values=c(col4), breaks=NULL) +geom_errorbar(limitssap, width=0.01, colour="grey70") + geom_point(aes( y=(sapTPA/sapscale)+ TPA3start1, col=as.factor(col3)), size=1.5, pch=8) + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0.15,0), "cm")) +scale_colour_manual(values=c(col9), name=NULL, breaks=c("2"), label=c("(Saplings)"))+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21) 

#Radial growth plot
limits <- aes(ymax=(x5+errorx5)+x5start1, ymin=(x5-errorx5)+x5start1)
P3 <- ggplot(F7, aes(x=x3, y=(x5+x5start1), fill=as.factor(col3)))  +geom_hline(y=x5start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0.27,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black"))+ scale_y_continuous(name=expression(paste("Radial growth (mm yr" ^"-1",")")), limits=c(DIAmin, DIAmax)+x5start1, breaks=c(-0.05, 0, 0.05)+x5start1, labels=c(-0.05, 0, 0.05)+round(x5start1,2)) + scale_x_continuous(name=expression(paste("   Precipitation JJA 5-yr\nmean anomaly (mm day"^"-1",")")), breaks=c(-0.5, 0, 0.5)) +scale_fill_manual(values=c(col4)) + theme(axis.title.x = element_text(vjust=-0.85))

 
#set up 5-year TEMP plots using ggplot
F8 <- subset(F6, v==2) # & x3>0.2)
v=2
col8 <- c(ct,"grey35", wt)
#NPP plot
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
T1 <- ggplot(F8, aes(x=x3, y=(y3+y3start1), fill=as.factor(col3)))  +geom_hline(y=y3start1, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x =element_blank())+ scale_y_continuous(element_blank(), limits=c(Cymin, Cymax)+y3start1, breaks=c(1.1, 1.4, round(y3start1, 2), 2.0, 2.3)) + scale_x_continuous(element_blank(),breaks=c(-0.5, 0, 0.5, 1.0))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Cool 5-yr period", "Warm 5-yr period"), values=c(col8))

#stem density
limits <- aes(ymax=(TPA+errorTPA)+TPA3start1, ymin=(TPA-errorTPA)+TPA3start1)
T2 <- ggplot(F8,aes(x=x3, y=(TPA+TPA3start1), fill=as.factor(col3)))    +geom_hline(y=TPA3start1, linetype=3) +geom_vline(linetype=3) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm"))  + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_blank())+ theme(axis.text.x =(element_blank())) + scale_y_continuous(element_blank(), limits=c(TPAymin, TPAymax)+TPA3start1, breaks=c(2.0, 4.4866, 7.0)) + scale_x_continuous(element_blank())  +scale_fill_manual(values=c(col8), breaks=NULL) +geom_errorbar(limitssap, width=0.01, colour="grey70") + geom_point(aes( y=(sapTPA/sapscale)+ TPA3start1, col=as.factor(col3)), size=1.5, pch=8) + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_colour_manual(values=c(col8), name=NULL, breaks=c("2"), label=c("(Saplings)"))+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21)  

#Radial growth plot
limits <- aes(ymax=(x5+errorx5)+x5start1, ymin=(x5-errorx5)+x5start1)
T3 <- ggplot(F8, aes(x=x3, y=(x5+x5start1), fill=as.factor(col3)))   +geom_hline(y=x5start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y =element_blank())+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(element_blank(), limits=c(DIAmin, DIAmax)+x5start1, breaks=c(-0.05, 0, 0.05)+x5start1) + scale_x_continuous(name=("Temperature JJA 5-yr\nmean anomaly (C)")) +scale_fill_manual(values=c(col8))


#set up 1-year soil moist. plots using ggplot
F9 <- subset(F6, v==3)
v=3
col9 <- c(de,"grey35", pe)
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
#NPP plot
R1 <- ggplot(F9,aes(x=x3, y=(y3+y3start1), fill=as.factor(col3))) +geom_hline(y=y3start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.02, colour="grey50")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x =(element_blank())) + scale_y_continuous(element_blank(), limits=c(Cymin, Cymax)+y3start1, breaks=c(1.1, 1.4, round(y3start1, 2), 2.0, 2.3)) + scale_x_continuous(element_blank()) + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c( "Dry summer", "Wet summer"), values=c(col9))

#Stem density plot
limits <- aes(ymax=(TPA+errorTPA)+TPA3start1, ymin=(TPA-errorTPA)+TPA3start1)
limitssap <- aes(ymax=(sapTPA/sapscale +errorsapTPA/sapscale)+TPA3start1, ymin=(sapTPA/sapscale-errorsapTPA/sapscale)+TPA3start1)
R2 <- ggplot(F9,aes(x=x3, y=(TPA+TPA3start1), fill=as.factor(col3)))   +geom_hline(y=TPA3start1, linetype=3) +geom_vline(linetype=3)+theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x =(element_blank()))+ scale_y_continuous(element_blank(), limits=c(TPAymin, TPAymax)+TPA3start1, breaks=c(2.0, 4.4866, 7.0)) + scale_x_continuous(element_blank()) +scale_fill_manual(values=c(col9), breaks=NULL) +geom_errorbar(limitssap, width=0.01, colour="grey70") + geom_point(aes( y=(sapTPA/sapscale)+ TPA3start1, col=as.factor(col3)), size=1.5, pch=8) + theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_colour_manual(values=c(col9), name=NULL, breaks=c("2"), label=c("(Saplings)"))+ geom_errorbar(limits, width=0.02, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21)

#radial growth plot
limits <- aes(ymax=(x5+errorx5)+x5start1, ymin=(x5-errorx5)+x5start1)
R3 <- ggplot(F9, aes(x=x3, y=(x5+x5start1), fill=as.factor(col3)))   +geom_hline(y=x5start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.02, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21)  +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_blank())+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(element_blank(), limits=c(DIAmin, DIAmax)+x5start1, breaks=c(-0.05, 0, 0.05)+x5start1) + scale_x_continuous(name=("Soil moisture JJA 1-yr\nmin/max (z-score)")) +scale_fill_manual(values=c(col9))

#set up plot window and plot all together
dev.new(width=7, height=8)
grid.arrange(P1,T1,R1,P2,T2,R2,P3,T3,R3, ncol=3)
#grid.text(0.045, unit(0.575,"npc") - unit(1,"line"), rot=90, label="(30)", gp=gpar(cex=0.75))

#######END Fig 3.  Climate Correlations######################################################


############# Figure S1 -  SAPLING mortality and recruitment



F7s <- subset(F6, v==1)
v=1
col4s <- c(dt,"grey35", pt)
limitssap <- aes(ymax=(sapMORTpct+ enderrorsapMORT), ymin=(sapMORTpct-enderrorsapMORT))
P1 <- ggplot(F7s, aes(x=x3, y=(sapMORTpct), colour=as.factor(col3)))  +geom_hline(y=0, linetype=3) +geom_vline(linetype=3) + geom_errorbar(limitssap, width=0.01, colour="grey50") + geom_point( size=1.5, pch=8) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) +theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) +  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ scale_y_continuous(name=expression(paste(Delta, " Sapling mortality (%)")), limits=c(sapmin, sapmax), breaks=c(-50,0,50, 100)) + scale_x_continuous(name=("Precipitation JJA 5-yr\nmean anomaly (mm/day)"), breaks=c(-0.5, 0, 0.5))  +scale_colour_manual( name=NULL, breaks=c("1", "3"), label=c("Dry 5-yr period", "Wet 5-yr period"), values=c(col4s))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) 



F8s <- subset(F6, v==2)
v=1
col5s <- c(ct,"grey35", wt)
limitssap <- aes(ymax=(sapMORTpct+ enderrorsapMORT), ymin=(sapMORTpct-enderrorsapMORT))
P2 <- ggplot(F8s, aes(x=x3, y=(sapMORTpct), colour=as.factor(col3)))  +geom_hline(y=0, linetype=3) +geom_vline(linetype=3) + geom_errorbar(limitssap, width=0.01, colour="grey50") + geom_point( size=1.5, pch=8) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) +theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) +  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ scale_y_continuous(name=expression(paste( Delta, " Sapling mortality (%)")), limits=c(sapmin, sapmax), breaks=c(-50,0, 50, 100)) + scale_x_continuous(name=("Temperature JJA 5-yr\nmean anomaly (C)"), breaks=c(-0.5, 0, 0.5, 1))  +scale_colour_manual( name=NULL, breaks=c("1", "3"), label=c("Cool 5-yr period", "Warm 5-yr period"), values=c(col5s))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) 


F9s <- subset(F6, v==3)
col6s <- c(de,"grey35", pe)
limitssap <- aes(ymax=(sapMORTpct+ enderrorsapMORT), ymin=(sapMORTpct-enderrorsapMORT))
P3 <- ggplot(F9s, aes(x=x3, y=(sapMORTpct), colour=as.factor(col3)))  +geom_hline(y=0, linetype=3) +geom_vline(linetype=3) + geom_errorbar(limitssap, width=0.01, colour="grey50") + geom_point( size=1.5, pch=8) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) +theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) +  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ scale_y_continuous(name=expression(paste(Delta, " Sapling mortality (%)")), limits=c(sapmin, sapmax), breaks=c(-50,0,50, 100)) + scale_x_continuous(name=("Soil moisture JJA 1-yr\nmin/max (z-score)"), breaks=c(-1, 0, 1))  +scale_colour_manual( name=NULL, breaks=c("1", "3"), label=c("Dry summer", "Wet summer"), values=c(col6s))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) 

dev.new(width=7, height=3)
grid.arrange(P1,P2,P3, ncol=3)




####Fig. S2. undisturbed Carbon climate plots ########################
#NOTE: run above analysis excuding disturbance


F6 <- as.data.table(F5)
setnames(F6, c("x3","y3","y5", "TPA","sapTPA", "x5", "DIAmortality", "errory3","errory5", "errorTPA","errorsapTPA", "errorx5","DIAmorterror", "n", "v", "col3", "y3start", "y5start", "TPA3start", "x5start", "sapTPA3start", "DIAmortstart", "sapMORTpct", "sapRECRpct", "enderrorsapMORT", "enderrorsapRECR", "SISPratio", "SISPerror", "y3pctdif"))
F6 <- subset(F6, n>5) # n>175) 

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
P1 <- ggplot(F7, aes(x=x3, y=(y5+y5start1), fill=as.factor(col3)))  +geom_hline(y=y5start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey50") + geom_point(colour="grey25", size=2.5, pch=21)  + theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(name=expression(paste(Delta, " NPP (Mg C ha"^"-1","yr" ^"-1",")")), limits=c(Cymin, Cymax)+y5start1, breaks=c(1.1, round(y5start1,2), 1.8 )) + scale_x_continuous(name=("Precipitation JJA 5-yr\nmean anomaly (mm/day)"), breaks=c(-0.5, 0, 0.5))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Wet 5-yr period", "Dry 5-yr period"), values=c(col4)) 

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



#############################################

####Fig. S3 & S4. deciduous and coniferous Carbon climate plots ########################
#NOTE: run above analysis with SISP>=300 (deciduous) & SISP<300 (coniferous)


F6 <- as.data.table(F5)
setnames(F6, c("x3","y3","y5", "TPA","sapTPA", "x5", "DIAmortality", "errory3","errory5", "errorTPA","errorsapTPA", "errorx5","DIAmorterror", "n", "v", "col3", "y3start", "y5start", "TPA3start", "x5start", "sapTPA3start", "DIAmortstart", "sapMORTpct", "sapRECRpct", "enderrorsapMORT", "enderrorsapRECR", "SISPratio", "SISPerror"))
F6 <- subset(F6, n>1) 

#calc plot y-scales
Cymax  <-  max(F6$y3+F6$errory3)
Cymin  <-  min(F6$y3-F6$errory3)

F7 <- subset(F6, v==1)
F7$sapTPA <- ifelse(F7$col3==2,NA, F7$sapTPA)
v=1
y3start1 <-(F6$y3start[1])
TPA3start1 <-(F6$TPA3start[1]) #mean(F2$TPAchg) #
sapTPA3start1 <-(F6$sapTPA3start[1])
x5start1 <- ((F6$x5start[1])*10)
col4 <- c(dt,"grey35", pt)
#NPP plot
limitssap <- aes(ymax=(y3/y3scale +errorsapTPA/sapscale)+TPA3start1, ymin=(sapTPA/sapscale-errorsapTPA/sapscale)+TPA3start1)
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
P1 <- ggplot(F7, aes(x=x3, y=(y3+y3start1), fill=as.factor(col3)))  +geom_hline(y=y3start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey30") + geom_point(colour="grey25", size=2.5, pch=21)  + theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA))  + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size=11, colour="black"))+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(name=expression(paste(Delta, " NPP (Mg C ha"^"-1","yr" ^"-1",")")), limits=c(Cymin, Cymax)+y3start1,breaks=c(1.0, 1.4, 1.8, 2.2)) + scale_x_continuous(name=("Precipitation JJA 5-yr\nmean anomaly (mm/day)"), breaks=c(-0.5, 0, 0.5))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Wet 5-yr period", "Dry 5-yr period"), values=c(col4)) 

#set up 5-year TEMP plots using ggplot
F8 <- subset(F6, v==2) # & x3>0.2)
v=2
col8 <- c(ct,"grey35", wt)
#NPP plot
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
T1 <- ggplot(F8, aes(x=x3, y=(y3+y3start1), fill=as.factor(col3)))  +geom_hline(y=y3start1, linetype=3)+geom_vline(linetype=3)+ geom_errorbar(limits, width=0.01, colour="grey30")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(element_blank(), limits=c(Cymin, Cymax)+y3start1, breaks=c(1.0, 1.4, 1.8, 2.2))   + scale_x_continuous(name=("Temperature JJA 5-yr\nmean anomaly (C)"),breaks=c(-0.5, 0, 0.5, 1.0))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Cool 5-yr period", "Warm 5-yr period"), values=c(col8))

#set up 1-year soil moist. plots using ggplot
F9 <- subset(F6, v==3)
v=3
col9 <- c(de,"grey35", pe)
limits <- aes(ymax=(y3+errory3)+y3start1, ymin=(y3-errory3+y3start1))
#NPP plot
R1 <- ggplot(F9,aes(x=x3, y=(y3+y3start1), fill=as.factor(col3))) +geom_hline(y=y3start1, linetype=3) +geom_vline(linetype=3)+ geom_errorbar(limits, width=0.02, colour="grey30")+ geom_point(colour="grey25", size=2.5, pch=21) +theme(axis.line = element_line( size = 0.35)) + theme(axis.ticks.margin = unit(0.06, "cm")) + theme(axis.ticks.length = unit(0.1, "cm")) + theme(legend.position = "none",plot.margin=unit(c(0,0,0,0), "cm")) +theme( panel.background = element_rect(fill=NA),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.background = element_rect(fill=NA)) + theme(axis.text.y = element_blank())+ theme(axis.text.x = element_text( hjust = 0.5, size=11, colour="black")) + scale_y_continuous(element_blank(), limits=c(Cymin, Cymax)+y3start1, breaks=c(1.0, 1.4, 1.8, 2.2))  + scale_x_continuous(name=("Soil moisture JJA 1-yr\nmin/max (z-score)"))+ theme(legend.title=element_blank(), legend.key=element_rect(fill="white"), legend.key.size = unit(0.35, "cm"), legend.position = c(0.5,0.9),plot.margin=unit(c(0,0,0,0), "cm")) +scale_fill_manual( name=NULL, breaks=c("1", "3"), label=c("Wet summer", "Dry summer"), values=c(col9))

dev.new(width=7, height=3)
grid.arrange(P1,T1,R1, ncol=3)



#############################################

#CLIMATE DATA
setnames(F8, c("TEMPchg","carbon1","carbon5", "stems","saplingstems", "radialgrow", "DIAmortality", "carbon1.2SE","carbon5.2SE", "stems.2SE","saplingstems.2SE", "radialgrow.2SE", "DIAmort.2SE", "n", "v", "col3", "carbon1.start", "carbon5.start", "stem.start", "radial.start", "saplingstem.start", "DIAmort.start"))

#warming results
print((F8[22])) #0.5 degree warming
print((F8[29])) #1.0 degree warming
print(c( "quantile", qt2[22])) #0.5 degree warming
print(c( "quantile", qt2[22])) #0.5 degree warming

# % difference from null
(1-((y3start1-F8$carbon1[22])/y3start1))*100 #percent 0.5 degree climate warming difference
(1-((y3start1-F8$carbon1[29])/y3start1))*100 #percent 1.0 degree climate warming difference

# % NPP change from saplings
((F8$carbon1[22]-F8$carbon5[22])/F8$carbon1[22])*100 #0.5 degree warming
((F8$carbon1[29]-F8$carbon5[29])/F8$carbon1[29])*100 #1.0 degree warming

##############END FIG. 3/S1#########################################################################################



####################################################################################################################
############## Extended data TABLE 3 - CLimate stats ############################################################################
#########
###

#ratio of mean conditions with Null weighted spatially by Condition
meanFunc2 <- function(x,i){(mean(x$variCOND[i], na.rm=TRUE)/mean((x$variNULL[i]), na.rm=TRUE))} 
#meanFunc3 <- function(x,i){(mean(x$deadcarbon5sumYR[i], na.rm=TRUE)/mean((x$mortality2[i]), na.rm=TRUE))} 

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

F2$Drysummer <- ifelse(F2$SWdrought5chg <= quantile(F2$SWdrought5chg, 0.17, na.rm=T), 1,0)
F2$Wetsummer <- ifelse(F2$SWpluvial5chg >= quantile(F2$SWpluvial5chg, 0.83, na.rm=T), 1,0)
F2$Dryperiod <- ifelse(F2$PREC5growchg <= quantile(F2$PREC5growchg, 0.17, na.rm=T), 1,0)
F2$Wetperiod <- ifelse(F2$PREC5growchg >= quantile(F2$PREC5growchg, 0.83, na.rm=T), 1,0)
F2$SPCDdiversity2 <- F2$SPCDdiversity-F2$PREV_SPCDdiversity

#set disturbance conditions to binary
F2$fire <- ifelse(F2$fire>0, 1,0)
F2$insect <- ifelse(F2$insect>0, 1,0)
F2$disease <- ifelse(F2$disease>0, 1,0)
F2$animal <- ifelse(F2$animal>0, 1,0)
F2$weather <- ifelse(F2$weather>0, 1,0)
F2$SPCDdiversity2 <- ifelse(F2$SPCDdiversity2>0, 1,0)

vari <- c("SWdrought5chg",  "SWpluvial5chg", "PREC5growchg", "TEMP5growchg") 
vari1 <- subset(F2, select=c(SWdrought5chg, SWpluvial5chg, PREC5growchg, TEMP5growchg)) 
vari2 <- c("fire", "insect", "disease", "animal", "weather", "Drysummer", "Wetsummer", "Dryperiod", "Wetperiod", "SPCDdiversity2")

###ADJUST TOP/bottom QUANTILE HERE
quantile1<-0.17 #0.035

F5b <- NULL
for(v2 in 1:10)
{
for(v in 1:4)
{
for(c in 1:2)
{
if(c==1 & v>2){f22 <- subset(F2, get(vari[v])<= quantile(vari1[[v]], quantile1, na.rm=TRUE)) }	#0.17 & 0.03 set slimate extreme bottom percentile
if(c==2 & v>2){f22 <- subset(F2, get(vari[v])>= quantile(vari1[[v]], 1-quantile1, na.rm=TRUE)) }	#0.83 & 0.965 set slimate extreme top percentile
if(c==1 & vari[v]=="SWdrought5chg"){f22 <- subset(F2, get(vari[v])<= quantile(vari1[[v]], quantile1, na.rm=TRUE) & SWpluvial5chg<= quantile(F2$SWpluvial5chg, 0.83, na.rm=TRUE)) }	#
if(c==2 &  vari[v]=="SWpluvial5chg"){f22 <- subset(F2, get(vari[v])>= quantile(vari1[[v]], 1-quantile1, na.rm=TRUE) & SWdrought5chg>= quantile(F2$SWdrought5chg, 0.17, na.rm=TRUE) ) }	
if(c==1 & vari[v]=="SWpluvial5chg"){next}
if(c==2 & vari[v]=="SWdrought5chg"){next}


#print(mean(f22$TEMP5growchg))
f23 <- subset(F2, get(vari[v])< quantile(vari1[[v]], 0.7, na.rm=TRUE) & get(vari[v])> quantile(vari1[[v]], 0.3, na.rm=TRUE))
F33 <- subset(f23, LATLON %in% (f22$LATLON)) # not neccesaary bc weighting will eliminate latlons not in f22!


f22$Csamp <- (ave(f22$DIAmean>0,f22$LATLON, FUN=function(x) sum(x, na.rm=TRUE)))/f22$nsampTOT #weighting factor of number of Cond data/total data in each grid cell
F33$Nsamp <- (ave(F33$DIAmean>0,F33$LATLON, FUN=function(x) sum(x, na.rm=TRUE)))/F33$nsampTOT #weighting factor of number of NULL data/total data in each grid cell

f24 <- subset(f22, !duplicated(LATLON), select=c(LATLON, Csamp))
F33 <- merge(F33, f24, all.x=TRUE, by="LATLON")
F33$CtoNratio <- F33$Csamp/F33$Nsamp

variCOND <- subset(f22, select=c("fire", "insect", "disease", "animal", "weather", "Drysummer", "Wetsummer", "Dryperiod", "Wetperiod", "SPCDdiversity2"))
variNULL <- subset(F33, select=c("fire", "insect", "disease", "animal", "weather", "Drysummer", "Wetsummer", "Dryperiod", "Wetperiod", "SPCDdiversity2", "CtoNratio"))

 #make dataframe for error bootstrapping
varNullR <- variNULL[[v2]]*variNULL$CtoNratio
ffboot <- as.data.frame(varNullR)
l1 <- length(variCOND[[v2]])
l2 <- length(variNULL[[v2]])
if(l1>l2){
	varNullR <- rep(NA, l1-l2)
	ffboot2 <- as.data.frame(varNullR)
	ffboot <- rbind(ffboot, ffboot2)
	ffboot$variCOND <- c(variCOND[[v2]])

	}
if(l1<=l2){
	ffboot$variCOND <- c(variCOND[[v2]], rep(NA, l2-l1))
	}
	
setnames(ffboot,c("varNullR","variCOND"), c("variNULL","variCOND"))

RatioBoot <- boot(ffboot, meanFunc2, R=1000)					
error <- round(sqrt(var(RatioBoot$t, na.rm=T))*100, 1)

#subset data based on condition to get actual ratio in Null climate and Conditional climates rather than Bootstrapped mean is approx. the same
NullVarRatiod <- (sum(variNULL[[v2]]*variNULL$CtoNratio)/length(variNULL[[v2]]))*100
CondVarRatiod <- (sum(variCOND[[v2]])/length(variCOND[[v2]]))*100
rateRatiod <- ((CondVarRatiod/NullVarRatiod)*100)-100
sig <- ifelse(abs(rateRatiod)>abs(error*1.96), "***", "")

Clim <- ifelse(vari[v]=="PREC5growchg" & c==2, "5-yr Wet period", NA)
Clim <- ifelse(vari[v]=="PREC5growchg" & c==1, "5-yr Dry period", Clim)
Clim <- ifelse(vari[v]=="TEMP5growchg" & c==1, "Cool period", Clim)
Clim <- ifelse(vari[v]=="TEMP5growchg" & c==2, "Warm period", Clim)
Clim <- ifelse(vari[v]=="SWpluvial5chg" & c==2, "Wet summer", Clim)
Clim <- ifelse(vari[v]=="SWdrought5chg" & c==1, "Dry summer", Clim)

n <- sum(variCOND[[v2]]==1)

f34 <- c(vari2[v2],Clim, n, round(CondVarRatiod,1), round(NullVarRatiod,1), round(rateRatiod,1), error, sig) #, round(Vgrow2,2), round(Vgrow,2), error2)
F5b <- rbind(F5b, f34)
print(c(v2, v, c))
}
}
}

F6b <- as.data.table(F5b)
setnames(F6b, c("Disturbance type", "Climate condition", "n", "%CondVarRatiod", "%NullVarRatiod", "%rateRatiod", "BootRate Error", "p < 0.05") )#, "Vgrow2","Vdead", "error2"))
print(na.omit(F6b))
#write.csv(na.omit(F6b))
#################################################################################################
##END MAIN FIGURES and climate conditions supplimentary figures








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
dev.new(width=8, height=5)
plot(c((F2$PREV_DIAmean),(F2$DIAmean)),c((F2$PREV_TPAsum),(F2$TPAsum)), type ="n", col="grey94", 
	main=c("Empirical Model of 5-year forest carbon change, DBH>=12.7 cm", "FIA survey 1997-2012, Eastern USA") 
	,ylim=c(60,155), xlim=c(23, 40), mgp=c(2,1,0)
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
arrows((F2$PREV_DIAmean),(F2$PREV_TPAsum),(F2$DIAmean),(F2$TPAsum), col= "grey85", length=0.065, angle=22,lwd=0.9)   #col= "grey85"

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


for(l in 1:4){
if(l==1){	  ##Reforestation
STOCKval<-c(5,15)
QMDvalues<-c(seq(21.7,30.7, 3), 35) #c(18,21,24,29, 45)  #
s1<-5
s2<-5
s3<-5
s4<-5
increment=1}

if(l==2){	  ##Early succession
STOCKval<-c(seq(0, 121, 121/6)) 
QMDvalues<-c(seq(12.7,21.7,1.8)) #c(seq(12,18.4,3.2))  #
s1<-30
s2<-10
s3<-10
s4<-5
increment=1}

if(l==3){	  ##middle
STOCKval<-c(seq(15, 85, 65/5),98,121) #121,121,121
QMDvalues<-c(seq(21.7,t3, 10.3/6))  #
s1<-30
s2<-10
s3<-20
s4<-2
increment=1}

if(l==4){	  ##moving average end
STOCKval<-c(5, 15, NA, 30, 38, NA, NA,20, 40, 55, 55, 60, 75, NA, NA,NA, 65, 85, 100, 85, 98, 121) #121, 121) #121,121,121
QMDvalues<-c(seq(t3-1,45,1.8))  #
s1<-100
s2<-50
s3<-40
s4<-5
increment=3}

S2<-NULL
g<-length(QMDvalues-1)
r<-length(STOCKval-1)
 for(s in 1:r)
 { S11r<-subset(F2, PREV_STOCKINGmid>= STOCKval[s] & PREV_STOCKINGmid <STOCKval[s+increment]) 	
	for(h in 1:g)
	{ 
	S11q<-subset(S11r, PREV_DIAmean>=(QMDvalues[h]) & PREV_DIAmean<(QMDvalues[h+increment]))
print(length(S11q$PREV_DIAmean))
	if(length(S11q$PREV_DIAmean)>s1 & s==1 | length(S11q$PREV_DIAmean)>s2 & s==2 | length(S11q$PREV_DIAmean)>s3  & s==3 | length(S11q$PREV_DIAmean)>s4  & s>=4)
		{				
	#arrows((mean(S11q$PREV_DIAmean, na.rm=TRUE)),(mean(S11q$PREV_TPAsum, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE)),col=("grey15"), length=0.085, angle=25,lwd=3+0.75) 
	arrows((mean(S11q$PREV_DIAmean, na.rm=TRUE)),(mean(S11q$PREV_TPAsum, na.rm=TRUE)),(mean(S11q$DIAmean, na.rm=TRUE)),(mean(S11q$TPAsum, na.rm=TRUE)), col= "grey75", length=0.08, angle=28,lwd=3)  # VectCols[s] col=s*1.8+10,  s*1.2+7 # s*2+10 #stockcol
	}
	}
	}
	}

############  END CARBON PLOT SET UP
#	
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
F2<-subset(F2, select=c(PLT_CN, PREV_STOCKINGmid, PREV_DIAmean, PREV_Csum, DIAmean, Csum, STOCKINGmid,  REMPER, DIAchg, Cchg, STOCKchg, TEMP5growchg, PREC5growchg, SWdrought5chg, SWpluvial5chg, cutting, STDORGCD, cut, cutDIAmean, TPAsum, NOdisturb, alldisturb, allmort, EHstocksum, NMHstocksum, SMHstocksum, LHstocksum)) #, PREV_SPCDdiversity, SPCDdiversity, PREV_SPGRPCD, SPGRPCD))# 
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

time<-150 #years forward modeling

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

col6<-c("ct",  "de",  "wt",  "dt", "pe", "pt",  "av", "ud", "ds", "lg", "lg2") 
labels<-c("5-yr cool period",  "dry summer",  "5-yr warm period",  "5-yr dry period", "wet summer", "5-yr wet period",  "average", "undisturbed", "disturbed", "average w/logging", "all mortality" ) 
col5<-c(ct,de,wt,dt, pe, pt, av, ud, ds,lg, lg2)




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
F5b<-subset(F5, samp<  (sum(F5$PLT_CN>0)*0.99)) #7200
F5c<-subset(F5, samp< (sum(F5$PLT_CN>0)*0.99))
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
if(c==7){F2<-subset(F5c)}
if(c==8){F2<-subset(F5b, NOdisturb==1)}
if(c==9){F2<-subset(F5b, alldisturb>0)}
if(c==10){F2<-subset(F5harvestb)}
if(c==11){F2<-subset(F5b, allmort>0)}	
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
F3d<-subset(F3b, select=c(CsumSTART, DIAmeanSTART, STOCKINGmidSTART, STOCKINGmid, DIAmean, Csum, PLT_CN , PREV_DIAmean, PREV_Csum, Cchg, DIAchg, REMPER, TPAsum, EHstocksum, NMHstocksum, SMHstocksum, LHstocksum))


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

F3c<-subset(F3c, select=c(CsumSTART, DIAmeanSTART, STOCKINGmidSTART, STOCKINGmid, DIAmean, Csum, PLT_CN, PREV_DIAmean, PREV_Csum, Cchg, DIAchg, REMPER, TPAsum, EHstocksum, NMHstocksum, SMHstocksum, LHstocksum))
F3b<-subset(F3b, select=c(CsumSTART, DIAmeanSTART, STOCKINGmidSTART, STOCKINGmid, DIAmean, Csum, PLT_CN, PREV_DIAmean, PREV_Csum, Cchg, DIAchg, REMPER, TPAsum, EHstocksum, NMHstocksum, SMHstocksum, LHstocksum))
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
EH<- mean(F3$EHstocksum, na.rm=T)
NMH<- mean(F3$NMHstocksum, na.rm=T)
SMH<- mean(F3$SMHstocksum, na.rm=T)
LH<- mean(F3$LHstocksum, na.rm=T)


arrows(x1,y1, x2, y2 ,col=("black"), length=0.07, angle=25,lwd=2) 
arrows(x1,y1, x2, y2 ,col=col5[c], length=0.07, angle=25,lwd=1.5) 
if(i==10 | i==20 | i==30 |i==40 )
{arrows(x1,y1, x2, y2 ,col=("salmon"), length=0.07, angle=25,lwd=2)	} #denote 100 years of growth with red vectors

#gather arrow data
Fbindall<-c(x1,y1, x2, y2, TPA2, c, Caccum, i, EH, NMH, SMH, LH)
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

Fbind<-c(x2, y2, TPA2, endCerror, endQMDerror, col6[c], Caccum, meancond)
Fbind2<-rbind(Fbind2, Fbind)
print(Fbind2)
}
}
}
#box()


#plot vectors on carbon mapping
Fbindall2<-as.data.table(Fbindall2)
setnames(Fbindall2, c("x1", "y1","x2", "y2", "TPA2", "c", "Caccum", "i", "EH", "NMH", "SMH", "LH"))


plot(Fbindall2$i*5, Fbindall2$NMH, type = "l", ylim=c(5,22))

lines((Fbindall2$i*5), Fbindall2$EH, col="green")
lines((Fbindall2$i*5), Fbindall2$LH, col="blue")
lines((Fbindall2$i*5), Fbindall2$SMH, col="red")



c1<-(1:11)
foreach(r=1:11) %do%{
fb2<-subset(Fbindall2, c==c1[r])
arrows(fb2$x1, fb2$y1, fb2$x2, fb2$y2,col="black", length=0.07, angle=25,lwd=2 )
arrows(fb2$x1, fb2$y1, fb2$x2, fb2$y2,col=col5[fb2$c], length=0.07, angle=25,lwd=1.5 )
Fbindall3<-subset(fb2, i==10 | i==20| i==30 )
arrows(Fbindall3$x1, Fbindall3$y1, Fbindall3$x2, Fbindall3$y2,col="salmon", length=0.07, angle=25,lwd=1.5 )
}

Fbind3<-as.data.table(Fbind2)
setnames(Fbind3,c("x2", "y2", "TPA2", "endCerror", "endQMDerror", "cond", "Caccum", "meancond"))

shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
	theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
	
	xy <- xy.coords(x,y)
	xo <- r*strwidth('A')
	yo <- r*strheight('A')
	for (i in theta) {
		text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
		  labels, col=bg, ... )
	}
	text(xy$x, xy$y, labels, col=col, ... )
}

labels3<-c("all forests w/ logging", "disturbed", "all mortality", "recurring wet summers", expression(paste("+0.37 mm day"^"-1","precip." )), expression(paste("+0.33"^"o","C warming" )),"all forests", expression(paste("-0.24 mm day"^"-1","precip." )), "recurring dry summers", expression(paste("-0.35"^"o","C cooling" )), "undisturbed")


y2<-as.numeric(as.character(Fbind3$y2))
x2<-as.numeric(as.character(Fbind3$x2))
y2<-c(72, 83, 88, seq(95,135, (40/6)), 148)
x2=(y2/14)+26.6
col55<-c(lg, ds, lg2, pe, pt, wt, av, dt, de, ct, ud)
w<-"white"
b<-"grey15"
colwb<-c(w, w, w, w, w, b, b, b, b, b, b)
xoff<-c(-5, -2.5, -2.5, 0, 0, 0, 0, 0, 0, 0, 1)
cond<-as.character(Fbind3$cond)
shadowtext(x2+xoff,y2, labels3, cex=0.8, pos=4, col= colwb, bg=col55) #, adj= -0.1
box()
##############END CARBON FORWARD PROJECTION
##
#
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
F2<-subset(F2, select=c(PLT_CN, PREV_STOCKINGmid, PREV_DIAmean, PREV_Csum, DIAmean, Csum, STOCKINGmid,  REMPER, DIAchg, Cchg, STOCKchg, TEMP5growchg, PREC5growchg, SWdrought5chg, SWpluvial5chg, cutting, STDORGCD, cut, cutDIAmean, TPAsum, NOdisturb, alldisturb, allmort)) #, PREV_SPCDdiversity, SPCDdiversity, PREV_SPGRPCD, SPGRPCD))# 
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

col6<-c("ct",  "de",  "wt",  "dt", "pe", "pt",  "av", "ud", "ds", "lg", "lg2") 
labels<-c("5-yr cool period",  "dry summer",  "5-yr warm period",  "5-yr dry period", "wet summer", "5-yr wet period",  "average", "undisturbed", "disturbed", "average w/logging", "all mortality" ) 
col5<-c(ct,de,wt,dt, pe, pt, av, ud, ds,lg, lg2)


qq2<-c(seq(0.3,0.6,0.2),rep(0.6, 500))
qq3<-c(seq(0.9,0.4,-0.2), rep(0.4, 500))
Fbindeq<-NULL
Fbindeq2<-NULL
Fbindeq3<-NULL

print("Run countdown:")
###
##for loop start

#ERROReqbind<-foreach(r=1:10, .combine=rbind) %dopar%{  #iterate model scenario 10 times and rbind results
for(r in 1:10) {  #iterate model scenario 10 times and rbind results #CAUTION THIS TAKES A LOONNG TIME  ~20 minutes - consider foreach parrallel proccessing

for(c in 1:11){   #run each scenario

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
if(c==9){F2<-subset(F5, alldisturb>0)}
if(c==10){F2<-subset(F5harvest)}
if(c==11){F2<-subset(F5, allmort>0)}	

F2$samp<-sample(1:sum(F2$PLT_CN>0, na.rm=T),sum(F2$PLT_CN>0, na.rm=T),replace=FALSE)
F2<-subset(F2, samp<=(7582*0.95)) 

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
cond<-as.data.table(c("ct",  "de",  "wt",  "dt", "pe", "pt",  "av", "ud", "ds", "lg", "lg2" ) )
Feq3<-cbind(Feq2, cond)
setnames(Feq3,c("ySD", "cond"))

#######END PREJECTION ERROR ANALYSIS######################################
#
#
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
 


F2<-subset(RM5mergedfullcut1991.11.17.2014, PREVcarbon1sum>0 & PREVSTOCKING5mid>=60 ) 
#Convert variables to standard names and metric units
F2$DIAmean<-(F2$DIA5meanalive* 2.54)			#mean tree diameter by plot at resurvey (cm)
F2$PREV_DIAmean<-(F2$PREVDIA5meanalive *2.54)	#mean tree diameter by plot at initial survey (cm)
F2$TPAsum<-(F2$endsumTPA /0.404686)			#number of trees per hectare by plot at resurvey
F2$PREV_TPAsum<-(F2$startsumTPA /0.404686)	#number of trees per hectare by plot at initial survey (cm)
F2$PREV_STOCKINGmid<-F2$PREVSTOCKING5mid	#calculated relative stocking by plot at initial survey
 
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

F2$average<-1
F2$nsampTOT<-ave(F2$DIAmean>0,F2$LATLON, FUN=function(x) sum(x, na.rm=TRUE))
F2$NOdisturb<-ifelse((F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)==0,1,NA) 
F2$alldisturb<-ifelse((F2$insect+F2$fire+F2$weather+F2$disease+F2$animal)>0 & (F2$unknowndamage+F2$vegetation)==0, 1, NA)
F2$allmort<-ifelse((F2$insect+F2$fire+F2$weather+F2$disease+F2$animal+F2$unknowndamage+F2$vegetation)>0, 1, NA)

vari<-c( "TEMP5growchg", "SWdrought5chg",  "TEMP5growchg",  "PREC5growchg", "SWpluvial5chg", "PREC5growchg", "average", "NOdisturb", "alldisturb", "allmort", "cutting") 

cond<-c("ct",  "de",  "wt",  "dt", "pe", "pt",  "av", "ud", "ds", "lg2", "lg" ) 
vari1<-subset(F2, select=c(noquote(vari))) 
#vari2<-c("fire")

F5<-subset(F2, cut==0 & STDORGCD==0 )
F5harvest<-subset(F2, STDORGCD==0 ) #, cutDIAmean>8.7 | is.na(cutDIAmean)==T )

qq3[i]<-0.4
qq2[i]<-0.6
F5b<-NULL
#for(v2 in 1:1)
#{
for(c in 1:11)
{
#same adjustment as above for largest effect with stable eq points
if(c==2){F2<-subset(F5, SWdrought5chg<=quantile(F5$SWdrought5chg, qq3[i]+0.1, na.rm=T) & SWpluvial5chg<=quantile(F5$SWpluvial5chg, 0.85, na.rm=T))}
if(c==1){F2<-subset(F5, TEMP5growchg<=quantile(F5$TEMP5growchg, qq3[i]-0.15, na.rm=T))}
if(c==3){F2<-subset(F5, TEMP5growchg>=quantile(F5$TEMP5growchg, qq2[i]+0.1,  na.rm=T))} 
if(c==4){F2<-subset(F5, PREC5growchg<=quantile(F5$PREC5growchg, qq3[i]+0.05, na.rm=T))}
if(c==5){F2<-subset(F5, SWpluvial5chg>=quantile(F5$SWpluvial5chg, qq2[i]+0.15, na.rm=T) & SWdrought5chg>=quantile(F5$SWdrought5chg, 0.15, na.rm=T))}
if(c==6){F2<-subset(F5, PREC5growchg>=quantile(F5$PREC5growchg, qq2[i]+0.15, na.rm=T)) }
if(c==7){F2<-subset(F5)}
if(c==8){F2<-subset(F5, NOdisturb==1)}
if(c==9){F2<-subset(F5, alldisturb==1)}
if(c==11){F2<-subset(F5harvest)}
if(c==10){F2<-subset(F5, allmort>0)}	

v<-c

if(v<7 )
	{
	f22<-F2	
	f23<-subset(F5, get(vari[v])< quantile(vari1[[v]], 0.7, na.rm=TRUE) & get(vari[v])> quantile(vari1[[v]], 0.3, na.rm=TRUE))
	F33<-subset(f23, LATLON %in% (f22$LATLON)) # not neccesaary bc weighting will eliminate latlons not in f22!
	}
if(v>=7 & v<11 )
	{
	if(c==2){next}
	f22<-subset(F2) 
	f23<-subset(F5, is.na(get(vari[v]))==T)
	F33<-subset(f23, LATLON %in% (f22$LATLON)) # not neccesaary bc weighting will eliminate latlons not in f22!
	}
if(v==11 )
	{
	f22<-subset(F2) 
	f23<-subset(F5harvest, cutting==0) #, is.na(get(vari[v]))==T)
	F33<-subset(f23, LATLON %in% (f22$LATLON)) # not neccesaary bc weighting will eliminate latlons not in f22!
	}
	
f22$Csamp<-(ave(f22$DIAmean>0,f22$LATLON, FUN=function(x) sum(x, na.rm=TRUE)))/f22$nsampTOT #weighting factor of number of Cond data/total data in each grid cell
F33$Nsamp<-(ave(F33$DIAmean>0,F33$LATLON, FUN=function(x) sum(x, na.rm=TRUE)))/F33$nsampTOT #weighting factor of number of NULL data/total data in each grid cell

f24<-subset(f22, !duplicated(LATLON), select=c(LATLON, Csamp))
F33<-merge(F33, f24, all.x=TRUE, by="LATLON")
F33$CtoNratio<- 1 #F33$Csamp/F33$Nsamp  #condition to null grid cell sampling ratio for weighted comparisons

F33$mortDIAmean <-F33$mortDIAmean*F33$CtoNratio

if(v<11)
{
f22b<-subset(f22, mortDIAmean>0) #
Vgrow<-mean(f22b$mortDIAmean, na.rm=T)*2.54 #
F33b<-subset(F33, mortDIAmean>0) #
Vgrow2<- mean(F33b$mortDIAmean, na.rm=T)/mean(F33b$CtoNratio, na.rm=T)*2.54 #
error2<-std.error(f22$mortDIAmean*2.54, na.rm=T)

}
if(v==11)
{
f22b<-subset(f22, cutDIAmean>0) 
Vgrow<-mean(f22b$cutDIAmean, na.rm=T)*2.54 #
F33b<-subset(F33, mortDIAmean>0) #
Vgrow2<-mean(F33b$mortDIAmean, na.rm=T)/mean(F33b$CtoNratio, na.rm=T)*2.54 
error2<-std.error(f22$cutDIAmean*2.54, na.rm=T)


}

Vgrow3<-Vgrow-Vgrow2 #mean(F2$mortDIAmean, na.rm=T)+mean(F2$cutDIAmean, na.rm=T)  #
n<-round(sum(f22b$mortality2/6.01)/sum(f22$PLT_CN>0), 1)*100

f34<-c(cond[c], n, round(Vgrow2,2), round(Vgrow,2), round(Vgrow3, 2), error2) #vari2[v2],
F5b<-rbind(F5b, f34)
print(c(c))
} 

F6b<-as.data.table(F5b)
setnames(F6b, c( "cond", "n","Vgrow2","Vdead", "difference", "error2")) #"Disturbance type",
print(na.omit(F6b))
#write.csv(na.omit(F6b))

Fbind4<-as.data.table(Fbind3)
setnames(Fbind4,c("x2", "y2mean", "TPA2", "endCerror", "endQMDerror", "cond", "Caccum", "meancond"))

FF<-as.data.table(merge(Fbind4, F6b, by="cond"))

FF<-as.data.table(merge(FF, Feq3, by="cond"))

FF2<-subset(FF) #, cond!="ct") 
y<-as.numeric(as.character(FF2$y2mean))
x<-as.numeric(as.character(FF2$difference))+as.numeric(as.character(FF2$Vdead[1]))
x[is.na(x)]<-as.numeric(as.character(FF2$Vdead[1]))
xerror<-(as.numeric(as.character(FF2$error2)))*2
yerror<-FF2$ySD*2
cond<-as.character(FF2$cond)


labels<-c( "all forests", expression(paste("-0.35"^"o","C cooling" )), "recurring\ndry summers", "disturbed",  expression(paste("-0.24 mm day"^"-1","precip." )), "all forests w/ logging", "all mortality", "recurring wet summers", expression(paste("+0.37 mm day"^"-1","precip." )), "undisturbed", expression(paste("+0.33"^"o","C warming" )))


#palette(colorRampPalette(c("black","black","black","black", "grey60", "grey60", "grey60", "grey60"))( 100)) ## (n)
dev.new(width=5, height=5) #width=6, height=5)
plot(x,y, type="n", main=c("Forest carbon equilibrium\n compared to tree morality size"),xlab=NA, ylab=NA, axes=FALSE, ylim=c(70,165), xlim=c(20, 28))

axis(side=1, tck=-0.01, labels=NA, lwd=0.75,  at=c( 20, 22, 24, 26,28, 28.5))
axis(side=2, tck=-0.01, labels=NA, lwd=0.75)
axis(side=1, lwd=0, line= -0.7, at=c(  20, 22, 24, 26,28))
axis(side=2, lwd=0, line= -0.7)

mtext(side=1, "Tree mortality mean diameter (cm)", line=1.2)
mtext(side=2, expression(paste("Carbon steady-state (Mg ha"^"-1",")")), line=1.2) 

segments(x-xerror, y, x+xerror, y, lwd=0.5, col="grey25")
segments(x, y-yerror, x, y+yerror, lwd=0.5, col="grey25")
abline(lm(y~x), col="grey5", lwd=0.95)
#points(x,y, pch=19)
points(x,y, pch=21, cex=0.99, bg= "black", lwd=0.25) #, col= c(av, ct, de, ds, dt, lg, lg2, pe, pt, ud, wt))
text(x+c(-0.1,0.1,-0.3,-0.2,0.3,-.25, 0.25,0.05,0.4,0,0.3),y+c(-0.25,1.25,-1.25,-1,0.75,0,0,-2.5,-1.25,0,0.25), labels, pos=c(2,4,2,2,4,2,3,2,4,4,4), col="grey25", cex=0.8, offset=0.4)

summary(lm(y~x))


#####################################################################################
#THE END






