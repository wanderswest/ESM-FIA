#EMPIRICAL SUCCESSION MAPPING Feb. 2015
#TRAVIS ANDREWS - TDA210@lehigh.edu

####CAUTION for() LOOPS AHEAD####


library(ncdf)
library(chron)
library(sp)
library(data.table)
library(plyr)
library(reshape2)
options(scipen=999) #remove sci notation

##open script in R directly   > source("/Volumes/m-z/tda210/USFS/Nullpaper2.R")

#####################################################################################################################################################################
############ RESURVEY GROWTH DATA PLOTS ###
#State_TREE_GRM_ESTN.CSV files must be downloaded from the FIA website (http://apps.fs.fed.us/fiadb-downloads/datamart.html) to a local drive
RMfiles<-list.files("/Users/robertbooth/Documents/TDA210/USFSFIA/GRM", pattern="*ESTN.CSV", full.names=TRUE) #get file names



S3RM<-data.table(NULL) #set up space for rbind of state data
t<-length(RMfiles)
print("script starting")  #for loop that opens state data, calculates variables, then adds it to data from previous states. 
for(i in 1:t)
{

     		RM <-subset(as.data.table(read.csv(RMfiles[i], header = TRUE, sep = ",", quote="\"", dec=".")), ESTN_TYPE=="AL" 
          	& LAND_BASIS=="TIMBERLAND", select=c(PLT_CN, INVYR, REMPER,TPAGROW_UNADJ,TPAREMV_UNADJ, 
			TPAMORT_UNADJ, ANN_NET_GROWTH, REMOVALS, MORTALITY, DIA_BEGIN, DIA_END, COMPONENT, DIA_BEGIN_RECALC)) 
		
			RM$SURVIVORTPA<- ifelse(RM$COMPONENT=="SURVIVOR", RM$TPAGROW_UNADJ,0) 
  	if(sum(RM$SURVIVORTPA)>1000) #states must have at least ~150 measured trees, else move to next state
		{ 
		RM$start<-RM$INVYR-RM$REMPER
		RM$CUT1TPA<- ifelse(RM$COMPONENT=="CUT1", RM$TPAGROW_UNADJ,0) 
		RM$CUT2TPA<- ifelse(RM$COMPONENT=="CUT2", RM$TPAGROW_UNADJ,0) 
		RM$CUT<-RM$CUT2TPA+RM$CUT1TPA
		RM$DIVERSION1TPA<- ifelse(RM$COMPONENT=="DIVERSION1", RM$TPAGROW_UNADJ,0) 
		RM$DIVERSION2TPA<- ifelse(RM$COMPONENT=="DIVERSION2", RM$TPAGROW_UNADJ,0) 
		RM$INGROWTHTPA<- ifelse(RM$COMPONENT=="INGROWTH", RM$TPAGROW_UNADJ,0) 
		RM$MORTALITY1TPA<- ifelse(RM$COMPONENT=="MORTALITY1", RM$TPAGROW_UNADJ,0) 
		RM$MORTALITY2TPA<- ifelse(RM$COMPONENT=="MORTALITY2", RM$TPAGROW_UNADJ,0)
		RM$MORTALITYTPA<-RM$MORTALITY1TPA+RM$MORTALITY2TPA
		RM$REVERSION1TPA<- ifelse(RM$COMPONENT=="REVERSION1", RM$TPAGROW_UNADJ,0) 
		RM$REVERSION2TPA<- ifelse(RM$COMPONENT=="REVERSION2", RM$TPAGROW_UNADJ,0) 
		RM$REDIV<-RM$REVERSION2TPA+RM$REVERSION1TPA+RM$DIVERSION2TPA+RM$DIVERSION1TPA
		RM$REDIV<-ave(RM$REDIV, RM$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE))
		RM$cut<-ave(RM$CUT, RM$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE))

		RM<-subset(RM, REDIV==0) #exclude plots with reversion or diversion
		RM$VolCNG<-ave(RM$ANN_NET_GROWTH, RM$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE)) #ANN_NET_GROWTH negitive in mortality and partial damaged
		RM$VolCNGmean<-ave(RM$ANN_NET_GROWTH, RM$PLT_CN, FUN=function(x) mean(x, na.rm=TRUE))
		RM$VolCNGmax<-ave(RM$ANN_NET_GROWTH, RM$PLT_CN, FUN=function(x) max(x, na.rm=TRUE))
		RM$VolCNGmin<-ave(RM$ANN_NET_GROWTH, RM$PLT_CN, FUN=function(x) min(x, na.rm=TRUE))
		RM$VolCNGmid<-ave(RM$ANN_NET_GROWTH, RM$PLT_CN, FUN=function(x) median(x, na.rm=TRUE))
		RM$ANN_GROW<-ifelse(RM$ANN_NET_GROWTH>=0,RM$ANN_NET_GROWTH,NA)
		RM$VolCNGgrow<-ave(RM$ANN_GROW, RM$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE))
    	RM$ANN_DIE<-ifelse(RM$ANN_NET_GROWTH<=0,RM$ANN_NET_GROWTH,NA)
    	RM$VolCNGdie<-ave(RM$ANN_DIE, RM$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE))
    	
    	RM$start1tpa<-RM$SURVIVORTPA+RM$MORTALITYTPA #initial number of trees is current survivors plus those that died during the resurvey period.
		RM$startsumTPA<-ave(RM$start1tpa, RM$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE)) 
		RM$end1tpa<-RM$SURVIVORTPA+RM$INGROWTHTPA #final number of trees is current survivors plus new trees that cross the 5 inch measurement threshold
		RM$endsumTPA<-ave(RM$end1tpa, RM$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE))

		##Re-measure plot mean DBH begining and end calcs: NOTE ingrowth/reversion and mortality/diversion do not have DIA_BEGIN or DIA_END respectively	
		RM$DIAbeginmean<-ave(RM$DIA_BEGIN, RM$PLT_CN, FUN=function(x) mean(x, na.rm=TRUE)) 
		RM$DIAendmean<-ave(RM$DIA_END, RM$PLT_CN, FUN=function(x) mean(x, na.rm=TRUE))
		RM$RADbeginmean<-ave(RM$DIA_BEGIN/2, RM$PLT_CN, FUN=function(x) mean(x, na.rm=TRUE)) 
		RM$RADendmean<-ave(RM$DIA_END/2, RM$PLT_CN, FUN=function(x) mean(x, na.rm=TRUE))

		RM$DIAgrow<-RM$DIA_END-RM$DIA_BEGIN
   		RM$DIAgrowmean<-ave(RM$DIAgrow, RM$PLT_CN, FUN=function(x) mean(x, na.rm=TRUE))
        RM$DIAgrow10<-ifelse(RM$DIA_BEGIN>10, RM$DIAgrow, NA)
   		RM$DIAgrow10mean<-ave(RM$DIAgrow10, RM$PLT_CN, FUN=function(x) mean(x, na.rm=TRUE))
        RM$DIAgrow15<-ifelse(RM$DIA_BEGIN>15, RM$DIAgrow, NA)
   		RM$DIAgrow15mean<-ave(RM$DIAgrow15, RM$PLT_CN, FUN=function(x) mean(x, na.rm=TRUE))
		RM$ingrowth2<-ave(RM$INGROWTHTPA, RM$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE))
		RM$mortality2<-ave(RM$MORTALITYTPA, RM$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE))
		RM$mortDIA<-ifelse(RM$MORTALITYTPA>0, RM$DIA_BEGIN, NA)
		RM$mortDIAmean<-ave(RM$mortDIA, RM$PLT_CN, FUN=function(x) mean(x, na.rm=TRUE))
		RM$mortDIAmax<-ave(RM$mortDIA, RM$PLT_CN, FUN=function(x) max(x, na.rm=TRUE))
		RM$cutDIA<-ifelse(RM$CUT>0, RM$DIA_BEGIN, NA)
		RM$cutDIAmean<-ave(RM$cutDIA, RM$PLT_CN, FUN=function(x) mean(x, na.rm=TRUE))
		RM$cutDIAmax<-ave(RM$cutDIA, RM$PLT_CN, FUN=function(x) max(x, na.rm=TRUE))

	RM1<-subset(RM, !duplicated(PLT_CN))  #simplify datatable to averaged variables for each plot in one row       
		
RM2<-subset(RM1, start>0, select=c(PLT_CN, INVYR, REMPER ,startsumTPA, VolCNG, VolCNGmean, VolCNGmin, VolCNGmax, VolCNGmid, VolCNGgrow, VolCNGdie, endsumTPA,  ingrowth2, mortality2, DIAbeginmean, DIAendmean, RADbeginmean, RADendmean, DIAgrowmean,DIAgrow10mean, DIAgrow15mean, cut, mortDIAmean, mortDIAmax, cutDIAmean, cutDIAmax ))

print(RMfiles[i])
print(Sys.time())
S3RM<-rbind(S3RM, RM2)
		}
}
RM5<-S3RM
#write.csv(RM5, file = "/Volumes/m-z/tda210/USFS/RM.9.10.2014.csv") #save data



##############DATA#############################DATA##############################DATA#################################

treefiles<-list.files("/Users/robertbooth/Documents/TDA210/USFSFIA", pattern="*TREE.CSV", full.names=TRUE)
condfiles<-list.files("/Users/robertbooth/Documents/TDA210/USFSFIA", pattern="*COND.CSV", full.names=TRUE)
plotfiles<-list.files("/Users/robertbooth/Documents/TDA210/USFSFIA", pattern="*PLOT.CSV", full.names=TRUE)

S3<-data.table(NULL)
t<-length(treefiles)

#load remeasured tree data 
#RM5<-as.data.table(read.csv("/Volumes/m-z/tda210/USFS/RM.9.10.2014.csv", header = TRUE, sep = ",", quote="\"", dec="."))

for(i in 1:t)
{
	print(c("Starting:", treefiles[i]))
	PLOT <- subset(as.data.table(read.csv(plotfiles[i],
                 header = TRUE, sep = ",", quote="\"", dec=".")), LON>(-95), select=c(CN, LAT, LON,
                 ELEV, ECOSUBCD, STATECD,PREV_PLT_CN, REMPER, DESIGNCD, KINDCD, MICROPLOT_LOC))
            setnames(PLOT, "CN", "PLT_CN")

	if(sum(PLOT$REMPER, na.rm=TRUE)>10) #proceed if state has been resurveyed
	{	
	TREE <- subset(as.data.table(read.csv(treefiles[i], header = TRUE, sep = ",", quote="\"", dec=".")), select=c(CN, PREV_TRE_CN, PLT_CN, INVYR,CONDID, DIA, TPA_UNADJ, HT, AGENTCD, DAMTYP1, DECAYCD, SPCD, STOCKING, STATUSCD, STANDING_DEAD_CD, SPGRPCD, CR, CCLCD, RECONCILECD, DIACHECK, DIAHTCD, PREVDIA, PREV_STATUS_CD, P2A_GRM_FLG, CARBON_AG, CARBON_BG, PLOT, TREE))

	COND <- subset(as.data.table(read.csv(condfiles[i], header = TRUE, sep = ",", quote="\"", dec=".")), select=c(PLT_CN,CONDID, FORINDCD, STDAGE, STDSZCD, FLDSZCD, FLDTYPCD,SITECLCD, SLOPE, PHYSCLCD, GSSTKCD,ALSTK,ALSTKCD, GSSTK, DSTRBCD1,DSTRBYR1,DSTRBCD2,DSTRBYR2, TRTCD1, TRTYR1,TRTCD2,TRTYR2,TRTCD3,TRTYR3,PRESNFCD, DSTRBCD3, DSTRBYR3, BALIVE, CARBON_DOWN_DEAD,STDORGCD, SISP))

	COND$CONmax<-ave(COND$CONDID, COND$PLT_CN,FUN=function(x) max(x, na.rm=TRUE)) 
	COND<-subset(COND, CONmax<2)  ## remove all plots with more than 1 condition 
	PC<-merge(COND,PLOT,all.x=TRUE, by=c("PLT_CN"))
	TREE$CONmax<-ave(TREE$CONDID, TREE$PLT_CN,FUN=function(x) max(x, na.rm=TRUE))
	TREE$RECONCILECD[is.na(TREE$RECONCILECD)]<-0
	TREE$PREVDIAna<-(TREE$PREVDIA)
	TREE$PREVDIAna[is.na(TREE$PREVDIAna)]<-0
	TREE$PREVSTATna<-TREE$PREV_STATUS_CD
	TREE$PREVSTATna[is.na(TREE$PREVSTATna)]<-0
	TREE$PLTnb<-ifelse(TREE$STATUSCD==1, 1,0)
	TREE$PLTn<-ave(TREE$PLTnb, TREE$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE))
	TREE$PrevPLTnb<-ifelse(TREE$PREV_STATUS_CD==1, 1,0)
	TREE$PrevPLTn<-ave(TREE$PrevPLTnb, TREE$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE))
	TREE$STATUSCDmax<-ifelse(TREE$STATUSCD==3, 3,0)
	TREE$STATUSCDmax<-ave(TREE$STATUSCDmax, TREE$PLT_CN, FUN=function(x) sum(x, na.rm=TRUE))
	
	TREE1<-subset(TREE, CONmax<2 & INVYR<2014) #remove edge effects 
	TCPb<-merge(TREE1, PC, by=c("PLT_CN")) 
	TCPb<-subset(TCPb) #, DESIGNCD==1)
	
	#connect PREV_CN for each tree prior to subset
	TREE3<-subset(TCPb, select=c(CN, DIA, HT, STOCKING, SPCD, TPA_UNADJ, AGENTCD, DSTRBCD1, STDAGE, MICROPLOT_LOC, SISP, CARBON_AG, CARBON_BG))
	PREVnames<-names(TREE3)
	PREVnamesnew<-paste("PREV_TRE_",PREVnames,sep="")
	setnames(TREE3, PREVnames,PREVnamesnew) 
	TREE3$PREV_TRE_CN<-as.numeric((TREE3$PREV_TRE_CN))
	TCPb$PREV_TRE_CN<-as.numeric((TCPb$PREV_TRE_CN))
	TCP<-merge(TREE3,TCPb,  all.x=TRUE, all.y=TRUE, by="PREV_TRE_CN") #all.x=TRUE, all.y=TRUE
	TCP<-subset(TCP, REMPER>3 & REMPER<9.5 & STATUSCDmax!=3 & STATUSCD!=0 & DESIGNCD==1) #remove plots that have not been resurveyed and those surveyed over unusually long or short periods or were harvested,  STATUSCD!=0:remove trees not resurveyed due to RECONCILECD = 5-9
	
	TCP$PREV_TRE_STDAGE<-ifelse(is.na(TCP$PREV_TRE_STDAGE)==TRUE, TCP$STDAGE-TCP$REMPER, TCP$PREV_TRE_STDAGE) #fill in missing with current age minus remeasurement period
	
	##DIAmean of DIA>5
	TCP$DIA5alive<-ifelse(TCP$DIA>=5 & TCP$STATUSCD==1 & TCP$P2A_GRM_FLG!="N", TCP$DIA, NA)
	TCP$DIA5meanalive<-ave(TCP$DIA5alive,TCP$PLT_CN,FUN=function(x) mean(x,na.rm=TRUE))
	TCP$PREVDIA5alive<-ifelse(TCP$PREVDIA>=5 & TCP$PREV_STATUS_CD==1 & TCP$P2A_GRM_FLG!="N", TCP$PREVDIA, NA)
	TCP$PREVDIA5meanalive<-ave(TCP$PREVDIA5alive,TCP$PLT_CN,FUN=function(x) mean(x,na.rm=TRUE))	
	TCP$BA5meanalive<-ave(((TCP$DIA5alive/2)^2)*pi,TCP$PLT_CN,FUN=function(x) mean(x,na.rm=TRUE))
	TCP$PREVBA5meanalive<-ave(((TCP$PREVDIA5alive/2)^2)*pi,TCP$PLT_CN,FUN=function(x) mean(x,na.rm=TRUE))
	
	#Trees per acre (TPA) for DIA>5
	TCP$TPA_UNADJ5<-ifelse(TCP$DIA5alive>=5, TCP$TPA_UNADJ, NA)
	TCP$TPAsum5<-ave(TCP$TPA_UNADJ5,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
	TCP$PREVTPA_UNADJ5<-ifelse(TCP$PREVDIA5alive>=5, TCP$PREV_TRE_TPA_UNADJ, NA)
	TCP$PREVTPAsum5<-ave(TCP$PREVTPA_UNADJ5,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
#Trees per acre (TPA) for DIA>10
	TCP$TPA_UNADJ10<-ifelse(TCP$DIA5alive>10, TCP$TPA_UNADJ, NA)
	TCP$TPAsum10<-ave(TCP$TPA_UNADJ10,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
	TCP$PREVTPA_UNADJ10<-ifelse(TCP$PREVDIA5alive>10, TCP$PREV_TRE_TPA_UNADJ, NA)
	TCP$PREVTPAsum10<-ave(TCP$PREVTPA_UNADJ10,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
#Trees per acre (TPA) for DIA>5 & DIA<10
	TCP$TPA_UNADJ510<-ifelse(TCP$DIA5alive>=5 & TCP$DIA5alive<=10, TCP$TPA_UNADJ, NA)
	TCP$TPAsum510<-ave(TCP$TPA_UNADJ510,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
	TCP$PREVTPA_UNADJ510<-ifelse(TCP$PREVDIA5alive>=5 & TCP$PREVDIA5alive<=10, TCP$PREV_TRE_TPA_UNADJ, NA)
	TCP$PREVTPAsum510<-ave(TCP$PREVTPA_UNADJ510,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 

#Trees per acre (TPA) for DIA>15
	TCP$TPA_UNADJ15<-ifelse(TCP$DIA5alive>=15, TCP$TPA_UNADJ, NA)
	TCP$TPAsum15<-ave(TCP$TPA_UNADJ15,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
	TCP$PREVTPA_UNADJ15<-ifelse(TCP$PREVDIA5alive>=15, TCP$PREV_TRE_TPA_UNADJ, NA)
	TCP$PREVTPAsum15<-ave(TCP$PREVTPA_UNADJ15,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 

#Trees per acre (TPA) for DIA<5 
	TCP$TPA_UNADJ05<-ifelse(TCP$PREVDIA<5 & TCP$STATUSCD==1 , TCP$TPA_UNADJ, NA) #use prev dia to account for saplings that cross 5inch threshold
	TCP$TPA_UNADJ05<-ifelse(TCP$DIA<5 & TCP$STATUSCD==1 , TCP$TPA_UNADJ, TCP$TPA_UNADJ05) #account for new saplings ingrowth
	TCP$TPAsum05<-ave(TCP$TPA_UNADJ05,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
	TCP$PREVTPA_UNADJ05<-ifelse(TCP$PREVDIA<5 & TCP$PREV_STATUS_CD==1, TCP$PREV_TRE_TPA_UNADJ, NA)
	TCP$PREVTPAsum05<-ave(TCP$PREVTPA_UNADJ05,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
	TCP$sapMort<-ifelse(TCP$PREVDIA<5 & TCP$PREV_STATUS_CD==1 & TCP$STATUSCD==2 , 1, NA) 
	TCP$sapMortsum<-ave(TCP$sapMort,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
	TCP$sapRecruit<-ifelse(TCP$DIA<5 & is.na(TCP$PREV_STATUS_CD)==T & TCP$STATUSCD==1 , 1, NA) 
	TCP$sapRecruitsum<-ave(TCP$sapRecruit,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
	
	#Stocking of plots for trees with DIA>5
	TCP$STOCKING5<-ifelse(TCP$DIA5alive>0, TCP$STOCKING, NA)
	TCP$STOCKING5mid<-ave(TCP$STOCKING5,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 	
	TCP$PREVSTOCKING5<-ifelse(TCP$PREVDIA5alive>0, TCP$PREV_TRE_STOCKING, NA)
	TCP$PREVSTOCKING5mid<-ave(TCP$PREVSTOCKING5,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 

	#Aboveground and belowground Carbon (excluding foliage) for DIA>5 live trees  per acre then dead trees	
	TCP$carbon5<-ifelse(TCP$STATUSCD==1 & TCP$DIA>=5 & TCP$P2A_GRM_FLG!="N", TCP$CARBON_AG*TCP$TPA_UNADJ+TCP$CARBON_BG*TCP$TPA_UNADJ, NA)
	TCP$carbon5sum<-ave(TCP$carbon5,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))
	TCP$PREVcarbon5<-ifelse(TCP$PREV_STATUS_CD==1 & TCP$PREVDIA>=5 & TCP$P2A_GRM_FLG!="N", TCP$PREV_TRE_CARBON_AG*TCP$PREV_TRE_TPA_UNADJ +TCP$PREV_TRE_CARBON_BG*TCP$PREV_TRE_TPA_UNADJ, NA)
	TCP$PREVcarbon5sum<-ave(TCP$PREVcarbon5,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))
	
	TCP$carbon1<-ifelse(TCP$STATUSCD==1 & TCP$P2A_GRM_FLG!="N", TCP$CARBON_AG*TCP$TPA_UNADJ+TCP$CARBON_BG*TCP$TPA_UNADJ, NA)
	TCP$carbon1sum<-ave(TCP$carbon1,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))
	TCP$PREVcarbon1<-ifelse(TCP$PREV_STATUS_CD==1 &TCP$P2A_GRM_FLG!="N", TCP$PREV_TRE_CARBON_AG*TCP$PREV_TRE_TPA_UNADJ +TCP$PREV_TRE_CARBON_BG*TCP$PREV_TRE_TPA_UNADJ, NA)
	TCP$PREVcarbon1sum<-ave(TCP$PREVcarbon1,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))

	TCP$deadcarbon5<-ifelse(TCP$PREV_STATUS_CD==1 & TCP$STATUSCD==2 , TCP$CARBON_AG*TCP$TPA_UNADJ+TCP$CARBON_BG*TCP$TPA_UNADJ, NA)
	TCP$deadcarbon5sumYR<-(ave(TCP$deadcarbon5,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)))/TCP$REMPER
	TCP$deadcarbon5T<-ifelse(TCP$PREV_STATUS_CD==1 & TCP$STATUSCD==2 , TCP$CARBON_AG+TCP$CARBON_BG, NA)
	TCP$deadcarbon5TREE<-(ave(TCP$deadcarbon5T,TCP$PLT_CN,FUN=function(x) mean(x,na.rm=TRUE)))
	TCP$deadcarbon5CT<-ifelse(TCP$PREV_STATUS_CD==1 & TCP$STATUSCD==2 & TCP$AGENTCD==80, TCP$CARBON_AG+TCP$CARBON_BG, NA)
	TCP$deadcarbon5CUTTREE<-(ave(TCP$deadcarbon5CT,TCP$PLT_CN,FUN=function(x) mean(x,na.rm=TRUE)))

#disturbances
	TCP$DEADsum<-ave(TCP$STANDING_DEAD_CD,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) 
	TCP$insect<-ifelse(TCP$AGENTCD==10 & TCP$PREVDIA>=5 & TCP$PREV_STATUS_CD==1,1,0)
	TCP$insect<-ave(TCP$insect,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))
	TCP$PREV_insect<-ifelse(TCP$PREV_TRE_AGENTCD==10,1,0)
	TCP$PREV_insect<-ave(TCP$PREV_insect,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))
	TCP$disease<-ifelse(TCP$AGENTCD==20 & TCP$PREVDIA>=5 & TCP$PREV_STATUS_CD==1,1,0)
	TCP$disease<-ave(TCP$disease,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))	
	TCP$PREV_disease<-ifelse(TCP$PREV_TRE_AGENTCD==20,1,0)
	TCP$PREV_disease<-ave(TCP$PREV_disease,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))	
	TCP$fire<-ifelse(TCP$AGENTCD==30 & TCP$PREVDIA>=5 & TCP$PREV_STATUS_CD==1,1,0)
	TCP$fire <-ave(TCP$fire,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))	
	TCP$PREV_fire<-ifelse(TCP$PREV_TRE_AGENTCD==30,1,0)
	TCP$PREV_fire <-ave(TCP$PREV_fire,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))	
	TCP$animal<-ifelse(TCP$AGENTCD==40 & TCP$PREVDIA>=5 &TCP$PREV_STATUS_CD==1,1,0)
	TCP$animal <-ave(TCP$animal,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))
	TCP$PREV_animal<-ifelse(TCP$PREV_TRE_AGENTCD==40,1,0)
	TCP$PREV_animal <-ave(TCP$PREV_animal,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))		
	TCP$weather<-ifelse(TCP$AGENTCD==50 & TCP$PREVDIA>=5 &TCP$PREV_STATUS_CD==1,1,0)
	TCP$weather<-ave(TCP$weather,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))	
	TCP$PREV_weather<-ifelse(TCP$PREV_TRE_AGENTCD==50,1,0)
	TCP$PREV_weather<-ave(TCP$PREV_weather,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))	
	TCP$vegetation<-ifelse(TCP$AGENTCD==60 & TCP$PREVDIA>=5 &TCP$PREV_STATUS_CD==1,1,0) #suppression, competition, vines
	TCP$vegetation <-ave(TCP$vegetation,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))	
	TCP$PREV_vegetation<-ifelse(TCP$PREV_TRE_AGENTCD==60,1,0) #suppression, competition, vines
	TCP$PREV_vegetation <-ave(TCP$PREV_vegetation,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))	
	TCP$unknowndamage<-ifelse(TCP$AGENTCD==70 & TCP$PREVDIA>=5 &TCP$PREV_STATUS_CD==1,1,0)
	TCP$unknowndamage <-ave(TCP$unknowndamage,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))	
	TCP$PREV_unknowndamage<-ifelse(TCP$PREV_TRE_AGENTCD==70,1,0)
	TCP$PREV_unknowndamage <-ave(TCP$PREV_unknowndamage,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))	
	TCP$cutting<-ifelse(TCP$AGENTCD==80 & TCP$PREVDIA>=5, 1,0)
	TCP$cutting <-ave(TCP$cutting,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE)) #should be zero!
	TCP$PREV_cutting<-ifelse(TCP$PREV_TRE_AGENTCD==80,1,0)
	TCP$PREV_cutting <-ave(TCP$PREV_cutting,TCP$PLT_CN,FUN=function(x) sum(x,na.rm=TRUE))
	TCP$SPCDdiversity<-ave(TCP$SPCD,TCP$PLT_CN,FUN=function(x) length(unique(x)))
	TCP$PREV_SPCDdiversity<-ave(TCP$PREV_TRE_SPCD,TCP$PLT_CN,FUN=function(x) length(unique(x)))

##disturbance categories
	TCP$Undisturbed<-ifelse((TCP$cutting==0 & TCP$PREV_cutting==0 & TCP$PREV_TRE_DSTRBCD1<1 & TCP$DSTRBCD1<5  & TCP$PREV_insect==0 & TCP$PREV_fire==0 & TCP$PREV_disease==0 & TCP$PREV_weather==0 & TCP$PREV_vegetation==0 & TCP$PREV_unknowndamage==0  & TCP$insect==0 & TCP$fire==0  & TCP$weather==0 & TCP$vegetation==0 & TCP$unknowndamage==0 & TCP$disease==0 & TCP$STDAGE>TCP$PREV_TRE_STDAGE), 1, NA)

	TCP$UndisturbedNOVEG<-ifelse((TCP$cutting==0 & TCP$PREV_cutting==0 & TCP$PREV_TRE_DSTRBCD1<1 & TCP$DSTRBCD1<5 & TCP$PREV_insect==0 & TCP$PREV_fire==0 & TCP$PREV_disease==0 & TCP$PREV_weather==0  & TCP$PREV_unknowndamage==0  & TCP$insect==0 & TCP$fire==0  & TCP$weather==0 & TCP$unknowndamage==0 & TCP$disease==0 & TCP$STDAGE>TCP$PREV_TRE_STDAGE), 1, NA)

	TCP$UndisturbedNOVEGnoPREV<-ifelse((TCP$cutting==0 &  TCP$DSTRBCD1<5 & TCP$insect==0 & TCP$fire==0  & TCP$weather==0 & TCP$unknowndamage==0 & TCP$disease==0 & TCP$animal==0), 1, NA)

	TCP$UndisturbedNOVEGnoPREVnoUN<-ifelse((TCP$cutting==0 & TCP$animal==0 & TCP$insect==0 & TCP$fire==0  & TCP$weather==0 & TCP$disease==0), 1, NA) #& TCP$DSTRBCD1<5 &

	TCP$minor.disturb<-ifelse((TCP$cutting==0 & TCP$PREV_cutting==0 & TCP$STDAGE>TCP$PREV_TRE_STDAGE & TCP$PREV_TRE_DSTRBCD1>1 | TCP$DSTRBCD1>1 | TCP$PREV_insect>0 & TCP$PREV_insect<5 | TCP$PREV_fire>0 & TCP$PREV_fire<5 | TCP$PREV_disease>0 & TCP$PREV_disease<5 | TCP$PREV_weather>0 & TCP$PREV_weather<5 | TCP$PREV_unknowndamage>0  & TCP$PREV_unknowndamage<5 | TCP$insect>0 & TCP$insect<5 | TCP$fire>0 & TCP$fire<5 | TCP$disease>0 & TCP$disease<5 | TCP$weather>0 & TCP$weather<5 | TCP$unknowndamage>0  & TCP$unknowndamage<5), 1, NA)

	TCP$major.disturb<-ifelse((TCP$cutting>0 | TCP$PREV_cutting>0 | TCP$STDAGE<TCP$PREV_TRE_STDAGE | TCP$PREV_TRE_DSTRBCD1>1 | TCP$DSTRBCD1>1 | TCP$PREV_insect>=5 | TCP$PREV_fire>=5 | TCP$PREV_disease>=5 | TCP$PREV_weather>=5 | TCP$PREV_unknowndamage>=5 | TCP$insect>=5 | TCP$fire>=5 | TCP$disease>=5 | TCP$weather>=5 | TCP$unknowndamage>=5), 1, NA)

	TCP$major.prev.disturb<-ifelse((TCP$cutting==0  & TCP$fire==0 & TCP$disease==0 & TCP$weather==0 & TCP$unknowndamage==0 & TCP$insect==0 & TCP$DSTRBCD1<1 & TCP$PREV_cutting>5 | TCP$PREV_TRE_DSTRBCD1>1 | TCP$PREV_insect>=5 | TCP$PREV_fire>=5 | TCP$PREV_disease>=5 | TCP$PREV_weather>=5 |  TCP$PREV_unknowndamage>=5 ), 1, NA)

	TCP$minor.prev.disturb<-ifelse((TCP$cutting==0 & TCP$PREV_cutting==0 & TCP$STDAGE>TCP$PREV_TRE_STDAGE & TCP$PREV_TRE_DSTRBCD1>1 |    TCP$PREV_insect>0 & TCP$PREV_insect<5 | TCP$PREV_fire>0 & TCP$PREV_fire<5 | TCP$PREV_disease>0 & TCP$PREV_disease<5 | TCP$PREV_weather>0 & TCP$PREV_weather<5 | TCP$PREV_unknowndamage>0  & TCP$PREV_unknowndamage<5 & TCP$insect==0 & TCP$fire==0 & TCP$disease==0 & TCP$weather==0 & TCP$unknowndamage==0 &TCP$DSTRBCD1==0 ), 1, NA)
	
	S1<-NULL
S1<-subset(TCP, !duplicated(PLT_CN) & PREV_PLT_CN>0, select=c(INVYR, PLT_CN, STOCKING5mid, PREVSTOCKING5mid, DIA5meanalive, PREVDIA5meanalive, BA5meanalive, PREVBA5meanalive, TPAsum5, PREVTPAsum5, TPAsum10, PREVTPAsum10, TPAsum510, PREVTPAsum510, TPAsum05, PREVTPAsum05, sapMortsum, sapRecruitsum, TPAsum15, PREVTPAsum15, LAT, LON, ELEV, ECOSUBCD, STATECD, PREV_PLT_CN, REMPER, STDAGE, PREV_TRE_STDAGE, SLOPE, PHYSCLCD, BALIVE, CARBON_DOWN_DEAD, STDORGCD, SISP, PREV_TRE_SISP, carbon5sum, PREVcarbon5sum, carbon1sum, PREVcarbon1sum, deadcarbon5sumYR, deadcarbon5TREE, deadcarbon5CUTTREE, UndisturbedNOVEG, minor.disturb, major.disturb, major.prev.disturb, minor.prev.disturb, DEADsum, SPCDdiversity, PREV_SPCDdiversity, UndisturbedNOVEGnoPREVnoUN, UndisturbedNOVEGnoPREV, animal, insect, fire, weather, disease, cutting, PREV_cutting, unknowndamage, vegetation, DSTRBCD1, PREV_TRE_DSTRBCD1))

print("S1 done")		
S3<-rbind(S3, S1)
	}		
print(treefiles[i])
print(Sys.time())
}




##############ENVIRONMENTAL DATA############

#######################################################################
a<-b<-2 #0.5 degree grid
print("starting environmental data")

# A note from the data providers: GHCN Gridded V2 data provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, from their Web site at http://www.esrl.noaa.gov/psd/ in any documents or publications using these data. We would also appreciate receiving a copy of the relevant publications. This will help PSD to justify keeping the GHCN Gridded V2 data set freely available online in the future. Thank you!


#######################################################################################
####### SOIL MOISTURE http://www.esrl.noaa.gov/psd/data/gridded/data.cpcsoil.html
##  - lat lon centered 0.5 grid. (i.e. bottom right grid corner -.25N and +.25E)
nc<-open.ncdf("/Users/robertbooth/Documents/TDA210/USFSFIA/soilw.mon.mean.V2.nc")
	# get latitudes and longitudes for NORTH AMERICA
	lat <- get.var.ncdf(nc,"lat",verbose=F)
	lat<-lat[70:129]
	nlat <- dim(lat)
	lon <- get.var.ncdf(nc,"lon")
	lon<-lon[360:719]
	nlon <- dim(lon)
	# get the data and attributes North AMERICA
	soilw.array <- get.var.ncdf(nc,"soilw",start= c(360,70,1), count= c(360,60,780)) 
	close.ncdf(nc)
	soilw.vec.long <- as.vector(soilw.array)
	soilw.mat <- matrix(soilw.vec.long, nrow=nlon*nlat*12, ncol=65, byrow=FALSE) #nrow=nlon*nlat*65.25
	lonlatyr <- expand.grid(lon,lat,1:12) #1948:2013),
	soilw.df02 <- data.table(cbind(lonlatyr,soilw.mat))
	setnames(soilw.df02,c("lon","lat","Month", 1948:2012)) 

SWsub1<-subset(soilw.df02, lon>230 & lon<294 & Month>5 & Month<9)
SWsub.ml1<-reshape(SWsub1, direction="long", varying=list(names(SWsub1)[4:68]), v.names="soilw",idvar=c("lon","lat","Month"), timevar="Year",times=1948:2012) 
SWsub.ml<-subset(SWsub.ml1, Year>=1985)

#SWsub.ml$lat<-(SWsub.ml$lat-0.25) #; ADDed a for Lat adjustment
#SWsub.ml$lon<-(SWsub.ml$lon+0.25)-360 #;  ADDed b for lon adj.

SWsub.ml$iLAT<-SWsub.ml$lat-0.25
SWsub.ml$iLON<-(SWsub.ml$lon+0.25)-360 
SWsub.ml$LATLONYR<-floor((SWsub.ml$iLAT* 10000000000000)+(SWsub.ml$iLON*(-100000000)) + (SWsub.ml$Year) + (SWsub.ml$Month/10)) #water year
SWsub.ml$LATLON<-as.numeric((SWsub.ml$iLAT*100000000)+(SWsub.ml$iLON*-1000)) 
SWsub.ml$soilw<-ave(SWsub.ml$soilw, SWsub.ml$LATLONYR, FUN=function(x) mean(x, na.rm=TRUE))  

SWsub.final<-subset(SWsub.ml, !duplicated(LATLONYR), select=c("soilw", "LATLONYR", "Year", "iLAT", "iLON", "LATLON"))
SWsub.final$soilwZ<-ave(SWsub.final$soilw, SWsub.final$LATLON, FUN=function(x) scale(x))
SWsub.final$soilwZ[is.na(SWsub.final$soilwZ)]<-0



###new stufff 

SWsub.final$SWdrought1<-ifelse(SWsub.final$soilwZ<0, SWsub.final$soilwZ, 0)   #  ave(SWsub.final$soilwZ, SWsub.final$LATLONYR, FUN=function(x) min(x, na.rm=TRUE))  
SWsub.final$SWpluvial1<-ifelse(SWsub.final$soilwZ>0, SWsub.final$soilwZ, 0)   # ave(SWsub.final$soilwZ, SWsub.final$LATLONYR, FUN=function(x) max(x, na.rm=TRUE))  

SWsub<-subset(SWsub.final, !duplicated(LATLONYR), Year>=1985, select=c("SWdrought1", "SWpluvial1", "LATLONYR", "Year", "iLAT", "iLON"))
SWfinal<-SWsub
SWfinal$SWdrought5<-0 #SWfinal$SWdrought1 #set resurvey year = no effect 
SWfinal$SWpluvial5<-0 #SWfinal$SWpluvial1

for(i in 1:4)
	{
	SWsub2<-subset(SWsub, !duplicated(LATLONYR), select=c("SWdrought1", "SWpluvial1", "LATLONYR"))
	SWsub2$LATLONYR<-SWsub2$LATLONYR+i
	setnames(SWsub2, c("SWdrought1","SWpluvial1"), c("SWdrought2","SWpluvial2"))
	SWfinal<-merge(SWfinal, SWsub2, by="LATLONYR")
	SWfinal$SWdrought5<-ifelse(SWfinal$SWdrought5<=SWfinal$SWdrought2, SWfinal$SWdrought5, SWfinal$SWdrought2)
	SWfinal$SWpluvial5<-ifelse(SWfinal$SWpluvial5>=SWfinal$SWpluvial2,SWfinal$SWpluvial5, SWfinal$SWpluvial2)
	SWfinal<-subset(SWfinal, select=-c(SWdrought2, SWpluvial2))
	}

SWfinal<-subset(SWfinal, !duplicated(LATLONYR) & Year>=1991, select=c("SWdrought5", "SWpluvial5", "LATLONYR"))
SW2<-merge(SWsub.final, SWfinal, by="LATLONYR")


#################### 



#SW2<-subset(SWsub.final, Year>=1991)
#SW2$SWdroughtminLL<-ave(SW2$soilwZ, SW2$LATLON, FUN=min)

SW2$SWdrought<-ifelse(SW2$soilwZ<(-1.5), SW2$soilwZ, NA)
SW2$droughtYRLL<-ave(SW2$SWdrought, SW2$LATLON, FUN=function(x) min(x,na.rm=TRUE)) #finding worst drought in latlon
SW2$droughtYR<-ifelse(SW2$soilwZ==SW2$droughtYRLL, SW2$Year, NA)
SW2$droughtYRLL2<-ave(SW2$droughtYR, SW2$LATLON, FUN=function(x) mean(x,na.rm=TRUE)) 
SW2$droughtYRLL2[is.na(SW2$droughtYRLL2)]<-0

SW2$SWpluvial<-ifelse(SW2$soilwZ>1, SW2$soilwZ, NA) #0.8
SW2$pluvialYRLL<-ave(SW2$SWpluvial, SW2$LATLON, FUN=function(x) max(x,na.rm=TRUE)) 
SW2$pluvialYR<-ifelse(SW2$soilwZ==SW2$pluvialYRLL, SW2$Year, NA)
SW2$pluvialYRLL2<-ave(SW2$pluvialYR, SW2$LATLON, FUN=function(x) mean(x,na.rm=TRUE)) 
SW2$pluvialYRLL2[is.na(SW2$pluvialYRLL2)]<-0

SW3<-subset(SW2, select=c(soilw, soilwZ, droughtYRLL2, pluvialYRLL2, droughtYRLL, pluvialYRLL,LATLONYR, LATLON, Year, SWdrought5, SWpluvial5)) #
SW3$LLn<-ave(SW3$soilwZ, SW3$LATLON,FUN=function(x) length(x)) 

#################### Soil Water TRENDS
T8b<-data.table(NULL)
T8b<-subset(SW3, soilw>0 , select=c(Year, soilw, LATLON)) #TPAsumSapTree
	T8b$LLn<-ave(T8b$Year, T8b$LATLON,FUN=function(x) length(x)) 
T8b<-subset(T8b, LLn>7)
myLMall<- function(T8b) 
	{
	lm<-lm(soilw~Year, data=T8b) #, na.exclude)   
		out<-c(lm$coefficients[1],
		lm$coefficients[2], 
		length(lm$model$INVYR), 
		summary(lm)$coefficients[2,2], 
		summary(lm)$fstatistic[1], 
		summary(lm)$r.squared,
		pf(summary(lm)$fstatistic[1], summary(lm)$fstatistic[2],
		summary(lm)$fstatistic[3], lower.tail = FALSE))	
	names(out)<-c("SWintercept", "SWslope","SWn", "SWse", "SW.F.stat", "SW.r.squared", "SW.p.value")
	return(out)
	}
 
TPAfunSW<-ddply(T8b, "LATLON", myLMall)  
S4sw1<-merge(SW3, TPAfunSW, all.x=TRUE, by="LATLON")

#open file if not run above
#S3<-read.csv("/Volumes/m-z/tda210/USFS/FullSAPtrees.2.28.2014.csv")
TCP1<-S3
TCP1$iLAT<-((ceiling(TCP1$LAT*a))/a)-0.5 #+0.25  #round down to integer then add .25 ; switched /* removed ).25 addition
	TCP1$iLON<-((floor(TCP1$LON*b))/b)+0.5 #-0.25  #round down to integer then subtract .25; switched /*  ).25 addition
	TCP1$LATLONYR<-(TCP1$iLAT*10000000000000) +(TCP1$iLON*(-100000000)) +TCP1$INVYR
	TCP1$LATLON<-as.factor(TCP1$iLAT*100000000+ TCP1$iLON*-1000) #LATLON identifyier

###merge Soil water data with TREE groups
S4sw<-merge(TCP1, S4sw1,all.x=TRUE, by="LATLONYR") 
#S4sw<-merge(S4RM, S4sw1,all.x=TRUE, by="LATLONYR") 

S4sw$SWdrought<-ifelse(S4sw$INVYR>S4sw$droughtYRLL2 & (S4sw$INVYR-(S4sw$REMPER))<= S4sw$droughtYRLL2, "YES.drought", "NO.drought")
S4sw$SWpluvial<-ifelse(S4sw$INVYR>S4sw$pluvialYRLL2 & (S4sw$INVYR-(S4sw$REMPER))<= S4sw$pluvialYRLL2, "YES.pluvial", "NO.pluvial")

print("done SW")


#####PRECIP DATA from Prec0.5#######################################################
  #load prec data (http://www.esrl.noaa.gov/psd/data/gridded/data.precl.html) - precip.mon.mean.0.5x0.5.nc
nc<-open.ncdf("/Users/robertbooth/Documents/TDA210/USFSFIA/precip.mon.mean.0.5x0.5.nc")
	# get latitudes and longitudes for NORTH AMERICA
	lat <- get.var.ncdf(nc,"lat",verbose=F)
	lat<-lat[70:129]
	nlat <- dim(lat)
	lon <- get.var.ncdf(nc,"lon")
	lon<-lon[360:719]
	nlon <- dim(lon)
	# get the data and attributes North AMERICA
	precip.array <- get.var.ncdf(nc,"precip",start= c(360,70,1), count= c(360,60,768)) 
	close.ncdf(nc)
	precip.vec.long <- as.vector(precip.array)
	precip.mat <- matrix(precip.vec.long, nrow=nlon*nlat*12, ncol=64, byrow=FALSE) #nrow=nlon*nlat*65.25
	lonlatyr <- expand.grid(lon,lat,1:12) #1948:2013),
	precip.df02 <- data.table(cbind(lonlatyr,precip.mat))
	setnames(precip.df02,c("lon","lat","Month", 1948:2011)) 

Psub<-subset(precip.df02, lon>230 & lon<294) # Month>4 & Month<10 #lat>10 & lat<65 & lon>230 & lon<294
Psub$twelve<-Psub[[67]]
setnames(Psub, "twelve", "2012")
Psub.ml1<-reshape(Psub, direction="long", varying=list(names(Psub)[4:68]), v.names="Precip",idvar=c("lon","lat","Month"), timevar="Year",times=1948:2012) 
Psub.ml<-subset(Psub.ml1, Year>=1985)

Psub.ml$iLAT<-Psub.ml$lat-0.25
Psub.ml$iLON<-Psub.ml$lon+0.25-360
Psub.ml$LATLON<-(Psub.ml$iLAT*100000000)+(Psub.ml$iLON*(-1000))
Psub.ml$LATLONYR<-floor((Psub.ml$iLAT*10000000000000)+(Psub.ml$iLON*(-100000000))+(Psub.ml$Year) + (Psub.ml$Month/10)) #water year

Psub.ml$Prec4<-ave(Psub.ml$Precip, Psub.ml$LATLONYR, FUN=function(x) mean(x,na.rm=TRUE))   #changed sum to mean to account for more grids in 2x2 setup
Psub.ml$Prec<-(Psub.ml$Prec4) # 12 months in year when using mean -> not neccesary when using mean stream discharge

Psub.ml$Zprec<-ave(Psub.ml$Prec, Psub.ml$LATLON, FUN=function(x) scale(x))
Psub.ml$Zprec[is.na(Psub.ml$Zprec)]<-0

###new stufff
Psub.ml$PRECwgrow<-ifelse(Psub.ml$Month>5 & Psub.ml$Month<9, Psub.ml$Precip, NA)

Psub.ml$PREC1<-ave(Psub.ml$Precip, Psub.ml$LATLONYR, FUN=function(x) mean(x, na.rm=TRUE))  
Psub.ml$PREC1grow<-ave(Psub.ml$PRECwgrow, Psub.ml$LATLONYR, FUN=function(x) mean(x, na.rm=TRUE))  

PRECsub<-subset(Psub.ml, !duplicated(LATLONYR), Year>=1991, select=c("PREC1", "PREC1grow", "LATLONYR", "Year", "iLAT", "iLON"))
PRECfinal<-PRECsub
PRECfinal$PREC5sum<-PRECfinal$PREC1
PRECfinal$PREC5growsum<-PRECfinal$PREC1grow

for(i in 1:4)
	{
	PRECsub2<-subset(Psub.ml, !duplicated(LATLONYR), select=c("PREC1", "PREC1grow", "LATLONYR"))
	PRECsub2$LATLONYR<-PRECsub2$LATLONYR+i
	setnames(PRECsub2, c("PREC1","PREC1grow"), c("PREC2","PREC2grow"))
	PRECfinal<-merge(PRECfinal, PRECsub2, by="LATLONYR")
	PRECfinal$PREC5sum<-PRECfinal$PREC5sum+PRECfinal$PREC2
	PRECfinal$PREC5growsum<-PRECfinal$PREC5growsum+PRECfinal$PREC2grow
	PRECfinal<-subset(PRECfinal, select=-c(PREC2,PREC2grow))
	}

PRECfinal$PREC5grow<-PRECfinal$PREC5growsum/5
PRECfinal$PREC5<-PRECfinal$PREC5sum/5

Psub.final<-subset(Psub.ml, !duplicated(LATLONYR) & Year>=1991, select=c("Prec", "Zprec", "LATLONYR", "LATLON"))
Psub.final<-merge(Psub.final, PRECfinal, by="LATLONYR")


#################### PRECIP TRENDS and extremes
T8b<-data.table(NULL)
T8b<-subset(Psub.final, Prec>0 , select=c( "Year", "Zprec", "LATLON")) #TPAsumSapTree
	T8b$LLn<-ave(T8b$Year, T8b$LATLON,FUN=function(x) length(x)) 
T8b<-subset(T8b, LLn>7)
myLMall<- function(T8b) 
	{
	lm<-lm(Zprec~Year, data=T8b) #, na.exclude)   
		out<-c(lm$coefficients[1],
		lm$coefficients[2], 
		length(lm$model$Year), 
		summary(lm)$coefficients[2,2], 
		summary(lm)$fstatistic[1], 
		summary(lm)$r.squared,
		pf(summary(lm)$fstatistic[1], summary(lm)$fstatistic[2],
		summary(lm)$fstatistic[3], lower.tail = FALSE))	
	names(out)<-c("PRECintercept", "PRECslope","PRECn", "PRECse", "PREC.F.stat", "PREC.r.squared", "PREC.p.value")
	return(out)
	}
 
TPAfunPREC<-ddply(T8b, "LATLON", myLMall)  
PREC<-merge(Psub.final, TPAfunPREC, by="LATLON")
S4swPREC<-merge(S4sw, PREC, all.x=TRUE, by="LATLONYR")

S4swPREC$PRECdrought<-ifelse(S4swPREC$INVYR>S4swPREC$droughtYRLL2 & (S4swPREC$INVYR-(S4swPREC$REMPER))<= S4swPREC$droughtYRLL2, "YES.PRECdrought", "NO.PRECdrought")
S4swPREC$PRECpluvial<-ifelse(S4swPREC$INVYR>S4swPREC$pluvialYRLL2 & (S4swPREC$INVYR-(S4swPREC$REMPER))<= S4swPREC$pluvialYRLL2, "YES.PRECpluvial", "NO.PRECpluvial")

print("done PREC")

#redo so no overlap
S4swPREC$pluvial.trend<-ifelse(S4swPREC$PRECslope>0 & S4swPREC$PREC.p.value<0.5,  "pluvial.trend", NA) #& S4swPREC$PREC.p.value<0.2, S4swPREC$PRECslope>0.005
S4swPREC$drought.trend<-ifelse(S4swPREC$PRECslope<(0) & S4swPREC$PREC.p.value<0.5 ,  "drought.trend", NA) #(S4swPREC$PRECslope<(-0.02) 
S4swPREC$SW.pluvial[S4swPREC$SWpluvial=="YES.pluvial"]<-"SW.pluvial"
S4swPREC$SW.drought[S4swPREC$SWdrought=="YES.drought"]<-"SW.drought"
S4swPREC$No.trend<-ifelse(is.na(S4swPREC$drought.trend) & is.na(S4swPREC$pluvial.trend), "No.trend",NA)
S4swPREC$No.extreme<-ifelse(is.na(S4swPREC$SW.pluvial) & is.na(S4swPREC$SW.drought),"No.extreme", NA)
S4swPREC$Stable.climate<-ifelse(S4swPREC$No.trend=="No.trend" & S4swPREC$No.extreme=="No.extreme", "Stable.climate",NA)


##############################################################################################################
####### TEMPERATURE http://www.esrl.noaa.gov/psd/data/gridded/data.ghcncams.html
##  - lat lon centered 0.5 grid. (i.e. bottom right grid corner -.25N and +.25E)
nc<-open.ncdf("/Users/robertbooth/Documents/TDA210/USFSFIA/air.mon.mean.nc")
	# get latitudes and longitudes for NORTH AMERICA
	lat <- get.var.ncdf(nc,"lat",verbose=F)
	lat<-lat[70:129]
	nlat <- dim(lat)
	lon <- get.var.ncdf(nc,"lon")
	lon<-lon[360:719]
	nlon <- dim(lon)
	# get the data and attributes North AMERICA
	temp.array <- get.var.ncdf(nc,"air",start= c(360,70,1), count= c(360,60,780)) 
	close.ncdf(nc)
	tempw.vec.long <- as.vector(temp.array)
	tempw.mat <- matrix(tempw.vec.long, nrow=nlon*nlat*12, ncol=65, byrow=FALSE) #nrow=nlon*nlat*65.25
	lonlatyr <- expand.grid(lon,lat,1:12)
	tempw.df02 <- data.table(cbind(lonlatyr,tempw.mat))
	setnames(tempw.df02,c("lon","lat","Month", 1948:2012)) 

TEMPsub1<-subset(tempw.df02, lon>230 & lon<294)
TEMPsub.ml1<-reshape(TEMPsub1, direction="long", varying=list(names(TEMPsub1)[4:68]), v.names="tempw",idvar=c("lon","lat","Month"), timevar="Year",times=1948:2012) 
TEMPsub.ml<-subset(TEMPsub.ml1, Year>=1985)

TEMPsub.ml$iLAT<-TEMPsub.ml$lat-0.25
TEMPsub.ml$iLON<-(TEMPsub.ml$lon+0.25)-360 
TEMPsub.ml$LATLONYR<-floor((TEMPsub.ml$iLAT* 10000000000000)+(TEMPsub.ml$iLON*(-100000000)) + (TEMPsub.ml$Year)) 
#TEMPsub.ml$LATLONYRgrow<-ifelse(TEMPsub.ml$Month>4 & TEMPsub.ml$Month<10, TEMPsub.ml$LATLONYR, NA)

TEMPsub.ml$tempwgrow<-ifelse(TEMPsub.ml$Month>5 & TEMPsub.ml$Month<9, TEMPsub.ml$tempw, NA)

TEMPsub.ml$temp1<-ave(TEMPsub.ml$tempw, TEMPsub.ml$LATLONYR, FUN=function(x) mean(x, na.rm=TRUE))  
TEMPsub.ml$temp1grow<-ave(TEMPsub.ml$tempwgrow, TEMPsub.ml$LATLONYR, FUN=function(x) mean(x, na.rm=TRUE))  

TEMPsub<-subset(TEMPsub.ml, !duplicated(LATLONYR), select=c("temp1", "temp1grow", "LATLONYR", "Year", "iLAT", "iLON"))
TEMPfinal<-TEMPsub
TEMPfinal$temp5sum<-TEMPfinal$temp1
TEMPfinal$temp5growsum<-TEMPfinal$temp1grow

for(i in 1:4)
	{
	TEMPsub2<-subset(TEMPsub.ml, !duplicated(LATLONYR), select=c("temp1", "temp1grow", "LATLONYR"))
	TEMPsub2$LATLONYR<-TEMPsub2$LATLONYR+i
	setnames(TEMPsub2, c("temp1","temp1grow"), c("temp2","temp2grow"))
	TEMPfinal<-merge(TEMPfinal, TEMPsub2, by="LATLONYR")
	TEMPfinal$temp5sum<-TEMPfinal$temp5sum+TEMPfinal$temp2
	TEMPfinal$temp5growsum<-TEMPfinal$temp5growsum+TEMPfinal$temp2grow
	TEMPfinal<-subset(TEMPfinal, select=-c(temp2,temp2grow))
	}

TEMPfinal$temp5grow<-TEMPfinal$temp5growsum/5
TEMPfinal$temp5<-TEMPfinal$temp5sum/5

TEMP3<-subset(TEMPfinal, select=c(LATLONYR, temp5grow, temp5))

S4swPRECtemp<-merge(S4swPREC, TEMP3, all.x=TRUE, by="LATLONYR")

print("done TEMP")


RM5merged<-merge(RM5, S4swPRECtemp, all.x=TRUE, by=c("PLT_CN", "INVYR", "REMPER"))
RM5merged<-subset(RM5merged, PREVSTOCKING5mid>0 & DIAbeginmean>0)
write.csv(RM5merged, file = "/Volumes/m-z/tda210/USFS/RM5mergedfullcut1991.11.17.2014.csv")
print("Done and Done")




