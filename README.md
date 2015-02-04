# ESM-FIA: Empirical Succession Mapping - R code

Forests are amazing! 

Explore this R code and associated filtered US Forest Service Forest Inventory and Analysis (FIA) data to build the Empirical Succession Mapping and perform analyses. Current analyses:

 Figure 1. Empirical Succession Mapping
 
 Figure 2. 5-year self thinning relationships to disturbance
 
 Figure 3. Influence of climate conditions (temperature & moisture anomalies)
 
 Figure 4. Forward modeling of eastern US forest carbon steady-states

Associated manuscript submitted for review. Feel free to contact me with any questions! 
 
 ~TD Andrews~   


Brief description of "Null2015.csv" variables (FIA nomenclature maintained when possible; "starting" and "PREV" refer to first survey, "ending" refers to resurvey): 

Refer to FIA manual for further details:
http://www.fia.fs.fed.us/library/database-documentation/current/ver6.0/FIADB%20User%20Guide%20P2_6-0-1_final.pdf 

"PLT_CN"  Unique plot identifier
"INVYR" Year of most recent plot census (roughly)             
"REMPER" Census interval             
"startsumTPA" starting trees per acre         
"endsumTPA" ending tree per acre        
"ingrowth2" number of trees grew past 5 in. diameter threshold         
"mortality2" number of trees greater than 5 in. that died during the census interval        
"DIAbeginmean" starting mean diameter        
"DIAendmean" ending mean diameter        
"DIAgrowmean" mean diameter growth       
"cut" indicates trees were harvested    
"mortDIAmean" mean diameter of trees that died during the census interval       
"cutDIAmean" mean diameter of trees that were harvested during the census interval         
"LATLONYR" Latitude, longitude, inventory year combined identifier          
"STOCKING5mid" ending plot stocking sum       
"PREVSTOCKING5mid" starting plot stocking sum  
"DIA5meanalive" ending mean diameter (repeat)      
"PREVDIA5meanalive" starting mean diameter (repeat)         
"BA5meanalive" ending basal area DIA>5in.  
"PREVBA5meanalive" starting basal area DIA>5in.  
"TPAsum5"  ending trees per acre (repeat)             
"PREVTPAsum5" starting trees per acre (repeat)         
"PREVTPAsum10" starting trees per acre DIA>10 (repeat)       
"TPAsum05"  ending sapling stem density - trees per acre  
"PREVTPAsum05" starting sapling stem density - trees per acre        
"sapMortsum" number of saplings that die during the census interval        
"sapRecruitsum" saplings recruited during census interval      
"TPAsum15" ending trees per acre DIA>15          
"PREVTPAsum15" starting trees per acre DIA>15         
"LAT" latitude (fuzzed)               
"LON" longitude (fuzzed)                
"STATECD" US state identifier           
"PREV_PLT_CN"  Previous plot unique identifier     
"STDAGE" estimated stand age            
"PREV_TRE_STDAGE" previous stand age estimate  
"STDORGCD"  planted forest indicator          
"SISP" Dominant species group (coniferous < 300 >= deciduous)              
"PREV_TRE_SISP"  Previous dominant species group (coniferous < 300 >= deciduous; poorly known)   
"carbon5sum" ending above and below ground vegetation carbon DIA>5        
"PREVcarbon5sum" starting above and below ground vegetation carbon DIA>5    
"carbon1sum"  ending above and below ground vegetation carbon DIA>1  
"PREVcarbon1sum"  starting above and below ground vegetation carbon DIA>1     
"deadcarbon5sum"  carbon in trees that die during the census interval DIA>5   
"SPCDdiversity" species richness    
"PREV_SPCDdiversity"previous species richness          
"animal" tree mortality from attributed to animals                   
"insect"  tree mortality from attributed to insects          
"fire" tree mortality from attributed to fire            
"weather"  tree mortality from attributed to weather disturbance           
"disease"  tree mortality from attributed to diseases          
"cutting"  tree mortality from attributed to harvesting          
"unknowndamage"  tree mortality from attributed to unknown causes   
"vegetation" tree mortality from attributed to vegetation        
"DSTRBCD1" large area disturbances           
"PREV_TRE_DSTRBCD1" previous large area disturbances  
"NPP" equal to carbon1sum +   deadcarbon5sum + sapling mortality carbon from previous sapling carbon              
"PREVNPP" equal to PREVcarbon1sum          
"SWdrought5" minimum z-scored soil moisture during the four years prior to INVYR-1       
"SWpluvial5"  maximum z-scored soil moisture during the four years prior to INVYR-1  
"PREC5grow" June, July, August mean precipitation for the five years prior to INVYR           
"PREC5" water-year mean precipitation for the five years prior to INVYR                     
"temp5grow"  June, July, August mean temperature for the five years prior to INVYR          
"temp5" water-year mean temperature for the five years prior to INVYR
