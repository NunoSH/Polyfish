
#This script aims to explain the coding procedures of the major steps from the methodology presented at the research paper entitled:
#"An approach to map and quantify the fishing effort of polyvalent passive gear fishing fleets using geospatial data" 
# Published on ICES Journal of Marine Sciences, in 2023. 

# Script written and developed by Nuno Sales Henriques
# with contribution of: Tommaso Russo, Antonio Parisi and Ramiro Magno 

# ----- // ----- #

# At this stage of the analysis, the AIS data has gone through the data processing described on section 3.2 of the  aforementioned paper:
# 3.2.1 - Data cleaning
# 3.2.2 - Split data into fishing trips (trackid)
# 3.2.3 - Remove inadequate fishing trips
# 3.2.4 - Interpolation (1 datapoint per minute)

#AIS dataset format:
#AISdata = AISdata[,c("trackid", "mmsi", "LAT", "LON", "DATE", "SPE", "HEA")]

#Variables description:
  # Trackid -> ID of each fishing trips
  # mmsi -> Vessel identifier
  # LAT -> Latitude
  # LON -> Longitude
  # DATE -> AIS datapoint timestamp (in Chronological Object format)
  # SPE -> Speed over ground (km/h)
  # HEA -> Course over ground (degrees)
#####

#######################################################
#### Section 3.3 - Set up Classification Variables ####
#######################################################

####    Past and Future OVERLAP Variables   ####  

#USING FISHY package
remotes::install_github("ramiromagno/fishy")
library(fishy)

#--------------------------# 
###     PAST OVERLAP     ###
#--------------------------#

start_past_time_range= -0.06 # 90 minutes past AIS ping i, in Chronological Object format
end_past_time_range= -20.16 # 20 days past AIS ping i, in Chronological Object format

START=Sys.time()
overlap_past<-pfind_overlap(AISdata, time_range = c(end_past_time_range, start_past_time_range),
                            speed_range = c(7,Inf), time_col = "DATE",
                            speed_col = "SPE",
                            lon_col = "LON",
                            lat_col = "LAT",
                            by="mmsi")
END=Sys.time()

#--------------------------# 
###    FUTURE OVERLAP    ###
#--------------------------#

start_future_time_range= 0.06 # 90 minutes ahead AIS ping i, in Chronological Object format
end_future_time_range= 20.16 # 20 days ahead AIS ping i, in Chronological Object format

START=Sys.time()
overlap_future<-pfind_overlap(AISdata, time_range = c(start_future_time_range, end_future_time_range),
                              speed_range = c(0,7), time_col = "DATE",
                              speed_col = "SPE",
                              lon_col = "LON",
                              lat_col = "LAT",
                              by="mmsi")
END=Sys.time()


###   Add binary variables PO ->PAST_OVERLAP and FO -> FUTURE_OVERLAP ###

#PAST OVERLAP (PO)
buf_dist_past=0.235 #buffer distance from the hauling datapoints to the deployment datapoints 
AISdata=AISdata %>% add_column(past_overlap= NA)
AISdata$past_overlap = ifelse(AISdata$min_past_d <= buf_dist_past,1,0)
#assign 0 to past_overlap to datapoints with speed that do not correspond to hauls
AISdata$past_overlap[AISdata$SPE > 7]<-0 


#FUTURE OVERLAP (FO)
buf_dist_fut= 0.15 
AISdata=AISdata %>% add_column(future_overlap= NA)
AISdata$future_overlap = ifelse(AISdata$min_fut_di <= buf_dist_fut,1,0)
#assign 0 to overlap_future to datapoints with speed that do not correspond to deployments
AISdata$future_overlap[AISdata$SPE <= 7]<-0
#####


###########################################
#### Section 3.4 - Data classification ####
###########################################

#### 3.4.1 - Data classification with HMM  ####
install.packages("momentuHMM")
library(momentuHMM)


#PrepData to be carried in MomentuHMM by trackid -> ID
names(AISdata)[names(AISdata)== "trackid"]<-"ID" 

#need to order data by ID and date
AISdata <- AISdata[with(AISdata, order(ID, DATE)),]


#use the "PrepData" function to calculate the step length and turning angle between AIS datapoints 

AISdata_prep=prepData(AISdata, type = "LL", coordNames = c("LON", "LAT"))


#HMM parameters

#set the probability distribution for each of the covariates
dist= list(step="gamma",  past_overlap="bern", future_overlap="bern")

#the order of states' parameters is: Drift - Nav - Haul - Deploy -  Slow_nav  


#Set the priors
Par0= list(step= c(0.004, 0.224, 0.045, 0.233, 0.054, 0.004, 0.055, 0.439, 0.037, 0.037),
           overlap_past=c(0.28, 0.2, 0.99, 0.34, 0.08),
           overlap_future=c(0.01, 0.02, 0.18, 0.99, 0.01))

#Fix Distribution parameters
fixPar0= list(step= c(0.002, 0.224, 0.045, 0.233, 0.054, 0.004, 0.055, 0.439, 0.037, 0.037),
              overlap_past=c(0.84, 0.001, 0.95, 0.001, 0.015),
              overlap_future=c(0.01, 0.001, 0.001, 0.97, 0.001))


#fit the model
t_start=Sys.time()
m=fitHMM(data = AISdata_prep, nbStates = 5, 
         Par0 = Par0, dist = dist, estAngleMean = list(angle=T),
         fixPar= fixPar0,
         formula = ~state4(future_overlap) + state3(past_overlap),
         stateNames = c( "Slow_nav", "Nav", "Haul", "Deploy", "Slow_nav"))
t_end=Sys.time()
print(m)
plot(m)
hist(m$data$angle)

states_pred=viterbi(m)
#from this point on, it is irrelevant for us if a vessel is drifting or slowly navigating, 
#therefore we combined the "Drif" and "Slow_nav" into the same group/state:

replace(states_pred, states_pred== 5, 1)

# 
AISdata_predicted=data.frame(AISdata_prep, states_pred)
write.csv(test_set_predicted, "_predicted_8.csv")
#####

#### 3.4.2 - Clean false positives  ####

#   SMOOTH THE TRACKS    #

states_pred_sm= zoo::rollapply(states_pred, 13, function(x) names(which.max(table(x))), partial=TRUE)

AISdata_predicted$sts_smooth=states_pred_sm
###


#    CLEAN FALSE POSITIVES    #


vessel=unique(AISdata_predicted$mmsi)

#clean state 3 (Hauling)
clean_st3=NULL

for(k in 1:length(vessel)){
  vessel_k= AISdata_predicted[which(AISdata_predicted$mmsi == vessel[k]),]
  
  #subset predicted state 3 and state 4 AIS_datapoints
  ais_st3=subset(vessel_k, sts_smooth == "3")
  ais_st4=subset(vessel_k, sts_smooth == "4")
  #calculate distances between state 3 and state 4
  dist_ais = spDists(x = data.matrix(ais_st3[, c("y", "x")]),
                     y = data.matrix(ais_st4[, c("y", "x")]),
                     longlat = T)
  #get the minimum distance (state 3 - state 4)
  min_dist = apply(dist_ais, 1, min)
  ais_st3=data.frame(ais_st3,min_dist)
  #classify state 3 datapoints in state 1 if min_dist > buf_dist_past from state 4 (0.235km)
  ais_st3$sts_smooth= ifelse(ais_st3$min_dist>0.235,1,3)
  ais_st3=ais_st3[,c(1:15)]
  vessel_k=rbind(vessel_k[!vessel_k$sts_smooth == "3",], ais_st3)
  clean_st3=rbind(clean_st3, vessel_k)
}



#Clean state 4 (Deploying)
clean_st3_st4=NULL

for(k in 1:length(vessel)){
  vessel_k= clean_st3[which(clean_st3$mmsi == vessel[k]),]
  
  #subset predicted state 3 and state 4 AIS_datapoints
  ais_st4=subset(vessel_k, sts_smooth == "4")
  ais_st3=subset(vessel_k, sts_smooth == "3")
  #calculate distances between st3 and st4
  dist_ais = spDists(x = data.matrix(ais_st4[, c("y", "x")]),
                     y = data.matrix(ais_st3[, c("y", "x")]),
                     longlat = T)
  #get the minimum distance (state 3 - state 4)
  min_dist = apply(dist_ais, 1, min)
  ais_st4=data.frame(ais_st4,min_dist)
  #classify state 4 datapoints in state 2 if min_dist > buf_dist_past from state 3 (0.15km)
  ais_st4$sts_smooth= ifelse(ais_st4$min_dist>0.15,2,4)
  ais_st4=ais_st4[,c(1:15)]
  vessel_k=rbind(vessel_k[!vessel_k$sts_smooth == "4",], ais_st4)
  clean_st3_st4=rbind(clean_st3_st4, vessel_k)
}

#####

#######################################################
#### Section 3.5 - Footprint and effort assessment ####
#######################################################

#### Match Deployment and Haul datapoints ####

# working functions:

#### #1 SpDist version in C (much faster)  ####

#SpDistN1 in C: is faster than the normal SpDist

spDistsN1 <- function(pt_x, pt_y, pts_x, pts_y, longlat = TRUE) {
  n <- as.integer(length(pts_x))
  dists <- vector(mode = "double", length = n)
  
  # The subsetting at 6 indicates that we want to return the sixth argument
  # passed to `.C()`, i.e. `dists`.
  res <- .C("fishy_dists", pts_x, pts_y, pt_x, pt_y, n, dists, longlat, PACKAGE = "fishy")[[6]]
  
  if (any(!is.finite(res))) {
    nAn <- which(!is.finite(res))
    dx <- abs(pts_x[nAn] - pt_x)
    dy <- abs(pts_y[nAn] - pt_y)
    if (all((c(dx, dy) < .Machine$double.eps^0.5)))
      res[nAn] <- 0
    else stop(paste("non-finite distances in spDistsN1"))
  }
  
  res
}


#### #2 Find_haul function ####
# for each Deployment ping, the function finds the closes(in time) 
# hauling ping within a specified time range and distance 
# Can only be used for one vessel at a time

find_haul <-
  function(df,
           time_range,
           state_deploy,
           state_haul,
           buff_dist,
           time_col = 'DATE',
           state_col = 'sts_smooth',
           lon_col = 'x',
           lat_col = 'y') {
    
    
    dt <- data.table::as.data.table(df)
    dt[['..id..']] <- seq.int(nrow(dt))
    data.table::setkeyv(dt, c(time_col))
    
    deploy = dt[dt[[state_col]] == state_deploy, ]
    n_row = unique(deploy$..id..)
    
    #t_lower <- as.data.frame(dt[ , time_col, with = FALSE])[, 1] + time_range[1]
    #t_upper <- as.data.frame(dt[ , time_col, with = FALSE])[, 1] + time_range[2]
    
    haul_p_id <- rep(NA_integer_, length = length(n_row))
    haul_trip_ID <- rep(NA_integer_, length = length(n_row))
    haul_DATE <- rep(NA_integer_, length = length(n_row))
    haul_x <- rep(NA_integer_, length = length(n_row))
    haul_y <- rep(NA_integer_, length = length(n_row))
    haul_state <- rep(NA_integer_, length = length(n_row))
    depl_haul_distance <- rep(NA_real_, length = length(n_row))
    
    
    # i in deploy$..id..
    for(i in 1:length(n_row)){
      ping_i = deploy[which(deploy$..id.. == n_row[i]),]
      
      dt_subset <-
        dt[dt[[time_col]] >= time_range[1] + ping_i$DATE & # t_lower[i] & #
             dt[[time_col]] < time_range[2] + ping_i$DATE & #  t_upper[i] &
             dt[[state_col]] == state_haul & !(dt[["..id.."]] %in% haul_p_id) ]
      
      
      if (nrow(dt_subset) > 0) {
        
        distances <- spDistsN1(pt_x = ping_i[[lon_col]],
                               pt_y = ping_i[[lat_col]],
                               pts_x = dt_subset[[lon_col]],
                               pts_y = dt_subset[[lat_col]],
                               longlat = TRUE)
        dt_subset$dist = distances
        dt_subset=dt_subset[dt_subset$dist <= buff_dist]
        iii <- dt_subset[which.min(dt_subset$DATE),]
        
        
        if (nrow(iii) == 0) next
        if (length(iii) > 0) {
          haul_p_id[i] <- iii$..id..#[iii]
          haul_trip_ID[i] <- iii$ID
          haul_DATE[i] <- iii$DATE
          haul_x[i] <- iii$x
          haul_y[i] <- iii$y
          haul_state[i] <- iii$sts_smooth
          depl_haul_distance[i] <- iii$dist
          
        }
      }
    }
    
    return(cbind(deploy, haul_p_id,haul_trip_ID, haul_DATE, haul_x, haul_y, haul_state, depl_haul_distance))
    
  }

# Example
# #define time range
start_future_time_range = 0.06 #(90 minutes)
end_future_time_range = 20.16 #(20 days)

example1= find_haul( vess, time_range = c(start_future_time_range, end_future_time_rang),
                 state_deploy = 4, state_haul =3,
                 buff_dist = 0.15,
                 state_col = "sts_smooth",
                 time_col = "DATE",
                 lon_col = "x",
                 lat_col = "y")

#####

#### #3 Pfind_haul function #####

## Find_haul with parallel:
## vessel id (mmsi) is the factor to split 
## the data into chucks to be iterated separately

pfind_haul <- function(df,
                       time_range,
                       state_deploy,
                       state_haul,
                       buff_dist,
                       time_col = 'DATE',
                       state_col = 'sts_smooth',
                       lon_col = 'x',
                       lat_col = 'y',
                       by = 'mmsi') {
  
  #dt <- data.table::as.data.table(df)
  data.table::setkeyv(dt, c(by, time_col))
  
  j_expression <- substitute(
    find_haul(
      df = .SD,
      time_range,
      state_deploy,
      state_haul,
      buff_dist,
      time_col,
      state_col,
      lon_col,
      lat_col
    ),
    env = list(
      time_range = time_range,
      state_deploy = state_deploy,
      state_haul = state_haul,
      buff_dist =buff_dist,
      time_col = time_col,
      state_col = state_col,
      lon_col = lon_col,
      lat_col = lat_col,
      by = by
    )
  )
  
  dt[, j = eval(j_expression), by = by]
  
}
# Example:

# #define time range
start_future_time_range = 0.06 #(90 minutes)
end_future_time_range = 20.16 #(20 days)
 
example2= pfind_haul(df = dt, time_range = c(start_future_time_range, end_future_time_range),
                  state_deploy = 4,
                  state_haul = 3,
                  buff_dist = 0.15,
                  state_col = "sts_smooth",
                  time_col = "DATE",
                  lon_col = "x",
                  lat_col = "y",
                  by= "mmsi")

#####

#### Match Deployment AIS datapoints with Hauling AIS datapoints using function "pfind_haul"####

AISdata = AISdata[,c(1, 4:16)]

#select only hauling and deployment datapoints
fishing = AISdata[AISdata$sts_smooth == 3 | AISdata$sts_smooth == 4,]

#switch to data.table
dt <- data.table::as.data.table(fishing)


# #define time range
start_future_time_range = 0.06 #(90 minutes)
end_future_time_range = 20.16 #(20 days)

strt= Sys.time()
final_eff = pfind_haul(dt, time_range = c(start_future_time_range, end_future_time_range),
                       state_deploy = 4,
                       state_haul = 3,
                       buff_dist = 0.15,
                       state_col = "sts_smooth",
                       time_col = "DATE",
                       lon_col = "x",
                       lat_col = "y",
                       by= "mmsi")
end=Sys.time()

final_eff <- final_eff[complete.cases(final_eff),]


#tidy the data
colnames(final_eff)[c(2,3, 15:22)] = c("Depl_trip_ID", "Depl_DATE", "Depl_p_id", "Haul_p_id",
                                       "Haul_trip_ID", "Haul_DATE", "Haul_x", "Haul_y", "Haul_state", "depl_haul_dist")

#####

#### Identification of complete (deployment + hauling) fishing events ####

#order the data by vessel and DATE
final_eff= final_eff[with(final_eff, order(mmsi, Depl_DATE)),]

#switch data from data.table to data.frame
eff=as.data.frame(final_eff)

vessels <- unique(eff$mmsi)
final_eff_clean = NULL

for(v in 1:length(vessels)){
  vessel_v = eff[which(eff$mmsi == vessels[v]),]
  
  # When difference between consecutive DATE > 0.0105 (15min aprox - chronological object)  => 1
  vessel_v$gap=c(0,diff(vessel_v$Depl_DATE)> 0.0105)
  
  #group by gap: create groups when consecutive datapoints are separated for more than 15 minutes (deploy_group)
  vessel_v$deploy_group=cumsum(vessel_v$gap) + 1
  
  #remove groups with less than 15 datapoints (assumption that deployments are no shorter than 15 minutes)
  vessel_v_clean = vessel_v %>%
    group_by(deploy_group) %>% filter(n() >= 15)
  
  #creat final group (fish_event) with the ID of the deployment trip, the hauling trip, and the group created above (deploy_group).
  vessel_v_clean$fishing_event= paste(vessel_v_clean$Depl_trip_ID, "_", vessel_v_clean$Haul_trip_ID, "_", vessel_v_clean$deploy_group)
  
  #remove groups with less than 15 datapoints. (assumption, again, that no deployment takes less than 15 minutes)
  vessel_v_clean= vessel_v_clean %>%
    group_by(fishing_event) %>% filter(n() >= 15)
  
  final_eff_clean=rbind(final_eff_clean, vessel_v_clean)
  
}
#####

##### Add SOAK TIME #### 
final_eff_clean$soak_time= final_eff_clean$Haul_DATE - final_eff_clean$Depl_DATE

#####

