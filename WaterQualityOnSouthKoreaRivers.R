rm(list = ls())

#======================================================================#
### Packages used
library(sf)
library(readr)
library(stringr)
library(MASS)
library(ggplot2)
library(spacetime)
library(tidyr)
library(dplyr)
library(purrr)
library(funHDDC)
library(lubridate)
library(nnet) 
library(pROC) 
library(caret) 
library(modelr)
library(ggh4x)
library(mice)



#======================================================================#
### Load Data
# region data
site_adjusted <- st_read("RSK_River_reproject_Region/RSK_River_reproject/site_adjusted.shp")

# station coordinate data
stations_ori <- st_read("Sites Point Shapefile/site_07_09.shp") 
stations_crs <- st_crs(stations_ori)

# river information data
edges <- st_read("SouthKoreaRivers/SouthKoreaRivers/RiversSouthKorea.shp")
edges <- st_transform(edges,stations_crs)

# landscape data 
wsh_pred_response <- st_read("wsh_pred_response - Copy/wsh_pred_response - Copy.shp")



#======================================================================#
### Data Wrangling
# check if we have NAs in the region column
sum(is.na(site_adjusted$sites_Rive)) 
summary(site_adjusted$sites_Rive) # we have two NA's for the region

# rename the site_adjusted table to stat_region
stat_region <- site_adjusted 

# there are two NAs, so fill the NAs by copying the previous value
stat_region <- na.locf(na.locf(stat_region),fromLast=TRUE)
#summary(stat_region$sites_Rive)

# plot the graph and founded that we have an outlier with region 3
stat_region %>% filter(sites_Rive == 3) %>% ggplot()+
  geom_sf(aes(color = sites_Rive, geometry = geometry)) #, geometry = c("station__1", "station__2"))) +
  viridis::scale_color_viridis(name = "region")

# guess the outlier's site using the planar coordinates 
region_3 <- stat_region %>% filter(sites_Rive == 3)
stat_region$sites_Rive[stat_region$station_co == "3101A70"] = 1 # change the assumed outlier (stationID 3101A70) value

# convert foreign object to an sf object 
stat_region <- st_as_sf(x = stat_region, 
                        coords = c("station__1", "station__2"),
                        crs = 4326) 



#======================================================================#
# visualize the stations color by region
ggplot(edges) +
  geom_sf(alpha = .1) +
  geom_sf(data = stat_region, size = 1, aes(color = sites_Rive), pch = 20) +
  coord_sf(crs = st_crs(stat_region)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  viridis::scale_color_viridis(name = "region") 



#======================================================================#
### Create Time Series List & Append the Region Columns for Each Station
# list all the files in the path
file.lists <- list.files(path ="Raw station files/flow_wqStations", pattern = ".csv", all.files = TRUE, full.names = TRUE)

# create an empty list named ts.list (time series list)
ts.list = list()

# use for loop to load the .csv files
for(i in 1:length(file.lists)) { 
  ts.list[[i]] <- read.csv(file.lists[[i]])
  
  # delete the second and fourth columns (we can't read these two columns)
  ts.list[[i]][ , c(2, 4)] <- list(NULL)
  
  # specify column names (rename)
  colnames(ts.list[[i]]) <- c('id', 'year', 'measurement_date', 'SS', 'TP', 'flow', 'station_id')
  
  # need to change the measurement_date type from character to time using lubridate
  ts.list[[i]]$measurement_date <- ymd(ts.list[[i]]$measurement_date)
  
  # remove the date part (only keep year-month)
  ts.list[[i]]$measurement_date <- format(ts.list[[i]]$measurement_date, "%Y-%m")
  
  # add the region column to each dataframe
  region_value = stat_region$sites_Rive[which(stat_region$Regression == unique(ts.list[[i]]$station_id))] # unique function == DISTINCT in sql (remove duplicates and show the unique values)
  ts.list[[i]]$region = rep(region_value, nrow(ts.list[[i]])) # repeat it with the number of rows (nrows)
}

# make a list by removing .csv from file names
listnames <- dir("/Raw station files/flow_wqStations/") %>% str_remove(pattern = ".csv") 

# name each dataframe with site names
names(ts.list) <- listnames 

# making ts.list into data.frame
df.ts.allst <- do.call(rbind,ts.list) %>% 
  as_tibble() %>% 
  group_by(station_id,year,measurement_date,region) %>% 
  summarise_at(vars(SS:flow),~median(.,na.rm=T)) %>% 
  mutate(date=lubridate::ymd(paste0(measurement_date,"-01"))) %>% # 
  ungroup()
















#======================================================================#
##### TP #####
#======================================================================#
### Data Imputation for TP
# missing data plots for TP
df.ts.allst %>%
  dplyr::select(station_id,region,date,TP) %>%
  filter(date>=ymd("2014-01-01")) %>% 
  pivot_wider(names_from = "date", values_from = TP) %>% 
  dplyr::select(region,`2014-01-01`:`2019-12-01`) %>%
  naniar::gg_miss_fct(fct = region) 

date_to_lab <- function(dd){
  datenams <- sort(ymd(unique(df.ts.allst$date))) 
  datenams <- datenams[datenams>=ymd("2014-01-01")] 
  paste0("t",which(datenams%in%ymd(dd))) 
}

ImpTP <- df.ts.allst %>% 
  dplyr::select(station_id,region,date,TP) %>% 
  filter(date>=ymd("2014-01-01"),region!=1) %>% 
  pivot_wider(names_from = "date", values_from = TP) %>% 
  rename_at(vars(`2014-01-01`:`2019-12-01`),date_to_lab) %>% 
  mice::mice() %>% mice::complete() 



#======================================================================#
### Functional Analysis for TP
data_TP <- ImpTP %>% 
  dplyr::select(-station_id,-region) %>% 
  mutate_all(~log(.+0.01)) %>% 
  as.matrix() %>% t() 
colnames(data_TP) <- ImpTP$station_id 

time_labels <- tibble(dates=df.ts.allst %>%
                        filter(date>=ymd("2014-01-01")) %>%
                        pull(date) %>% unique() %>% sort(),
                      time_nam=ImpTP %>% 
                        dplyr::select(-station_id,-region) %>% 
                        colnames())

basis<- create.bspline.basis(c(1,72), nbasis=36)                               
varsmooth<-smooth.basis(argvals=1:72,y=data_TP,fdParobj=basis)$fd                

# test for clusters 
set.seed(91213345)
clust_test <- funHDDC(varsmooth,K=3:6,model=c("AkjBQkDk"), 
                      init="kmeans",threshold=0.2)

cluster_label <- clust_test$class
plot(varsmooth,col=cluster_label,
     main=expression(log(TP)))

num_clust <- max(cluster_label) 

# number of stations assigned to each cluster
table(cluster_label) # 67 67 71

#get mean function for each cluster
meanvar_byclust<- lapply(1:num_clust,function(clust){                                 
  mean.fd(fd(varsmooth$coefs[,which(cluster_label==clust)],varsmooth$basis))
})

# plot data and cluster means
par(mfrow=c(2,ceiling(num_clust/2)))
for(k in 1:num_clust){
  matplot(data_TP[,cluster_label==k],
          type="l",col=grey(0.5),lwd=0.5, 
          xlab="month",
          ylab=ifelse(k==1,"tp median (e^y)",""), 
          main=paste("Cluster ",k))
  lines(meanvar_byclust[[k]],type="l",col=k+1,lwd=1.5)
}



#======================================================================#
# Multinomial regression for TP clusters against landscape feat.
datareg_TP <- tibble(station_id=colnames(data_TP),
                     clust=case_when(
                       cluster_label==1~"clust1",
                       cluster_label==2~"clust2",
                       cluster_label==3~"clust3")) %>% 
  left_join(wsh_pred_response %>% 
              as_tibble() %>% 
              dplyr::select(site,area_for:av_slope),
            by=c("station_id"="site")) %>% 
  filter(!is.na(area_for))

tp.formula <- as.formula(paste0("clust ~ ",
                                paste0(names(datareg_TP[,-(1:2)]),
                                       collapse="+"))) 

reg_TP <- multinom(tp.formula, data = datareg_TP) 
summary(reg_TP)
z<- summary(reg_TP)$coefficients/summary(reg_TP)$standard.errors 

# 2-tailed z test
(1 - pnorm(abs(z), 0, 1)) * 2 



#======================================================================#
### Test the model efficiency
# divide the dataset into the train and test data
cl1 <- sample(which(datareg_TP$clust == "clust1"), 0.33*67) # 22
cl2 <- sample(which(datareg_TP$clust == "clust2"), 0.33*67) # 22
cl3 <- sample(which(datareg_TP$clust == "clust3"), 0.34*71) # 24 --- total 68

train = c(cl1, cl2, cl3)
test = datareg_TP[-train,]

# create confusion matrix
reg_TP_probs = predict(reg_TP, test, type = "class") # same as the fitted.values
confusionMatrix(table(reg_TP_probs, test$clust), dnn = c("Predicted", " Actual"))

# plot the ROC curve and get its AUC 
test_roc = multiclass.roc(test$clust, as.numeric(unlist(reg_TP_probs)), 
                          plot = TRUE, print.auc = TRUE, main = "ROC Curve", legacy.axes = TRUE)



#======================================================================#
### Missing Stations Map for TP 
# get the missing data column
missing <- df.ts.allst %>% 
  dplyr::select(station_id,region,date,TP) %>% 
  filter(date>=ymd("2014-01-01")) %>% 
  pivot_wider(names_from = "date", values_from = TP) %>% 
  dplyr::select(station_id,`2014-01-01`:`2019-12-01`) %>% 
  pivot_longer(contains("-"),names_to="date",values_to="TP") %>% 
  group_by(station_id) %>% 
  summarise(TPmiss=mean(is.na(TP)))

missing <- missing %>% filter(missing$TPmiss > 0)
missing_coord <- missing %>% left_join(stat_region, by = c("station_id" = "station_co"))

missing_coord <- st_as_sf(x = missing_coord, 
                          coords = c("station__1", "station__2"), 
                          crs = 4326) 
ggplot(edges) +
  geom_sf(alpha = .1) +
  geom_sf(data = missing_coord, size = 1, aes(color = as.factor(TPmiss)), pch = 20) +
  coord_sf(crs = st_crs(missing_coord)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  viridis::scale_color_viridis(name = "TP missing stations (%)", discrete = T) 



#======================================================================#
### TP Multinomial Regression Visualization
datareg_TP %>% data_grid(area_for = seq_range(area_for, 20), .model = reg_TP) %>% 
  add_predictions(reg_TP,)

grid_allvars <- map_df(names(datareg_TP %>% dplyr::select(area_for:av_slope)),
                       function(.x){
                         out <- data_grid(datareg_TP,
                                          seq_range(!!as.name(.x),20), .model = reg_TP) %>% 
                           dplyr::select(!ends_with(.x))
                         names(out)[names(out)==paste0("seq_range(",.x,", 20)")] = .x
                         out %>% 
                           mutate(predictor=.x)
                       })

lev.nams <- c("area_for","area_urb","area_agr","area_soila",
              "area_soilb","area_soilc","av_elev","av_slope")
lab.nams <- c("forest","urban","agro","soil(a)",
              "soil(b)","soil(c)","elev","slope")

plotprob_data <- data.frame(grid_allvars,
                            grid_allvars %>% predict(reg_TP, ., type = "probs")) %>% 
  pivot_longer(clust1:clust3,names_to="cluster",values_to = "fitted_probs") %>% 
  mutate(thepredvalue=case_when(
    predictor=="area_for" ~ area_for,
    predictor=="area_urb" ~ area_urb,
    predictor=="area_agr" ~ area_agr,
    predictor=="area_soila" ~ area_soila,
    predictor=="area_soilb" ~ area_soilb,
    predictor=="area_soilc" ~ area_soilc,
    predictor=="av_elev" ~ av_elev,
    predictor=="av_slope" ~ av_slope)) %>% 
  mutate(predictor=factor(predictor,
                          levels=lev.nams,
                          labels=lab.nams)) %>% 
  dplyr::select(cluster,predictor,thepredvalue,fitted_probs)

plotprob_data %>% 
  ggplot(.,aes(y=fitted_probs,x=thepredvalue)) +
  geom_line(aes(col=cluster))+
  facet_grid2(predictor~cluster,scales="free",independent="all") + 
  theme_minimal(base_size=8) +
  theme(legend.position = "none") +
  xlab("predictor values") +
  ylab("fitted probablities")
















#======================================================================#
##### SS #####
#======================================================================#
# Data Imputation for SS
#missing data plots for SS
df.ts.allst %>%
  dplyr::select(station_id,region,date,SS) %>%
  filter(date>=ymd("2014-01-01")) %>% 
  pivot_wider(names_from = "date", values_from = SS) %>%
  dplyr::select(region,`2014-01-01`:`2019-12-01`) %>%
  naniar::gg_miss_fct(fct = region) 

ImpSS <- df.ts.allst %>% 
  dplyr::select(station_id,region,date,SS) %>% 
  filter(date>=ymd("2014-01-01"),region!=1) %>% 
  pivot_wider(names_from = "date", values_from = SS) %>% 
  rename_at(vars(`2014-01-01`:`2019-12-01`),date_to_lab) %>% 
  mice::mice() %>% mice::complete() 



#======================================================================#
### Functional Analysis for SS
# make SS matrix with months in the rows and stations in columns
data_SS <- ImpSS %>% 
  dplyr::select(-station_id,-region) %>% 
  mutate_all(~log(.+0.01)) %>% 
  as.matrix() %>% t()
colnames(data_SS) <- ImpSS$station_id

basisSS<- create.bspline.basis(c(1,72), nbasis=36)
varsmoothSS<-smooth.basis(argvals=1:72,y=data_SS,fdParobj=basisSS)$fd                

set.seed(91213345)
clust_testSS <- funHDDC(varsmoothSS,K=2:7,model=c("AkjBQkDk"),
                        init="kmeans",threshold=0.2)

cluster_labelSS <- clust_testSS$class
plot(varsmoothSS,col=cluster_labelSS,
     main=expression(log(SS)))

num_clustSS <- max(cluster_labelSS) 

# number of stations assigned to each cluster
table(cluster_labelSS) 

#get mean function for each cluster
meanvar_byclustSS<- lapply(1:num_clustSS,function(clust_ss){
  mean.fd(fd(varsmoothSS$coefs[,which(cluster_labelSS==clust_ss)], varsmoothSS$basis))
})

# plot data and cluster means
par(mfrow=c(1,num_clustSS))
for(k in 1:num_clustSS){
  matplot(data_SS[,cluster_labelSS==k],
          type="l",col=grey(0.5),lwd=0.5, 
          xlab="month",
          ylab=ifelse(k==1,"ss median",""), 
          main=paste("Cluster ",k))
  lines(meanvar_byclustSS[[k]],type="l",col=k+1,lwd=1.5)
}



#======================================================================#
### Multinomial regression for SS clusters against landscape feat.
datareg_SS <- tibble(station_id=colnames(data_SS),
                     clust=case_when(
                       cluster_labelSS==1~"clust1",
                       cluster_labelSS==2~"clust2")) %>%
  left_join(wsh_pred_response %>% 
              as_tibble() %>% 
              dplyr::select(site,area_for:av_slope),
            by=c("station_id"="site")) %>% 
  filter(!is.na(area_for))

ss.formula <- as.formula(paste0("clust ~ ",
                                paste0(names(datareg_SS[,-(1:2)]),
                                       collapse="+"))) 

reg_SS <- multinom(ss.formula, data = datareg_SS) 
summary(reg_SS)
z_ss <- summary(reg_SS)$coefficients/summary(reg_SS)$standard.errors 

# 2-tailed z test
(1 - pnorm(abs(z_ss), 0, 1)) * 2 

pp_ss <- prop.table(table(datareg_SS$clust)) 



#======================================================================#
# divide to train and test data
cl1_ss <- sample(which(datareg_TP$clust == "clust1"), 0.53*110) # 58.3 -> 58
cl2_ss <- sample(which(datareg_TP$clust == "clust2"), 0.47*96) # 45.12 -> 45

train_SS = c(cl1_ss, cl2_ss) 
test_SS = datareg_SS[-train_SS,]

# confusion matrix
reg_SS_probs = predict(reg_SS, test_SS, type = "class") 
confusionMatrix(table(reg_SS_probs, test_SS$clust,dnn = c("Predicted", " Actual")))

# plot the ROC curve and get its AUC
test_roc_ss = roc(test_SS$clust, as.numeric(unlist(reg_SS_probs)), 
                  plot = TRUE, print.auc = TRUE, main = "ROC Curve", legacy.axes = TRUE)



#======================================================================#
### Missing stations map for SS
# get the missing data column
missing_ss <- df.ts.allst %>% 
  dplyr::select(station_id,region,date,SS) %>% 
  filter(date>=ymd("2014-01-01")) %>% 
  pivot_wider(names_from = "date", values_from = SS) %>% 
  dplyr::select(station_id,`2014-01-01`:`2019-12-01`) %>% 
  pivot_longer(contains("-"),names_to="date",values_to="SS") %>% 
  group_by(station_id) %>% 
  summarise(SSmiss=mean(is.na(SS)))

# only filter the missing value stations
missing_ss <- missing_ss %>% filter(missing_ss$SSmiss > 0)

missing_coord_ss <- missing_ss %>% left_join(stat_region, by = c("station_id" = "station_co"))

missing_coord_ss <- st_as_sf(x = missing_coord_ss, 
                             coords = c("station__1", "station__2"), crs = 4326) 

ggplot(edges) +
  geom_sf(alpha = .1) +
  geom_sf(data = missing_coord_ss, size = 1, aes(color = as.factor(SSmiss)), pch = 20) +
  coord_sf(crs = st_crs(missing_coord_ss)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  viridis::scale_color_viridis(name = "SS missing stations (%)", discrete = T) 



#======================================================================#
### SS Multinomial Regression Visualization
datareg_SS %>% data_grid(area_for = seq_range(area_for, 20), .model = reg_SS) %>% 
  add_predictions(reg_SS,)
data_grid(datareg_SS, .model = reg_SS)

grid_allvars_SS <- map_df(names(datareg_SS %>% dplyr::select(area_for:av_slope)),
                          function(.x){
                            out <- data_grid(datareg_SS,
                                             seq_range(!!as.name(.x),20), .model = reg_SS) %>% 
                              dplyr::select(!ends_with(.x))
                            names(out)[names(out)==paste0("seq_range(",.x,", 20)")] = .x
                            out %>% 
                              mutate(predictor=.x)
                          })

plotprob_data_SS <- data.frame(grid_allvars_SS,
                               grid_allvars_SS %>% predict(reg_SS, ., type = "probs")) %>% 
  rename(clust1 = grid_allvars_SS.....predict.reg_SS.....type....probs..) %>% 
  mutate(clust2 = (1 - clust1)) %>% 
  pivot_longer(clust1:clust2,names_to="cluster",values_to = "fitted_probs") %>% 
  mutate(thepredvalue=case_when(
    predictor=="area_for" ~ area_for,
    predictor=="area_urb" ~ area_urb,
    predictor=="area_agr" ~ area_agr,
    predictor=="area_soila" ~ area_soila,
    predictor=="area_soilb" ~ area_soilb,
    predictor=="area_soilc" ~ area_soilc,
    predictor=="av_elev" ~ av_elev,
    predictor=="av_slope" ~ av_slope)) %>% 
  mutate(predictor=factor(predictor,
                          levels=lev.nams,
                          labels=lab.nams)) %>% 
  dplyr::select(cluster,predictor,thepredvalue,fitted_probs)

plotprob_data_SS %>% 
  ggplot(.,aes(y=fitted_probs,x=thepredvalue)) +
  geom_line(aes(col=cluster))+
  facet_grid2(predictor~cluster,scales="free",independent="all") + 
  theme_minimal(base_size=8) +
  theme(legend.position = "none") +
  xlab("predictor values") +
  ylab("fitted probablities")