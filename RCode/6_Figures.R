## --------------------------------------------------------------------------------------------------------- ##
##  Project Name: Functional reorganization of North American wintering avifauna                             ##
##                                                                                                           ##
##  Objective:     Here, we examined variation avian fauna in North America over time.                       ##
##                                                                                                           ##
##  Authors:       Juan P. Quimbayo                                                                          ##
##                 Stephen Murphy                                                                            ##
##                 Marta Jarzyna                                                                             ##
##                                                                                                           ##
##  Date:          2023-10-03                                                                                ##
##                                                                                                           ##
##  Notes:         1. This file is intended to provide a guide to the basic                                  ##
##                    project workflow, attempting to 'integrate_analysis' the steps                         ##
##                    necessary to conduct the analyses & visual outputs.                                    ##
##                                                                                                           ##
## --------------------------------------------------------------------------------------------------------- ##

rm (list=ls())

# FIGURES -------------------------
library (dplyr)
library (plyr)
library (maps)
library (mapproj)
library (RColorBrewer)
library (ggplot2)  
library (ggpmisc)
library (ggpubr)
library (cluster)
library (fastcluster)
library (factoextra)
library (reshape2)
library (RColorBrewer)
library (sp)
library (sf)
library (rgdal)
library (raster)
library (rasterVis)
library (viridis)
library (scatterpie)
library (scales)
library (betareg)
library (wesanderson)
library (ggforce)
library (tidyr)
library (ggmcmc)
library (tidybayes)
library (gggap)
library (brms)
library (bayesplot)

# 1. Calling data ---------------------------------------------------------
DB_FIndex <- readRDS("DataInter/DB_FIndex.RDS")
nullFRic <- readRDS("DataInter/nullFRic.rds")
nullFRic$ID_unique <- paste(nullFRic$Abbre,nullFRic$year, sep = "_")
nullFRic <- nullFRic [,c(5,3:4)]

DB_FIndex <- left_join(DB_FIndex, nullFRic)
DB_FIndex$ses_fric <- (DB_FIndex$fric - DB_FIndex$meanNull)/DB_FIndex$SDNull


# 2. Extracting slope values from functional indices ----------------------

 location_years <- split (DB_FIndex, DB_FIndex$Location)

slope_values_sric   <- list ()
slope_values_cfric  <- list ()
slope_values_feve   <- list ()
slope_values_fdiv   <- list ()
slope_values_fori   <- list ()
slope_values_Dim1   <- list ()
slope_values_Dim2   <- list ()

# Loop to extracting estimate value from each linear regression - function was "coefficients"

for (m in 1:length(location_years)){
  slope_values_sric[[m]]  <- summary(glm (sp_richn~year, family = poisson, data = location_years[[m]]))$coefficients[2,1]
  slope_values_cfric[[m]] <- summary(lm (ses_fric~year, data = location_years[[m]]))$coefficients[2,1]
  slope_values_feve[[m]]  <- summary(betareg (feve~year, data = location_years[[m]]))$coefficients$mean[2,1]
  slope_values_fdiv[[m]]  <- summary(betareg (fdiv~year, data = location_years[[m]]))$coefficients$mean[2,1]
  slope_values_fori[[m]]  <- summary(betareg (fori~year, data = location_years[[m]]))$coefficients$mean[2,1]
  slope_values_Dim1[[m]]  <- summary(glm (Dim.1~year, data = location_years[[m]]))$coefficients[2,1]
  slope_values_Dim2[[m]]  <- summary(glm (Dim.2~year, data = location_years[[m]]))$coefficients[2,1]
}

beta_slopes_f_indices <- data.frame(Abbrev=unique(DB_FIndex$Location),
                               coef_sric= do.call("rbind",slope_values_sric),
                               coef_cfric= do.call("rbind",slope_values_cfric),
                               coef_feve= do.call("rbind",slope_values_feve),
                               coef_fdiv= do.call("rbind",slope_values_fdiv),
                               coef_fori= do.call("rbind",slope_values_fori),
                               coef_Dim1= do.call("rbind",slope_values_Dim1),
                               coef_Dim2= do.call("rbind",slope_values_Dim2)) 



# 3. Add coordinates information ------------------------------------------
LocGeoCoor <- readRDS("DataRaw/Coords.rds")
colnames(LocGeoCoor) <- c("Abbrev","Location_Name","Lat","Long")
LocGeoCoor <- data.frame(LocGeoCoor)
Location.Coords <- LocGeoCoor [, c("Abbrev","Lat", "Long")]

# Join databases
beta_slopes_f_indices   <- left_join(beta_slopes_f_indices,LocGeoCoor)




# 4. Figure 1: Functional maps  ----------------------------------------------
# This function was created the functional maps  
map_function <- function (db, f_index, NameIndex, GridCellReso, HighColor, LowColor, MidColor){
  usa <- map_data("usa")
  db <- db
  db <- data.frame (x=db$Long,
                    y=db$Lat,
                    value=f_index)
  
  #Transform in raster data.frame
  r <- rasterFromXYZ(db, crs = CRS("+proj=longlat +datum=WGS84"), res = GridCellReso, digits = 0)
  
  
  # create an extent object for the extent of the United States
  usa_bbox <- c(-125.5, 24.5, -67, 49.5)
  usa_extent <- extent(usa_bbox)
  
  # crop the raster to the extent of the world map
  raster_cropped <- crop(r, extent(usa_extent))
  # convert the cropped raster to a data frame
  df <- data.frame(rasterToPoints(raster_cropped))
  
  
  # Create a raster plot of the variable/index
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = value))+
    geom_polygon(data=usa, aes(x=long, y=lat, group=group),
                 colour="black", fill=NA)+
    labs(fill=NameIndex) +
    scale_fill_gradient2 (low = LowColor, high = HighColor, 
                          mid = MidColor) +
    coord_quickmap() + 
    theme_void()
  
}


## 4.1 Species richness map --------------------------------------------------
Fig_1a <- map_function(db =slopes_f_indices_coord,
                       f_index = beta_slopes_f_indices$coef_sric,
                       NameIndex = expression(Delta~"SppR"),
                       GridCellReso=0.7,
                       HighColor = brewer.pal(9, "YlOrRd"),
                       LowColor = rev(brewer.pal(9, "Blues")),
                       MidColor = "white")

## 4.2 Functional richness map -------------------------------------------------

Fig_1b <- map_function(db =beta_slopes_f_indices,
                         f_index = beta_slopes_f_indices$lm_coef_sesfric,
                         NameIndex = expression(Delta~"SES.FRic"),
                         GridCellReso=1,
                         HighColor = brewer.pal(9, "YlOrRd"),
                         LowColor = rev(brewer.pal(9, "Blues")),
                         MidColor = "white")

## 4.3 Functional evenness map -------------------------------------------------

Fig_1c <- map_function(db =beta_slopes_f_indices,
                       f_index = beta_slopes_f_indices$beta_coef_feve,
                       NameIndex = expression(Delta~"FEve"),
                       GridCellReso=1,
                       HighColor = brewer.pal(9, "YlOrRd"),
                       LowColor = rev(brewer.pal(9, "Blues")),
                       MidColor = "white")

## 4.4 Functional divergence map -------------------------------------------------

Fig_1d <- map_function(db = slopes_f_indices_coord,
                       f_index = beta_slopes_f_indices$beta_coef_fdiv,
                       NameIndex = expression(Delta~"FDiv"),
                       GridCellReso=1,
                       HighColor = brewer.pal(9, "YlOrRd"),
                       LowColor = rev(brewer.pal(9, "Blues")),
                       MidColor = "white")

## 4.5 Functional originality map -------------------------------------------------

Fig_1e <- map_function(db =beta_slopes_f_indices,
                       f_index = beta_slopes_f_indices$beta_coef_fori,
                       NameIndex = expression(Delta~"FOri"),
                       GridCellReso=1,
                       HighColor = brewer.pal(9, "YlOrRd"),
                       LowColor = rev(brewer.pal(9, "Blues")),
                       MidColor = "white")



ggarrange(Fig_1a, Fig_1b, Fig_1c, Fig_1d, Fig_1e,
          align="hv",
          ncol=1, nrow=5,
          font.label = list(size=13, color="black",
                            face="bold", family = "sans"))



# 5. Functional temporal trends  ------------------------------------------

# Call biogeographical information
db_domains <- read.csv("DataInter/RegionalPool_coord.csv", header = T)
colnames(DB_FIndex)[15] <- 'Abbrev' 

db_domains_indices <- left_join(DB_FIndex, db_domains)

# Split dataframe per domain
domains_predict <- split (db_domains_indices, db_domains_indices$DomainName)

results_lm_domains_sric <- list()
results_lm_domains_fric <- list()
results_lm_domains_sesfric <- list()
results_lm_domains_feve <- list()
results_lm_domains_fdiv <- list()
results_lm_domains_fori <- list()

for (i in 1:length(domains_predict)){
    # Species richness
    lm_domain_sric <- lm (sp_richn~year, data = domains_predict[[i]])
    db <- as.data.frame (domains_predict[[i]]) 
    predict_domain_sric <- data.frame (lm_predict=predict(lm_domain_sric,db),
                                     year=db$year,
                                     Domain="Domain")
    results_lm_domains_sric[[i]] <- predict_domain_sric
  
  
    # Functional richness
    lm_domain_fric <- betareg (fric~year, data = domains_predict[[i]])
    db <- as.data.frame (domains_predict[[i]]) 
    predict_domain_fric <- data.frame (lm_predict=predict(lm_domain_fric,db),
                                  year=db$year,
                                  Domain="Domain")
    results_lm_domains_fric[[i]] <- predict_domain_fric
    
    # SES Functional richness
    lm_domain_sesfric <- lm (ses_fric~year, data = domains_predict[[i]])
    db <- as.data.frame (domains_predict[[i]]) 
    predict_domain_sesfric <- data.frame (lm_predict=predict(lm_domain_sesfric,db),
                                       year=db$year,
                                       Domain="Domain")
    results_lm_domains_sesfric[[i]] <- predict_domain_sesfric
    
    # Functional evenness
    lm_domain_feve <- betareg (feve~year, data = domains_predict[[i]])
    db <- as.data.frame (domains_predict[[i]]) 
    predict_domain_feve <- data.frame (lm_predict=predict(lm_domain_feve,db),
                                       year=db$year,
                                       Domain="Domain")
    results_lm_domains_feve[[i]] <- predict_domain_feve
    
    # Functional divergence
    lm_domain_fdiv <- betareg (fdiv~year, data = domains_predict[[i]])
    db <- as.data.frame (domains_predict[[i]]) 
    predict_domain_fdiv <- data.frame (lm_predict=predict(lm_domain_fdiv,db),
                                       year=db$year,
                                       Domain="Domain")
    results_lm_domains_fdiv[[i]] <- predict_domain_fdiv
    
    # Functional originality
    lm_domain_fori <- betareg (fori~year, data = domains_predict[[i]])
    db <- as.data.frame (domains_predict[[i]]) 
    predict_domain_fori <- data.frame (lm_predict=predict(lm_domain_fori,db),
                                       year=db$year,
                                       Domain="Domain")
    results_lm_domains_fori[[i]] <- predict_domain_fori
    
    
  }

# This function was created to shown the temporal trends including the NEON Domains  
linear_figure_function <- function (dbFull, findex, model_equation,dbPre, yNameAxis){
  
  lm_findex_fit <- model_equation
  predicted_df  <- data.frame (findex_pred=predict(lm_findex_fit,dbFull), 
                              year=dbFull$year, Location="Local") 
  
  ggplot (dbFull, aes(x=year, y=findex, shape=Abbrev)) +
    stat_smooth(method = lm, formula = y~x, se = F,fullrange=TRUE, colour="grey90", size=0.5)+
    geom_line (color="black",data =predicted_df, aes(x=year, y=findex_pred, shape=Location))+
    geom_line (color="#800000",data=dbPre[[1]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#DC143C",data=dbPre[[2]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#FF7F50",data=dbPre[[3]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#FFA500",data=dbPre[[4]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#BDB76B",data=dbPre[[5]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#808000",data=dbPre[[6]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#9ACD32",data=dbPre[[7]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#006400",data=dbPre[[8]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#32CD32",data=dbPre[[9]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#20B2AA",data=dbPre[[10]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#AFEEEE",data=dbPre[[11]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#4682B4",data=dbPre[[12]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#00008B",data=dbPre[[13]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#4B0082",data=dbPre[[14]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#8B008B",data=dbPre[[15]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#DA70D6",data=dbPre[[16]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#BC8F8F",data=dbPre[[17]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    labs (x="Years", y=yNameAxis) +
    theme_classic() +
    theme (axis.text.y  = element_text(size=10, angle=0, family = "sans"),
           axis.title.y  = element_text(size=13, angle=90, family = "sans"),
           axis.text.x   = element_text(size=10, angle=0, family = "sans"),
           axis.title.x  = element_text(size=13, angle=0, family = "sans"),
           legend.position = "na")
  
}

# Temporal trends per each functional index
Fig_1.1a <- linear_figure_function (dbFull = db_domains_indices,
                                  findex = db_domains_indices$sp_richn,
                                  model_equation = lm(sp_richn~year, data =db_domains_indices),
                                  dbPre = results_lm_domains_sric,
                                  yNameAxis = expression("SRic"))
Fig_1.1b <- linear_figure_function (dbFull = db_domains_indices,
                                  findex = db_domains_indices$ses_fric,
                                  model_equation = lm(ses_fric~year, data =db_domains_indices),
                                  dbPre = results_lm_domains_sesfric,
                                  yNameAxis = expression("SES.FRic"))
Fig_1.1c <- linear_figure_function (dbFull = db_domains_indices,
                                  findex = db_domains_indices$feve,
                                  model_equation = betareg(feve~year, data =db_domains_indices),
                                  dbPre = results_lm_domains_feve,
                                  yNameAxis = expression("FEve"))
Fig_1.1d <- linear_figure_function (dbFull = db_domains_indices,
                                  findex = db_domains_indices$fdiv,
                                  model_equation = betareg(fdiv~year, data =db_domains_indices),
                                  dbPre = results_lm_domains_fdiv,
                                  yNameAxis = expression("FDiv"))
Fig_1.1e <- linear_figure_function (dbFull = db_domains_indices,
                                  findex = db_domains_indices$fori,
                                  model_equation = betareg(fori~year, data =db_domains_indices),
                                  dbPre = results_lm_domains_fori,
                                  yNameAxis = expression("FOri"))


ggarrange(Fig_1.1a, Fig_1.1b, Fig_1.1c, Fig_1.1d, Fig_1.1e,
          align="hv",
          ncol=1, nrow=5, legend = "none",
          font.label = list(size=13, color="black",
                            face="bold", family = "sans"))


# 6. Figure 2: Cluster and map --------------------------------------------

complete_cluster <- left_join(beta_slopes_f_indices, db2)
cluster_matrix_full <- complete_cluster [,c("Abbrev","lm_coef_sesfric",
                                            "beta_coef_feve","beta_coef_fdiv",
                                            "beta_coef_fori", "slope_Dim.1",
                                            "slope_Dim.2")]

# Estimation how many clusters are identified
fviz_nbclust (cluster_matrix_full[,2:length(cluster_matrix_full)], 
              kmeans,
              method = 'silhouette') 

cluster_k_full  <- kmeans (cluster_matrix_full[,2:length(cluster_matrix_full)], 
                           centers = 5, 
                           nstart = 2000)

cluster_db_full <- complete_cluster
cluster_db_full <- cbind (cluster_db_full, k_cluster_full=cluster_k_full$cluster)

# Boxplot shown clusters 
fig_matrix_full <- cluster_db_full[, c("lm_coef_sesfric",
                                       "beta_coef_feve",
                                       "beta_coef_fdiv",
                                       "beta_coef_fori",
                                       "slope_Dim.1",
                                       "slope_Dim.2",
                                       "k_cluster_full")]

fig_matrix_full <- melt (fig_matrix_full, id=c("k_cluster_full"))
colnames(fig_matrix_full) <- c("k_cluster","findex","slope")
fig_matrix_full$k_cluster <- as.factor (fig_matrix_full$k_cluster)


Fig_2a <- ggplot (droplevels(subset(fig_matrix_full, slope<=0.04)), 
                  aes(x=factor(findex, levels = c("slope_Dim.2",
                                                  "slope_Dim.1",
                                                  "beta_coef_fori",
                                                  "beta_coef_fdiv",
                                                  "beta_coef_feve",
                                                  "lm_coef_sesfric")), 
                      y=slope))+
  geom_boxplot(aes(fill=k_cluster), outlier.color = NA)+
  geom_hline(yintercept = 0, linetype="dashed",color="black")+
  coord_flip()+
  labs(y="Slope value", x="Index", fill="Cluster")+
  facet_grid(k_cluster~.)+
  scale_x_discrete(labels=c("slope_Dim.2"=expression(Delta~"TraitDim2"),
                            "slope_Dim.1"=expression(Delta~"TraitDim1"),
                            "beta_coef_fori"=expression(Delta~"FOri"),
                            "beta_coef_fdiv"=expression(Delta~"FDiv"),
                            "beta_coef_feve"=expression(Delta~"FEve"),
                            "lm_coef_sesfric"=expression(Delta~"cFRic")))+
  scale_fill_viridis(discrete = T, option = "D")+
  #scale_fill_viridis_d(option = "turbo")+
  theme_classic()+
  theme (strip.background = element_blank(),
         strip.text = element_blank(),
         axis.text.y = element_text(size=13, angle=0, family = "sans"),
         axis.title.y = element_text(size=19, angle=90, family = "sans"),
         axis.text.x   = element_text(size=13, angle=0, family = "sans"),
         axis.title.x  = element_blank())


# Map using NEON Domains
df_domains <- ddply (db_domains,. (DomainName), summarise,
                     Lat=mean(Lat), Long=mean(Long))
# merge cluster_db and domains including species richness slope
cluster_domains_full <- left_join(cluster_db_full, db_domains)
cluster_domains_full$k_cluster_full <- paste ("k",cluster_domains_full$k_cluster_full, sep="_")
cluster_domains_full$presence <- 1
cluster_domains_full <- ddply (cluster_domains_full,. (DomainName,k_cluster_full),
                               summarise, nCluster=sum(presence))
cluster_domains_full <- dcast(DomainName~k_cluster_full, 
                              value.var = "nCluster",
                              data = cluster_domains_full)

cluster_domains_full [is.na(cluster_domains_full)] <- 0
cluster_domains_full <- left_join (cluster_domains_full, df_domains)
cluster_domains_full <- cluster_domains_full[,c("DomainName","Lat","Long","k_1","k_2","k_3",'k_4','k_5')]
cluster_domains_full$radius <- (apply (cluster_domains_full[,4:length(cluster_domains_full)],1,sum))/100

#CRS projection 
crs_projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

# upload data
neonDomains <- readOGR(dsn="NEONDomains_0/", layer="NEON_Domains")
neonDomains <- spTransform(neonDomains, CRSobj = crs_projection)

#setdiff(cluster_domains$DomainID,neonDomains@data$DomainID)
# Transform object under same map projection
cluster_db_full$k_cluster_full <- as.factor(cluster_db_full$k_cluster_full)

coord <- cbind(cluster_db_full$Long, cluster_db_full$Lat)
prj.coord  <- project(coord, proj = crs_projection)
coords_map <- cbind (prj.coord,slopes_f_indices_coord)
names(coords_map)[1:2] <- c("X.prj","Y.prj") 

mapDataInfo <- map_data (neonDomains)


Fig_2b <- ggplot () +
  geom_polygon(data=neonDomains, aes(x=long, y=lat, group=group, colour=group), 
               color="black", fill="grey90") +
  geom_scatterpie(data=cluster_domains_full,aes(x=Long,y=Lat, group=DomainName, r=radius),
                  alpha=1, color=NA,
                  cols =c("k_1",'k_2','k_3','k_4','k_5'))+
  geom_scatterpie_legend(cluster_domains_full$radius, x=-123, y=30)+
  lims(x=c(-130,-60), y=c(24,50)) +
  labs(fill="Cluster")+
  #scale_fill_viridis_d(option = "turbo")+
  scale_fill_viridis(discrete = T, option = "D")+
  theme_void()


ggarrange(Fig_2b,Fig_2a,
          align="hv", widths = c(2.3,0.7),
          ncol=2, nrow=1,
          font.label = list(size=13, color="black",
                            face="bold", family = "sans"))


# 7. Environmental factors maps -------------------------------------------


elevation  <- read.csv ("DataInter/elevation_data.csv", header = T, sep = ",")
climate_db <- read.csv ("DataInter/Environmental_data.csv", header = T, sep = ",")
setdiff(elevation$Location, climate_db$Location)

# Merging elevation and climate database
env_factors <- merge (climate_db, elevation[,c("Abbrev","Elevation")],by="Abbrev")
unique (is.na(env_factors$Abbrev))
unique (is.na(env_factors$Lat))
unique (is.na(env_factors$Location))
unique (is.na(env_factors$DomainID))
unique (is.na(env_factors$DomainName))
unique (is.na(env_factors$Elevation))

# IMPORTANT
# Some values of elevation column are missing.  
table (is.na (unique (env_factors$Elevation))) 
env_factors <- na.omit(env_factors)

# Extracting slope values from abiotic factors 
env_factors_years <- split (env_factors, env_factors$Abbrev)

slope_meanTemp  <- list ()
slope_maxTemp   <- list ()
slope_minTemp   <- list ()
slope_ppt       <- list ()

# Loop to extracting estimate value from each linear regression - function was "coefficients"

for (m in 1:length(env_factors_years)){
  slope_meanTemp[[m]]  <- summary(lm (meanTemp~Year, data = env_factors_years[[m]]))$coefficients[2,1]
  slope_maxTemp[[m]]   <- summary(lm (maxTemp~Year, data = env_factors_years[[m]]))$coefficients[2,1]
  slope_minTemp[[m]]   <- summary(lm (minTemp~Year, data = env_factors_years[[m]]))$coefficients[2,1]
  slope_ppt[[m]]       <- summary(lm (ppt~Year, data = env_factors_years[[m]]))$coefficients[2,1]
  
}
slopes_env_factors <- data.frame(Abbrev=unique(env_factors$Abbrev),
                                 coef_meanTemp= do.call("rbind",slope_meanTemp),
                                 coef_maxTemp= do.call("rbind",slope_maxTemp),
                                 coef_minTemp= do.call("rbind",slope_minTemp),
                                 coef_ppt= do.call("rbind",slope_ppt),
                                 elev= env_factors$Elevation)

rm (slope_meanTemp,slope_maxTemp,slope_minTemp,slope_ppt)

abioticFactor_coord <- left_join(Location.Coords,slopes_env_factors)
abioticFactor_coord <- na.omit(abioticFactor_coord)

# This function was created to shown the maps considering environmental factors 
map_abiotic <- function (db, f_index, NameIndex, GridCellReso, HighColor, LowColor, MidColor){
  usa <- map_data("usa")
  db <- db
  db <- data.frame (x=db$Long,
                    y=db$Lat,
                    value=f_index)
  
  #Transform in raster data.frame
  r <- rasterFromXYZ(db, crs = CRS("+proj=longlat +datum=WGS84"), res = GridCellReso, digits = 0)
  
  
  # create an extent object for the extent of the United States
  usa_bbox <- c(-125.5, 24.5, -67, 49.5)
  usa_extent <- extent(usa_bbox)
  
  # crop the raster to the extent of the world map
  raster_cropped <- crop(r, extent(usa_extent))
  # convert the cropped raster to a data frame
  df <- data.frame(rasterToPoints(raster_cropped))
  
  
  # Create a raster plot of the variable/index
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = value))+
    geom_polygon(data=usa, aes(x=long, y=lat, group=group),
                 colour="black", fill=NA)+
    labs(fill=NameIndex) +
    scale_fill_gradient2 (low = LowColor, high = HighColor, mid = MidColor) +
    coord_quickmap() + 
    theme_void()
  
}

Fig_S1a <- map_abiotic(db =abioticFactor_coord,
                       f_index = abioticFactor_coord$coef_maxTemp,
                       NameIndex = expression(Delta~"Max.Temp"),
                       GridCellReso = .7,
                       HighColor = brewer.pal(9, "YlOrBr"),
                       LowColor = rev(brewer.pal(9, "Blues")),
                       MidColor = "white")
Fig_S1b <- map_abiotic(db =abioticFactor_coord,
                       f_index = abioticFactor_coord$coef_meanTemp,
                       NameIndex = expression(Delta~"Mean.Temp"),
                       GridCellReso = .7,
                       HighColor = brewer.pal(9, "YlOrBr"),
                       LowColor = rev(brewer.pal(9, "Blues")),
                       MidColor = "white")

Fig_S1c <- map_abiotic(db =abioticFactor_coord,
                       f_index = abioticFactor_coord$coef_minTemp,
                       NameIndex = expression(Delta~"Min.Temp"),
                      GridCellReso = .7,
                      HighColor = brewer.pal(9, "YlOrBr"),
                      LowColor = rev(brewer.pal(9, "Blues")),
                      MidColor = "white")

Fig_S1d <- map_abiotic(db =abioticFactor_coord,
                       f_index = abioticFactor_coord$coef_ppt,
                       NameIndex = expression(Delta~"Precipitation"),
                       GridCellReso = .7,
                       HighColor = brewer.pal(9, "YlGnBu"),
                       LowColor = rev(brewer.pal(9, "YlOrBr")),
                       MidColor = "white")

Fig_S1e <- map_abiotic(db =abioticFactor_coord,
                       f_index = abioticFactor_coord$elev,
                       NameIndex = "Elevation",
                       GridCellReso = .7,
                       HighColor = brewer.pal(9, "Oranges"),
                       LowColor = rev(brewer.pal(9, "Greys")),
                       MidColor = "white")

ggarrange(Fig_S1a, Fig_S1d, 
          Fig_S1b, Fig_S1e, 
          Fig_S1c,
          align="hv",
          ncol=2, nrow=3,
          font.label = list(size=13, color="black",
                            face="bold", family = "sans"))



# 8. Figure 3: Trait dimensions maps --------------------------------------


res <- readRDS("DataInter/MFA_Results.rds")
cwm_full <- readRDS("DataInter/cwm_full.rds")

##### Plot MFA results #####
temporal_cwm_dim <- data.frame (cbind(Id=rownames(cwm_full),
                                      res$ind$coord[,1:2]))
temporal_cwm_dim <- separate (temporal_cwm_dim, col = Id, 
                              into=c("Abbrev","year"), sep="_")
rownames(temporal_cwm_dim) <- NULL
temporal_cwm_dim$Dim.1 <- as.numeric (temporal_cwm_dim$Dim.1)
temporal_cwm_dim$Dim.2 <- as.numeric (temporal_cwm_dim$Dim.2)
temporal_cwm_dim$year  <- as.numeric (temporal_cwm_dim$year)

# Slope Dim.1 and Dim.2 
dimensions_years <- split (temporal_cwm_dim, temporal_cwm_dim$Abbrev)

slope_Dim.1 <- list ()
slope_Dim.2 <- list ()

# Loop to extracting estimate value from each linear regression - function was "coefficients"
for (m in 1:length(dimensions_years)){
  slope_Dim.1[[m]] <- summary(lm (Dim.1~year, data = dimensions_years[[m]]))$coefficients[2,1]
  slope_Dim.2[[m]] <- summary(lm (Dim.2~year, data = dimensions_years[[m]]))$coefficients[2,1]
  
}
slopes_dimensions <- data.frame(Abbrev=unique(temporal_cwm_dim$Abbrev),
                                slope_Dim.1= do.call("rbind",slope_Dim.1),
                                slope_Dim.2= do.call("rbind",slope_Dim.2))
slopes_dimensions <- left_join(slopes_dimensions, Location.Coords)

# This function was created to shown the trait dimensions maps
map_function2 <- function (db, f_index, NameIndex, GridCellReso){
  usa <- map_data("usa")
  db <- db
  db <- data.frame (x=db$Long,
                    y=db$Lat,
                    value=f_index)
  
  #Transform in raster data.frame
  r <- rasterFromXYZ(db, crs = CRS("+proj=longlat +datum=WGS84"), res = GridCellReso, digits = 0)
  
  
  # create an extent object for the extent of the United States
  usa_bbox <- c(-125.5, 24.5, -67, 49.5)
  usa_extent <- extent(usa_bbox)
  
  # crop the raster to the extent of the world map
  raster_cropped <- crop(r, extent(usa_extent))
  # convert the cropped raster to a data frame
  df <- data.frame(rasterToPoints(raster_cropped))
  
  
  # Create a raster plot of the variable/index
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = value))+
    geom_polygon(data=usa, aes(x=long, y=lat, group=group),
                 colour="black", fill=NA)+
    labs(fill=NameIndex) +
    # scale_fill_gradient (low = LowColor, high = HighColor, mid = MidColor, alpha=1) +
    scale_fill_viridis(option = "H")+
    coord_quickmap() + 
    theme_void()
  
}

# Map using slopes Dim.1 and Dim.2
Fig_3a <- map_function(db =slopes_dimensions,
                       f_index = slopes_dimensions$slope_Dim.1,
                       NameIndex = expression(Delta~"Dim.1"),
                       GridCellReso=0.7,
                       HighColor = brewer.pal(9, "RdPu"),
                       LowColor = rev(brewer.pal(6, "Blues")),
                       MidColor = "white")


Fig_3b <- map_function (db = db2,
                        f_index = db2$slope_Dim.2,
                        NameIndex = expression(Delta~"Dim.2"),
                        GridCellReso=0.7, 
                        HighColor = brewer.pal(9, "RdPu"),
                        LowColor = rev(brewer.pal(6, "Blues")),
                        MidColor = "white")

# Linear correlations 
temporal_dim_domains <- left_join(temporal_cwm_dim, db_domains)
dim_domains <- split (temporal_dim_domains, temporal_dim_domains$DomainName)

results_lm_domains_dim1 <- list()
results_lm_domains_dim2 <- list()

for (i in 1:length(dim_domains)){
  # Dim.1
  lm_domain_dims <- lm (Dim.1~year, data = dim_domains[[i]])
  db <- as.data.frame (dim_domains[[i]]) 
  predict_domain_dim1 <- data.frame (lm_predict=predict(lm_domain_dims,db),
                                     year=db$year,
                                     Domain="Domain")
  results_lm_domains_dim1[[i]] <- predict_domain_dim1
  
  
  # Dim.2
  lm_domain_dims <- lm (Dim.2~year, data = dim_domains[[i]])
  db <- as.data.frame (dim_domains[[i]]) 
  predict_domain_dim2 <- data.frame (lm_predict=predict(lm_domain_dims,db),
                                     year=db$year,
                                     Domain="Domain")
  results_lm_domains_dim2[[i]] <- predict_domain_dim2
  
}

linear_figure_function2 <- function (dbFull, findex, model_equation,dbPre, yNameAxis){
  
  lm_findex_fit <- model_equation
  predicted_df  <- data.frame (findex_pred=predict(lm_findex_fit,dbFull), 
                               year=dbFull$year, Location="Local") 
  
  ggplot (dbFull, aes(x=year, y=findex, shape=Abbrev)) +
    stat_smooth(method = lm, formula = y~x, se = F,fullrange=TRUE, colour="grey90", size=0.5)+
    geom_line (color="black",data =predicted_df, aes(x=year, y=findex_pred, shape=Location))+
    geom_line (color="#800000",data=dbPre[[1]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#DC143C",data=dbPre[[2]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#FF7F50",data=dbPre[[3]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#FFA500",data=dbPre[[4]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#BDB76B",data=dbPre[[5]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#808000",data=dbPre[[6]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#9ACD32",data=dbPre[[7]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#006400",data=dbPre[[8]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#32CD32",data=dbPre[[9]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#20B2AA",data=dbPre[[10]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#AFEEEE",data=dbPre[[11]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#4682B4",data=dbPre[[12]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#00008B",data=dbPre[[13]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#4B0082",data=dbPre[[14]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#8B008B",data=dbPre[[15]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#DA70D6",data=dbPre[[16]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    geom_line (color="#BC8F8F",data=dbPre[[17]], aes(x=year, y=lm_predict, shape=Domain),linetype="dashed")+
    labs (x="Years", y=yNameAxis) +
    theme_classic() 
  
}

Fig_3.1a <- linear_figure_function2 (dbFull = temporal_dim_domains,
                                  findex = temporal_dim_domains$Dim.1,
                                  model_equation = lm(Dim.1~year, data =temporal_dim_domains),
                                  dbPre = results_lm_domains_dim1,
                                  yNameAxis = expression("Dim.1"))


Fig_3.1b <- linear_figure_function2 (dbFull = temporal_dim_domains,
                                   findex = temporal_dim_domains$Dim.2,
                                   model_equation = lm(Dim.1~year, data =temporal_dim_domains),
                                   dbPre = results_lm_domains_dim2,
                                   yNameAxis = expression("Dim.2"))

ggarrange(Fig_3a, Fig_3.1a,
          Fig_3b, Fig_3.1b,
          align="hv", widths = c(2,0.5),
          ncol=2, nrow=2,
          font.label = list(size=13, color="black",
                            face="bold", family = "sans"))


# Total contribution set traits
contrib1 <- facto_summarize(res, element = "group", result = "contrib", axes=1)
contrib1$name <- as.factor(c("Morphologic","Mass","Diet","FStrategy"))


cont_groupsd1 <- ggplot(contrib1, aes(name, contrib))+
  ylim (0, 50) +
  geom_bar(stat="identity", position = "identity",fill=c("#999999","#4682B4","#9999CC", "#66CC99"))+
  geom_hline(yintercept = 25, linetype = 2,color = "red")+
  theme_bw()+xlab("")+ylab("")+
  ggtitle("Contribution (%) to Dim1")

contrib2 <- facto_summarize(res, element = "group", result = "contrib", axes=2)
contrib2$name <- as.factor(c("Morphologic","Mass","Diet","FStrategy"))

cont_groupsd2 <- ggplot(contrib2, aes(name, contrib))+
  ylim (0, 50) +
  geom_bar(stat="identity", position = "identity",fill=c("#999999","#4682B4","#9999CC", "#66CC99"))+
  geom_hline(yintercept = 30, linetype = 2,color = "red")+
  theme_bw()+xlab("")+ylab("")+
  ggtitle("Contribution (%) to Dim2")

ggarrange(dynamic_dim1, dynamic_dim2,
  cont_groupsd1, cont_groupsd2,
  align="hv",ncol=2, nrow=2,
  font.label = list(size=13, color="black",face="bold", family = "sans"))

# by each trait (x=Dimension 1 ; y=Dimension 2)
morp  <- 7
mass  <- 1
diet  <- 10
fstra <- 7
imptraits <- data.frame(trait = rownames(res$quanti.var$coord),
                        Dim.1 = res$quanti.var$coord[,1],
                        Dim.2 = res$quanti.var$coord[,2],
                        col1=c(rep(1,morp),rep(2,mass),rep(3,diet),rep(4,fstra)),
                        typeTrait=c(rep("Morp",7),rep("Mass",1),
                                    rep("Diet", 10),rep("ForStr",7)))

imptraits$trait <- ifelse(imptraits$trait=="Beak.Length_Culmen","BLC", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Beak.Length_Nares","BLN", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Beak.Width","BWI", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Beak.Depth","BDE", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Tarsus.Length","TAL", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Hand.Wing.Index","HWI", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Tail.Length","TLE", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Mass","MASS", imptraits$trait)

imptraits$trait <- ifelse(imptraits$trait=="Diet.Inv","INV", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Diet.Vend","END", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Diet.Vect","ECT", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Diet.Vfish","FISH", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Diet.Vunk","UNK", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Diet.Scav","SCA", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Diet.Fruit","FRU", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Diet.Nect","NECT", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Diet.Seed","SEED", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="Diet.PlantO","PLAN", imptraits$trait)

imptraits$trait <- ifelse(imptraits$trait=="ForStrat.watbelowsurf","WB.SUR", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="ForStrat.wataroundsurf","WA.SUR", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="ForStrat.ground","GROUND", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="ForStrat.understory","UNDERS", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="ForStrat.midhigh","MHIGH", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="ForStrat.canopy","CANOPY", imptraits$trait)
imptraits$trait <- ifelse(imptraits$trait=="ForStrat.aerial","AERIAL", imptraits$trait)

unique (imptraits$trait)
p_trait1 <- ggplot(imptraits %>% arrange(col1, Dim.1),
                   aes(y=factor(trait, levels = c(# Morp.
                                                  "HWI","BLC","TAL","BLN", 
                                                  "BWI","BDE","TLE",       
                                                  # Mass
                                                  "MASS",                  
                                                  # Diet
                                                  "FRU","SEED","END","NECT",
                                                  "FISH","UNK","INV","ECT", 
                                                  "SCA","PLAN", 
                                                  # Niche
                                                  "UNDERS","MHIGH","GROUND","CANOPY",
                                                  "AERIAL","WA.SUR", 
                                                  "WB.SUR")),        
                       x=Dim.1, colour = as.factor(col1))) +
  geom_point(size=3, shape = 16, stroke = 1.5, alpha = 0.9, show.legend = FALSE) +
  geom_segment(aes(x=0, xend=Dim.1, y=trait, yend=trait))+
  scale_x_continuous(limits = c(-1,1)) +
  geom_vline(xintercept = 0, linetype = 3, size = .5, color = "black", alpha = 1) +
  scale_color_manual(values = c("#999999","#4682B4","#9999CC","#66CC99")) +
  labs (x="Dim 1", y="Trait") +
  theme_classic()+
  theme (legend.position = "na",
         axis.title.y = element_blank())


# gggap (plot=p_trait1, 
#        segments = c(1,8),
#        ylim=c(1,25))

p_trait2 <- ggplot(imptraits %>% arrange(col1, Dim.2),
                   aes(y=factor(trait, levels = c(# Morp.
                                                  "TLE","BDE","BWI","HWI",
                                                  "BLN","BLC","TAL",
                                                  # Mass
                                                  "MASS",
                                                  # Diet
                                                  "INV","ECT","FISH","UNK","SCA",
                                                  "NECT","PLAN","END","FRU","SEED",
                                                  # Niche
                                                  "WA.SUR","AERIAL","UNDERS","CANOPY",
                                                  "GROUND","MHIGH","WB.SUR")),
                       x=Dim.2, colour = as.factor(col1))) +
  geom_point(size=3, shape = 16, stroke = 1.5, alpha = 0.9, show.legend = FALSE) +
  geom_segment(aes(x=0, xend=Dim.2, y=trait, yend=trait))+
  scale_x_continuous(limits = c(-1,1)) +
  geom_vline(xintercept = 0, linetype = 3, size = .5, color = "black", alpha = 1) +
  scale_color_manual(values = c("#999999","#4682B4","#9999CC","#66CC99")) +
  labs (x="Dim 2", y="Trait") +
  theme_classic()+
  theme (legend.position = "na",
         axis.title.y = element_blank())


ggarrange(p_trait1, p_trait2,
          align="hv",ncol=1, nrow=2,
          font.label = list(size=13, color="black",face="bold", family = "sans"))



# 9. Principal component figure -------------------------------------------

tempPCA2 <- princomp(slopes_env_factors[,2:4], cor = F)
summary (tempPCA2)
tempPCA2$loadings


Temp_PRCOMP <- data.frame (loadings(tempPCA2))
Temp_PRCOMP <- as.data.frame(Temp_PRCOMP[1:3, ])

Temp_PRCOMP$Metric <- c("Mean Temp",
                        "Max. Temp",
                        "Min. Temp")
rownames(Temp_PRCOMP)<- NULL
Temp_PRCOMP <- melt (Temp_PRCOMP, id=("Metric"))
Temp_PRCOMP <- droplevels(Temp_PRCOMP[!Temp_PRCOMP$variable=="Comp.3", ])

PC_figure <- ggplot (Temp_PRCOMP, aes(y=value, fill=Metric,
                         x=variable))+
  geom_bar(stat = "identity", color="black", position = position_dodge())+
  scale_fill_manual(values = c("grey30","grey60","grey90"))+
  labs(x="Principal component axis", fill="Temperature\nmetric")+
  scale_x_discrete(breaks=c("Comp.1", "Comp.2"), labels=c("PC1 (81%)","PC2 (18%)"))+
  theme_classic2()+
  theme (axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text = element_text(size = 13))


# Figure 4: Bayesian models -----------------------------------------------

mRich <- readRDS("DataInter/mRich.rds")
mcFRic <- readRDS("DataInter/mSESFRic.rds")
mFEve <- readRDS("DataInter/mFEve.rds")
mFDiv <- readRDS("DataInter/mFDiv.rds")
mFOri <- readRDS("DataInter/mFOri.rds")
mDim1 <- readRDS("DataInter/mDim1.rds")
mDim2 <- readRDS("DataInter/mDim2.rds")

results_models <- rbind (mcmc_intervals_data(mRich, prob = 0.8, 
                    prob_outer = 0.95, 
                    pars=c("b_Comp.1",
                           "b_Comp.2", 
                           "b_coef_ppt_scale",
                           "b_elevation_scale")),
       mcmc_intervals_data(mcFRic, prob = 0.8, 
                           prob_outer = 0.95, 
                           pars=c("b_Comp.1",
                                  "b_Comp.2", 
                                  "b_coef_ppt_scale",
                                  "b_Elevation")),
       mcmc_intervals_data(mFEve, prob = 0.8, 
                           prob_outer = 0.95, 
                           pars=c("b_Comp.1",
                                  "b_Comp.2", 
                                  "b_coef_ppt_scale",
                                  "b_elevation_scale")),
       mcmc_intervals_data(mFDiv, prob = 0.8, 
                           prob_outer = 0.95, 
                           pars=c("b_Comp.1",
                                  "b_Comp.2", 
                                  "b_coef_ppt_scale",
                                  "b_elevation_scale")),
       mcmc_intervals_data(mFOri, prob = 0.8, 
                           prob_outer = 0.95, 
                           pars=c("b_Comp.1",
                                  "b_Comp.2", 
                                  "b_coef_ppt_scale",
                                  "b_elevation_scale")),
       mcmc_intervals_data(mDim1, prob = 0.8, 
                           prob_outer = 0.95, 
                           pars=c("b_Comp.1",
                                  "b_Comp.2", 
                                  "b_coef_ppt_scale",
                                  "b_Elevation")),
       mcmc_intervals_data(mDim2, prob = 0.8, 
                           prob_outer = 0.95, 
                           pars=c("b_Comp.1",
                                  "b_Comp.2", 
                                  "b_coef_ppt_scale",
                                  "b_Elevation")))

results_models <- data.frame (results_models)
results_models$index <- c (rep("SRic",4),rep("cFRic",4),rep("FEve",4),rep("FDiv",4),
                              rep("FOri",4),rep("Dim.1",4),rep("Dim.2",4))

results_models$parameter <- gsub("b_Elevation","b_elevation_scale",results_models$parameter)
scale(results_models$ll)  
  
  
ggplot(results_models,  aes(x = scale(m), y = factor(parameter, levels = c("b_Comp.1","b_Comp.2",
                                                                    "b_coef_ppt_scale",
                                                                    "b_elevation_scale")), 
                                    color = ifelse(scale(ll) < 0 & scale(hh) < 0, '#990033',    
                                            ifelse(scale(ll) > 0 & scale(hh) > 0, '#3300CC', 
                                                   'grey')))) +
  scale_y_discrete(labels=c("b_Comp.1"="PC.1",
                            "b_Comp.2"="PC.2",
                            "b_coef_ppt_scale"=expression(Delta~"PRT"),
                            "b_elevation_scale"="ELEV"))+
  scale_color_identity()+
  facet_grid(~factor(index, levels = c("SRic","cFRic","FEve","FDiv","FOri","Dim.1","Dim.2")), 
             scales = "free")+
  geom_linerange(aes(xmin = scale(l), xmax = scale (h)), size=2)+ # adds internal confidence  
  geom_linerange(aes(xmin = scale(ll), xmax = scale(hh))) + # adds outer confidence   #
  theme_bw()+
  theme (axis.title.y = element_blank(),
         axis.title.x = element_blank())