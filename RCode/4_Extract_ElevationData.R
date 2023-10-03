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

library(elevatr)
library(sf)
library(raster)

# 1. Call location coordinates data  --------------------------------------

geoCoorLoc <- read.csv("DataInter/RegionalPool_coord.csv", header = T, sep=",")


# 2. Projection map -------------------------------------------------------

ll_proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
elev <- elevatr::get_elev_point(pt_df, prj = ll_proj)

coordinates(geoCoorLoc) <- c('Long', 'Lat')
proj4string(geoCoorLoc) <- CRS(ll_proj)


# 3. Extracting elevation data  -------------------------------------------

elevation_data <- elevatr::get_elev_raster(locations = geoCoorLoc, z = 9, clip = "locations")
elevation_data <- data.frame (coordinates(geoCoorLoc), Abbrev=geoCoorLoc$Abbrev,
                                Location=geoCoorLoc$Location, DomainID=geoCoorLoc$DomainID,
                                DomainName=geoCoorLoc$DomainName,
                                extract(elevation_data, geoCoorLoc))
colnames (elevation_data)[7] <- "Elevation"


                                