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

# This code was created to extract environmental factors from PRISM Climate data repository within radio of 24 km


rm(list = ls())

# 1. Call location coordinates data ---------------------------------------

geoCoorLoc <- read.csv("DataInter/RegionalPool_coord.csv", header = T, sep=",")



# 2. Extracting data in a radio 24 km -------------------------------------

Ext_Abiotic_Factors <- function (WayFile, AbiFactor, YearsPer, MonthPer, coorDB, NameFactor) {
  library(tidyverse) # a suite of packages for wrangling and tidying data
  library(prism)     # package to access and download climate data
  library(raster)    # the climate data comes in raster files- this package helps process those
  library(stringr)
  library(magrittr)
  library(reshape)
  
  # First, set a file path where prism data will be stored
  options(prism.path = WayFile)
  
  
  # Now, select the type of data (mean temperature and precipitation for us) and date range
  get_prism_monthlys(type = AbiFactor, years = YearsPer, mo=MonthPer, keepZip = F, keep_pre81_months = F)

  # Grab the prism data and compile the files
  climate_data <- prism_archive_ls() %>% 
    pd_stack(.)
  
  # Extract project coordinates from raster stack
  climate_crs <- climate_data@crs@projargs
  
  # Convert these locations to format that can be matched to Prism climate data
  coordinates(coorDB) <- c('Long', 'Lat')
  proj4string(coorDB) <- CRS(climate_crs)
  
  
  # Join Climate & Data
  # Extract the  extracted from the raster stack for those sites 
  climate_location <- data.frame (coordinates(coorDB), Abbrev=coorDB$Abbrev,
                                  Location=coorDB$Location, DomainID=coorDB$DomainID,
                                  DomainName=coorDB$DomainName,
                                  raster::extract(climate_data, coorDB,
                                          buffer=24000,
                                          fun=mean))
  
  # Reshape data
  climate_location <- melt (climate_location, id=c("Long","Lat","Abbrev","Location","DomainID","DomainName"))
  
  # Split header column
  climate_location <- separate(climate_location, "variable",
                               into = c("Repository","AbioFactor","Lixo","Resol","Year","Lixo2"),
                               sep = "_")
  
  # Rename columns
  climate_location <- climate_location[,c("Long","Lat","Abbrev","Location","DomainID","DomainName",
                                          "Resol","Year","value")]
  climate_location <- separate(climate_location, "Year",
                               into = c("Year","Month"),
                               sep=4)
  
  
  colnames(climate_location)[10] <- NameFactor
  
  print (climate_location)
  
  }

meanTemp <- Ext_Abiotic_Factors (WayFile = "Prism_ClimateData/Temp/meanTemp",
                                 AbiFactor = "tmean",
                                 YearsPer = 1967:2018,
                                 MonthPer = c(1,11,12),
                                 coorDB = geoCoorLoc,
                                 NameFactor = "meanTemp")
maxTemp <- Ext_Abiotic_Factors (WayFile =  "Prism_ClimateData/Temp/maxTemp",
                                AbiFactor = "tmax",
                                YearsPer = 1967:2018,
                                MonthPer = c(1,11,12),
                                coorDB = geoCoorLoc,
                                NameFactor = "maxTemp")
minTemp <- Ext_Abiotic_Factors (WayFile = "Prism_ClimateData/Temp/minTemp",
                                AbiFactor = "tmin",
                                YearsPer = 1967:2018,
                                MonthPer = c(1,11,12),
                                coorDB = geoCoorLoc,
                                NameFactor = "minTemp")
precipitation <- Ext_Abiotic_Factors (WayFile = "Prism_ClimateData/Pres",
                                      AbiFactor = "ppt",
                                      YearsPer = 1967:2018,
                                      MonthPer = 1:12,
                                      coorDB = geoCoorLoc,
                                      NameFactor = "Prec")


# IMPORTANT --------------------------------------------------------
# All temperature data were extracted only during winter months (Nov,Dec,Jan).
# The precipitation data were extracted all months in each year.


# 3. Joint environmental data ---------------------------------------------

library (dplyr)
library (plyr)
meanTemp_db <- ddply (meanTemp,. (Long, Lat, Abbrev, Location, DomainID, DomainName, Resol, Year),
                      summarise, meanTemp=mean(meanTemp))

maxTemp_db  <- ddply (maxTemp,. (Long, Lat, Abbrev, Location, DomainID, DomainName, Resol, Year),
                      summarise, maxTemp=mean(maxTemp))
  
minTemp_db  <- ddply (minTemp,. (Long, Lat, Abbrev, Location, DomainID, DomainName, Resol, Year),
                      summarise, minTemp=mean(minTemp))

sumValues_Year <- function (x){
  x <- na.omit(x)
  x <- plyr::ddply(x,. (Long, Lat, Abbrev, Location, DomainID, DomainName, Resol, Year),
                   summarise, ppt=sum(Prec))
}
precipitation_db  <- sumValues_Year(precipitation)


# Integrating environmental factors

climate_db <- cbind(meanTemp_db,
                    maxTemp=maxTemp_db$maxTemp,
                    minTemp=minTemp_db$minTemp,
                    ppt=precipitation_db$ppt)

head (climate_db)

