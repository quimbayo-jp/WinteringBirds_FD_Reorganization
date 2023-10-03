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

# clean up
rm(list=ls())

# Libraries  
library (reshape2)
library (mFD)
library (FD)
library (gawdis)
library (plyr)
library (tidyr)
library (FactoMineR)
library (stringr)

# Specific functions 
source("RCode/QFunctionalSpace.R")

# 1. Calling data ------------------------------------
AbunYear  <-readRDS("DataRaw/commMats_byYear.rds")
TraitSPP  <-read.csv("DataInter/Traits_Birds.csv", stringsAsFactors = T, 
                     dec = ".", sep = ",", header = T, row.names = 1)
TraitSPP$Mass <- log10(TraitSPP$Mass)


# IMPORTANT
# Overall 296 bird species was observed over time. 
# Total sites: 1183
# Total years:52 (1967-2018)

# Name comparison between TraitDataset and AbundDataset
setdiff(rownames(TraitSPP), colnames(AbunYear[[1]][,2:length(AbunYear[[1]])]))

# 2. Exclude Nocturnal and pelagic species from Trait database ------------
TraitSPP <- droplevels (subset(TraitSPP, Nocturnal==0)) # Ten species were removed from this dataset. 
TraitSPP <- droplevels (subset(TraitSPP, PelagicSpecialist==0)) 
TraitSPP <- TraitSPP[,1:28] # Nocturnal and Pelagic.Specialist columns were removed
TraitSPP <- TraitSPP[,c(1:5,9:length(TraitSPP))]

# 3. Estimating functional indices per location in each year -----------
findices_estimated    <- list ()
traits_weights_values <- matrix(NA, ncol =length(AbunYear), nrow = 25)

for (j in 1:length(AbunYear)) {
  
  Abun_Year <- AbunYear[[j]]
  
  # Removing Nocturnal species from Abundance estimation dataset
  Abun_Year <- melt(Abun_Year, id=c("Circle"))
  colnames (Abun_Year) <- c("Abbrev","Name","Abun")
  Abun_Year <- droplevels(Abun_Year[Abun_Year$Name %in% rownames(TraitSPP), ])
  Abun_Year <- dcast (Abbrev~Name, value.var = "Abun", data = Abun_Year)
  
  # Preparing dataset
  rownames(Abun_Year) <- Abun_Year$Abbrev
  Abun_Year <- Abun_Year[,2:length(Abun_Year)]
  
  # Dissimilarity matrix using trait dataset
  # Each number represent one type of group (i.e. Quantitative, Ordinal)
  # Beak length = 1 (rep 4 times)
  # Tarsus.Length = 2
  # Hand.Wing.Index = 3
  # Tail.Length = 4
  # Mass = 5
  # Diet = 6 (rep 10 times)
  # Foraging strategies = 7 (rep 7 times)
  
  traitGroups <- c(rep(1,4),2,3,4,5,rep(6,10),rep(7,7))
  # The fuzzy groups are all those traits that are grouped in several colm
  fuzzyGroups <- c(1,6,7)
  
  ### The next step computing the multi-trait dissimilarity. 
  #The approach is based on minimizing the differences in the correlation between the 
  #dissimilarity of each trait, or groups of traits, and the multi-trait dissimilarity. 
  #This is done using either an analytic or a numerical solution, both available in 
  #the function.
  
  traits.dist <- gawdis(TraitSPP, w.type = "optimized", groups = traitGroups, fuzzy = fuzzyGroups, opti.maxiter=300) #300 is the default
  traits_weights <- attr(traits.dist, "weights")
  traits_weights_values[,j] <- traits_weights
  
  fspaces_quality_birds <- mFD::quality.fspaces(sp_dist = traits.dist,
                                                maxdim_pcoa = 10,
                                                deviation_weighting = "absolute",
                                                fdendro = "average")
  sp_faxes_coord_birds <- fspaces_quality_birds$"details_fspaces"$"sp_pc_coord"
  
  # Estimation Functional indices 
  fd_indices_birds     <- mFD::alpha.fd.multidim(sp_faxes_coord = sp_faxes_coord_birds[,c("PC1","PC2","PC3","PC4")],
                                                 asb_sp_w = data.matrix(Abun_Year),
                                                 scaling = TRUE,
                                                 check_input = TRUE,
                                                 details_returned = TRUE)
  
  
  # notify (title = "Functional indices estimation",
  #        msg = c("Finish year", j))
  
  #system("say finished one year!")
  findices_estimated[[j]] <- fd_indices_birds$functional_diversity_indices
  
}

# Security copy 
# copy_findices_estimated <- findices_estimated 
# saveRDS(findices_estimated, "DataInter/findices_estimated.RDS")

# Add Year and location in each list which contain the functional indices estimated
for (i in 1:length(findices_estimated)){
  findices_estimated[[i]]$year <- 1967+i
  findices_estimated[[i]]$Location <- rownames (findices_estimated[[i]])
}
DB_FIndex   <- do.call("rbind", findices_estimated)
DB_FIndex$ID_unique <- paste (DB_FIndex$Location, DB_FIndex$year, sep="_")

#saveRDS(DB_FIndex, "DataInter/DB_FIndex.rds")

head (traits_weights_values)
traits_weights_values <- data.frame (traits_weights_values)
traits_weights_values$trait <- names (traits_weights)
traits_weights_values$mean_value <- apply (traits_weights_values[,1:52],1,mean)
# traits_weights_values[,53:54]
# sum (traits_weights_values[19:25,54])

#saveRDS(DB_FIndex, "DataInter/DB_TraitsWeights.rds")


# 4. Estimating CWM -------------------------
cwm_results <- list ()

for (j in 1:length(AbunYear)) {
  
  Abun_Year <- AbunYear[[j]]
  
  # Removing Nocturnal species from Abundance estimation dataset
  Abun_Year <- melt(Abun_Year, id=c("Circle"))
  colnames (Abun_Year) <- c("Abbrev","Name","Abun")
  Abun_Year <- droplevels(Abun_Year[Abun_Year$Name %in% rownames(TraitSPP), ])
  
  Abun_Year$Name <- tolower(Abun_Year$Name) 
  Abun_Year$Name <- gsub (" ", "_", Abun_Year$Name)
  Abun_Year$Abun <- log (Abun_Year$Abun+1)
  Abun_Year <- Abun_Year[order(Abun_Year$Name, decreasing = F), ]
  Abun_Year <- dcast (Abbrev~Name, value.var = "Abun", data = Abun_Year)
  
  # Preparing dataset
  rownames(Abun_Year) <- Abun_Year$Abbrev
  Abun_Year <- Abun_Year[,2:length(Abun_Year)]
  
  # Traits organization data 
  TraitSPP0 <- TraitSPP
  TraitSPP0$Name <- rownames (TraitSPP0)
  rownames (TraitSPP0) <- NULL
  TraitSPP0$Name <- tolower(TraitSPP0$Name)
  TraitSPP0$Name <- gsub (" ", "_", TraitSPP0$Name)
  TraitSPP0 <- TraitSPP0[order(TraitSPP0$Name, decreasing = F), ]
  
  #setdiff(TraitSPP$Name, colnames(Abun_Year))
  
  rownames(TraitSPP0) <- TraitSPP0$Name
  TraitSPP0 <- TraitSPP0[,1:length(TraitSPP0)-1]
  
  # Estimating CWM
  Abun_Year0 <- data.matrix(Abun_Year)
  cwm_values <- functcomp(TraitSPP0,Abun_Year0, CWM.type = "all")
  
  cwm_results[[j]] <- cwm_values
  
  # notify (title = "CWM estimation",
  #         msg = c("Finish year", j))
  
}
for (i in 1:length(cwm_results)){
  cwm_results[[i]]$year <- 1967+i
  cwm_results[[i]]$Location <- rownames (cwm_results[[i]])
}
cwm_results   <- do.call("rbind", cwm_results)
rownames(cwm_results) <- NULL

# Estimating CWM separate per location -- The values are same
Abun_Location <- list ()
for (j in 1:length(AbunYear)) {
  
  Abun_Year <- AbunYear[[j]]
  
  # Removing Nocturnal species from Abundance estimation dataset
  Abun_Year <- melt(Abun_Year, id=c("Circle"))
  colnames (Abun_Year) <- c("Abbrev","Name","Abun")
  Abun_Year <- droplevels(Abun_Year[Abun_Year$Name %in% rownames(TraitSPP), ])
  Abun_Year$year <- 1967+j
  Abun_Location [[j]] <- Abun_Year
  
}
Abun_Location <- do.call("rbind",Abun_Location)
matrix_location <- split (Abun_Location, Abun_Location$Abbrev)
cwm_results_location <- list () 
for (m in 1:length(matrix_location)){
  
  Abun_Year <- matrix_location[[m]]
  Abun_Year$Id_LocYear <- paste (Abun_Year$Abbrev, Abun_Year$year, sep = "_")
  
  # Removing Nocturnal species from Abundance estimation dataset
  Abun_Year <- droplevels(Abun_Year[Abun_Year$Name %in% rownames(TraitSPP), ])
  
  Abun_Year$Name <- tolower(Abun_Year$Name) 
  Abun_Year$Name <- gsub (" ", "_", Abun_Year$Name)
  Abun_Year <- Abun_Year[order(Abun_Year$Name, decreasing = F), ]
  Abun_Year <- dcast (Id_LocYear~Name, value.var = "Abun", data = Abun_Year)
  
  # Preparing dataset
  rownames(Abun_Year) <- Abun_Year$Id_LocYear
  Abun_Year <- Abun_Year[,2:length(Abun_Year)]
  
  # Traits organization data 
  TraitSPP0 <- TraitSPP
  TraitSPP0$Name <- rownames (TraitSPP0)
  rownames (TraitSPP0) <- NULL
  TraitSPP0$Name <- tolower(TraitSPP0$Name)
  TraitSPP0$Name <- gsub (" ", "_", TraitSPP0$Name)
  TraitSPP0 <- TraitSPP0[order(TraitSPP0$Name, decreasing = F), ]
  
  #setdiff(TraitSPP$Name, colnames(Abun_Year))
  
  rownames(TraitSPP0) <- TraitSPP0$Name
  TraitSPP0 <- TraitSPP0[,1:length(TraitSPP0)-1]
  
  # Estimating CWM
  Abun_Year0 <- data.matrix(Abun_Year)
  cwm_values <- functcomp(TraitSPP0,Abun_Year0, CWM.type = "all")
  cwm_values$Abbrev <- rownames(cwm_values)
  rownames(cwm_values) <- NULL
  cwm_values <- separate(cwm_values, col =Abbrev, into = c("Abbrev","year"), sep = "_")
  
  
  cwm_results_location[[m]] <- cwm_values
  
  # notify (title = "CWM estimation",
  #         msg = c("Finish location", m))
  
}
cwm_results_location <- do.call ("rbind",cwm_results_location)
cwm_results_location <- cwm_results_location[order(cwm_results_location$Abbrev, decreasing = F),]

#saveRDS(cwm_results, "DataInter/cwm_results.rds")

# 5. MFA Analysis ------
cwm_order <- cwm_results
cwm_order$Id <- paste (cwm_order$Location, cwm_order$year, sep = "_")
rownames(cwm_order) <- cwm_order$Id
cwm_order <- cwm_order[,1:25]


Abun_LocationYear <- Abun_Location
Abun_LocationYear$Id <- paste(Abun_LocationYear$Abbrev, Abun_LocationYear$year, sep="_")
Abun_LocationYear <-  dcast(Id~Name, value.var = "Abun", data = Abun_LocationYear)
rownames(Abun_LocationYear) <- Abun_LocationYear$Id
Abun_LocationYear <- Abun_LocationYear[,2:length(Abun_LocationYear)]

cwm_full <- cbind (cwm_order,Abun_LocationYear)

# define columns in groups
morp  <- 7
mass  <- 1
diet  <- 10
fstra <- 7
n_spp <- ncol(Abun_LocationYear)

res <- MFA (cwm_full, group = c(morp, mass, diet, fstra, n_spp), 
            name.group = c("Morphologic","Mass","Diet","FStrategy", "Species"),
            num.group.sup = c(5), graph = F)
#saveRDS(cwm_full, "DataInter/cwm_full.rds")
#saveRDS(res, "DataInter/MFA_Results.rds")


# 6. Estimate species richness ----------
richness_db <- AbunYear
for (i in 1:length(richness_db)){
  richness_db[[i]]$year <- 1967+i
}
richness_db <- do.call("rbind", richness_db)
richness_db <- melt (richness_db, id=c("Circle","year"))
richness_db$value [richness_db$value<0.95] <- 0
richness_db <- droplevels (subset (richness_db, value>0))
richness_db$value [richness_db$value>=0.95] <- 1

tapply(richness_db$value, list (richness_db$Circle,richness_db$year), sum)

