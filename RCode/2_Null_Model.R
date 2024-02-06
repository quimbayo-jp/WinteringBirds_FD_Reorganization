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

library (picante)
library (gawdis)
library (reshape)
library (foreach)
library (plyr)
library (dplyr)
library (mFD)
library (tictoc)

# Null model -----------------------------------------
# 1. Calling data ------------------------------------
AbunYear  <-readRDS  ("DataRaw/commMats_byYear.rds")
TraitSPP  <-read.csv ("DataInter/Traits_Birds.csv", stringsAsFactors = T, 
                     dec = ".", sep = ",", header = T, row.names = 1)
TraitSPP$Mass <- log10(TraitSPP$Mass)

# IMPORTANT
# Overall 296 bird species was observed over time. 
# Total sites: 1183
# Total years:52 (1968-2019)

# Name comparison between TraitDataset and AbundDataset
setdiff(rownames(TraitSPP), colnames(AbunYear[[1]][,2:length(AbunYear[[1]])]))

# 2. Exclude Nocturnal and pelagic species from Trait database ------------
TraitSPP <- droplevels (subset(TraitSPP, Nocturnal==0)) # Ten species were removed from this dataset. 
TraitSPP <- droplevels (subset(TraitSPP, PelagicSpecialist==0)) 
TraitSPP <- TraitSPP[,1:28] # Nocturnal and Pelagic.Specialist columns were removed
TraitSPP <- TraitSPP[,c(1:5,9:length(TraitSPP))]

# 3. Estimating regional pool ------------------------------------------
domainsCoord <- read.csv("DataInter/RegionalPool_coord.csv", heade=T, sep =",")

localPool <- list ()

foreach (i = 1:length(AbunYear)) %do% {
  sppPool <- AbunYear[[i]]
  sppPool <- melt(sppPool, id=c("Circle")) 
  sppPool <- droplevels(subset(sppPool, value>0))
  sppPool$value <- ifelse(sppPool$value>0, 1, sppPool$value)
  
  localPool[[i]] <- sppPool
  
}


localPool <- do.call("rbind", localPool)
colnames (localPool) <- c("Abbrev","Name","Presence") 
localPool <- ddply (localPool,. (Abbrev, Name), summarise, Pres.Years=sum(Presence))
# Comment: This object "localPool" showed that some species only were observed in some years


RegionalPool <- left_join(domainsCoord, localPool)
RegionalPool$Presence <- 1

DomainPool <- ddply (RegionalPool,. (Abbrev,DomainName),
                     summarise, RegRichness=length(unique(Name)))

numberReps <- 999
null_fric  <- list()

tic("Start time")
for (m in 1:numberReps) {
  
# 4. Building Null Model ---------------------------------------------
null_findices <- list ()

  for (j in 1:length(AbunYear)) {
  # Database selected per each year
  Abun_Year <- AbunYear[[j]]
  
  # Removing Nocturnal species from Abundance estimation dataset
  Abun_Year <- melt(Abun_Year, id=c("Circle"))
  colnames (Abun_Year) <- c("Abbrev","Name","Abun")
  Abun_Year <- droplevels(Abun_Year[Abun_Year$Name %in% rownames(TraitSPP), ])
  Abun_Year <- droplevels (subset (Abun_Year, Abun>=1))
  
  # Joint abundance data and domain information 
  Abun_Year_Domain <- inner_join(Abun_Year, DomainPool, by=c("Abbrev"="Abbrev"))
  
  # Randomize via Independent swap null (holds row and column totals constant; richness and species occurrence frequencies across the gradient)
  Species_Domain  <- select(Abun_Year_Domain, Abbrev,DomainName,Name,Abun)
  Species_DomainU <- distinct(Species_Domain)
  domainIDs       <- as.factor(unique (Species_DomainU$DomainName))
  
  for (i in 1:99){
    # All these new matrix are only for only one year.
    # Thus this process should be repeat in all years
    
    comm_null_all <- list()
    
      for (k in 1:length(domainIDs)){
      comm_null_domain <- list()
      #isolate all data for a given domain
      comm_domain <- droplevels(subset(Species_Domain, DomainName == domainIDs[k])) 
      #make a sample matrix that domain can read (with species names and any numeric value to indicate presence)
      comm_domain_samp <- select(comm_domain, Abbrev, Abun, Name) #comm name, sp abundance, sp name
      #subset the data as needed to convert
      comm_domain_mat <- sample2matrix(comm_domain_samp) #convert to matrix format
      comm_domain_mat <- decostand(comm_domain_mat, method = "pa") #make binary (still the empirical matrix)
      comm_domain_mat_rdm <- randomizeMatrix(comm_domain_mat, null.model = "independentswap", iterations = 1000) #creates randomized occurrence matrix
      comm_domain_samp_rdm <- matrix2sample(comm_domain_mat_rdm)  #take the new randomized matrix to sample format
      colnames(comm_domain_samp_rdm) <- c("Abbre.null", "pres.null", "spName.null")
      
      #order the empirical and randomized data exactly the same and bind them together
      comm_domain_samp_rdm_order <-comm_domain_samp_rdm[order(comm_domain_samp_rdm$Abbre.null),]
      comm_domain_order <- comm_domain[order(comm_domain$Abbrev),]
      comm_domain_order <- droplevels(subset(comm_domain_order,Abun>0))
      comm_null_domain1 <- cbind(comm_domain_samp_rdm_order, comm_domain_order)
      #select and order column to match SpOcc
      comm_null_domain1 <- select(comm_null_domain1, Abbrev, DomainName, Name, Abun,spName.null)
      #join with species info
      comm_null_mo_k <- left_join(comm_null_domain1, Species_DomainU)
      
      comm_null_all <- rbind(comm_null_all, comm_null_mo_k)  #this step combines all the data for all the domain in the dataset = a complete random set of observations per bin per domain  It should have the same number of observations as the empiricial dataset (spocc).
      }
    
    #saveRDS(comm_null_all, file=paste0("Null_models_matrix/null_model_year_",j,".rds"))
    
    }
  
  # Calculate null functional richness ---- using the null assemblages
  null_abund <- cast(Abbrev~spName.null, value = "Abun",
                      data = comm_null_all, fill = 0)
  rownames(null_abund) <- null_abund$Abbrev
  null_abund <- null_abund[,2:length(null_abund)]
  
  null_traits <- TraitSPP
  setdiff(colnames(null_abund), rownames(null_traits))
  
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
  # Distance matrix
  traits_dist_null <- gawdis(null_traits, w.type = "optimized", groups = traitGroups, fuzzy = fuzzyGroups, opti.maxiter=3) #300 is the default
     
  # Quality functional space
  fspaces_quality_birds_null <- mFD::quality.fspaces(sp_dist = traits_dist_null,
                                                maxdim_pcoa = 10,
                                                deviation_weighting = "absolute",
                                                fdendro = "average")
  sp_faxes_coord_birds_null <- fspaces_quality_birds_null$"details_fspaces"$"sp_pc_coord"
  
  # Estimation Functional indices 
  fd_indices_birds_null   <- mFD::alpha.fd.multidim(sp_faxes_coord = sp_faxes_coord_birds_null[,c("PC1","PC2","PC3","PC4")],
                                                 asb_sp_w = data.matrix(null_abund),
                                                 scaling = TRUE,
                                                 check_input = TRUE,
                                                 details_returned = TRUE)
  
  null_findices[[j]] <- fd_indices_birds_null$functional_diversity_indices
  
  
  }
  # Add year in each list
  for (n in 1:length(null_findices)){
    null_findices[[n]]$year <- 1967+n
    null_findices[[n]]$Abbre <- rownames (null_findices[[n]])
  }
  
  #Merge all list that contain the functional indices estimated per year
  null_full_findices_db   <- do.call("rbind", null_findices)
  # Separated FRic per year
  FRIc_db <- cast (Abbre~year, value = "fric", data =null_full_findices_db)
  
  null_fric [[m]] <- FRIc_db
  
  # notify (title = "Null functional indices estimation",
  #         msg = c("Finish replicate", m))
  
}
toc()

# 5. Extracting Null.FRic in all years ----------------------------------------
null_fric_year <- null_fric

for (z in 1:length(null_fric_year)){
  null_fric_year[[z]]$null_sim <-paste0("nullSimu_",z)
}

# Combining all null values estimated of functional richness from null communities 
total_null_com <- do.call ("rbind",null_fric_year)

total_null_com <- as.data.frame(total_null_com)
total_null_com <- melt(total_null_com, id=c("Abbre","null_sim"))
colnames (total_null_com)[3]<- "year"
colnames (total_null_com)[4]<- "null_FRic"


year_matrix <- split (total_null_com, total_null_com$year)

null_values_FRic <- list ()
for (y in 1:length(year_matrix)){
  nullIndex <- year_matrix[[y]]
  nullIndex <- cast (Abbre+year~null_sim, value = "null_FRic", data = nullIndex)
  nullIndex$null_mean <- apply(nullIndex[,3:length(nullIndex)], 1, mean)
  
  null_values_FRic[[y]] <- nullIndex
}
null_values_FRic_mean <- do.call ("rbind",null_values_FRic)
null_values_FRic_mean <- as.data.frame(null_values_FRic_mean)

null_values_FRic_sd <- list ()
for (y in 1:length(year_matrix)){
  nullIndex <- year_matrix[[y]]
  nullIndex <- cast (Abbre+year~null_sim, value = "null_FRic", data = nullIndex)
  nullIndex$null_sd <- apply(nullIndex[,3:length(nullIndex)], 1, sd)
  
  null_values_FRic_sd[[y]] <- nullIndex
}
null_values_FRic_sd <- do.call ("rbind", null_values_FRic_sd)
null_values_FRic_sd <- as.data.frame(null_values_FRic_sd)

full_null_values_FRic <- left_join(null_values_FRic_mean,null_values_FRic_sd) 


