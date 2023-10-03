## --------------------------------------------------------------------------------------------------------- ##
##  Project Name:  Multiple axes of avian functional diversity in North America                              ##
##                                                                                                           ##
##  Objective:     Here, we examined variation avian fauna in North America over time.                       ##
##                                                                                                           ##
##  Authors:       Juan P. Quimbayo                                                                          ##
##                 Marta Jarzyna                                                                             ##
##                                                                                                           ##
##  Date:          2022-09-21                                                                                ##
##                                                                                                           ##
##  Notes:         1. This file is intended to provide a guide to the basic                                  ##
##                    project workflow, attempting to 'model analysis' the steps                             ##
##                    necessary to conduct the analyses & visual outputs.                                    ##
##                                                                                                           ##
## --------------------------------------------------------------------------------------------------------- ##
# clean up
rm(list=ls())

library (dplyr)
library (viridis)
library (plyr)
library (tidyverse)
library (mgcv)    
library (sjPlot)
library (brms)
library (car)
library (ggmcmc)
library (tidybayes)
library (ggdist)
library (ggpubr)
library (betareg)
library (factoextra)

# 1. Functional data ------------------------------------------------------

 
DB_FIndex <- readRDS ("DataInter/DB_FIndex.rds")
MFA_results <- readRDS ("DataInter/MFA_Results.rds")
MFA_results <- data.frame (MFA_results$ind$coord)
MFA_results$ID_unique <- rownames(MFA_results)
rownames(MFA_results) <- NULL


# 2. Null model data ------------------------------------------------------

nullFRic <- readRDS ("DataInter/nullFRic.rds")
nullFRic$ID_unique <- paste(nullFRic$Abbre,nullFRic$year, sep = "_")
nullFRic <- nullFRic [,c(5,3:4)]



# 3. Joint functional data ------------------------------------------------


DB_FIndex <- left_join(DB_FIndex, nullFRic)
DB_FIndex <- left_join(DB_FIndex, MFA_results)
DB_FIndex$ses_fric <- (DB_FIndex$fric - DB_FIndex$meanNull)/DB_FIndex$SDNull


# 4. Call environmental and elevation data --------------------------------

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


# 5. Estimating slope from environmental and elevation data ---------------

# Extracting slope values from environmental factors 
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
                               coef_ppt= do.call("rbind",slope_ppt))
slopes_env_factors <- merge (slopes_env_factors, elevation[,c("Abbrev","Elevation")],
                             by="Abbrev", all.x = T)

rm (slope_meanTemp,slope_maxTemp,slope_minTemp,slope_ppt)



# 6. Extracting slope values from functional indices ----------------------

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

slopes_f_indices <- data.frame(Abbrev=unique(DB_FIndex$Location),
                               coef_sric= do.call("rbind",slope_values_sric),
                               coef_cfric= do.call("rbind",slope_values_cfric),
                               coef_feve= do.call("rbind",slope_values_feve),
                               coef_fdiv= do.call("rbind",slope_values_fdiv),
                               coef_fori= do.call("rbind",slope_values_fori),
                               coef_Dim1= do.call("rbind",slope_values_Dim1),
                               coef_Dim2= do.call("rbind",slope_values_Dim2)) 


rm (slope_values_sric,
    slope_values_cfric,
    slope_values_feve,
    slope_values_fdiv,
    slope_values_fori,
    slope_values_Dim1,
    slope_values_Dim2)

# 7. Evaluating correlation among environmental predictors ------------------------------

# Add Histograms in function Pairs
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

# Add Coeficients of correlation
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


# 8. Re-scaling environmental data ----------------------------------------

rescale_variables <- function (x) {(x - mean(x))/(1*sd (x))}
slopes_env_factors$coef_meanTemp_scale <- rescale_variables(slopes_env_factors$coef_meanTemp)
slopes_env_factors$coef_maxTemp_scale <- rescale_variables(slopes_env_factors$coef_maxTemp)
slopes_env_factors$coef_minTemp_scale <- rescale_variables(slopes_env_factors$coef_minTemp)
slopes_env_factors$coef_ppt_scale <- rescale_variables(slopes_env_factors$coef_ppt)
slopes_env_factors$elevation_scale <- rescale_variables(slopes_env_factors$Elevation)
dim(slopes_env_factors)

#pdf("Output/Figures/Fig_S1.pdf", height = 15, width = 15, pointsize=30)
pairs (slopes_env_factors[,7:length(slopes_env_factors)], pch=13,
       diag.panel=panel.hist, lower.panel=panel.cor,
       panel=panel.smooth)
#dev.off()


# 9. Calculating Principal component analysis  ----------------------------

tempPCA <- princomp (slopes_env_factors[,2:4], cor = F)
summary (tempPCA)
tempPCA$rotation

biplot (tempPCA)

res.pca <- prcomp(slopes_env_factors[,2:4], scale = TRUE)
summary (res.pca)
res.pca$rotation

slopes_env_factors <- cbind (slopes_env_factors, res.pca$x[,1:2])
slopes_env_factors2 <- cbind (slopes_env_factors, tempPCA$scores[,1:2])


# 10. Joint all Functional indices and predictors -------------------------

FD_Index_EnvFactors <- left_join(slopes_f_indices, slopes_env_factors)
FD_Index_EnvFactors <- left_join(slopes_f_indices, slopes_env_factors2)


# 11. Joint data with biogeography information ------------------------------

db_domains <- read.csv("DataInter/RegionalPool_coord.csv", header = T)

FD_Envir_Domain <- left_join(FD_Index_EnvFactors,db_domains)
FD_Envir_Domain <- na.omit(FD_Envir_Domain)


# 12. Bayesian models -----------------------------------------------------


## Using Bayesian approach
### Testing Variation inflate factor (VIF)
vif (lm(coef_cfric~PC1+
          PC2+
          coef_ppt_scale+
          elevation_scale, 
        data=FD_Envir_Domain))


## 12.1. Species richness model --------------------------------------------
  mRich  <- brm (bf(coef_sric~PC1+
                    PC2+
                    coef_ppt_scale+
                    elevation_scale +
                    (1|DomainName)),
               prior = c(prior(normal(0, 10), 'b'),
                         prior(normal(0, 50), 'Intercept'),
                         prior(student_t(3, 0, 2.5), 'sd'),
                         prior(student_t(3, 0, 2.5), 'sigma')),
               sample_prior = TRUE,
               data=FD_Envir_Domain,
               family = gaussian(),
               chains = 4, iter = 20000, warmup = 18000, cores = 4,
               control = list(adapt_delta=0.98), seed = 123, thin = 1)


## 12.2. Functional richness model --------------------------------------------
  mcFRic  <- brm (bf(coef_cfric~PC1+
                    PC2+
                           coef_ppt_scale+
                    elevation_scale +
                           (1|DomainName)),
                    prior = c(prior(normal(0, 10), 'b'), 
                              prior(normal(0, 50), 'Intercept'), 
                              prior(student_t(3, 0, 2.5), 'sd'), 
                              prior(student_t(3, 0, 2.5), 'sigma')),
                    sample_prior = TRUE,
                    data=FD_Envir_Domain,
                    family = gaussian(),
                    chains = 4, iter = 20000, warmup = 18000, cores = 4,
                    control = list(adapt_delta=0.98), seed = 123, thin = 1)


## 12.3. Functional evenness model --------------------------------------------
  mFEve <- brm (bf(coef_feve~PC1+
                   PC2+
                           coef_ppt_scale+
                           elevation_scale +
                           (1|DomainName)),
                      prior = c(prior(normal(0, 10), 'b'), 
                                prior(normal(0, 50), 'Intercept'), 
                                prior(student_t(3, 0, 2.5), 'sd'), 
                                prior(student_t(3, 0, 2.5), 'sigma')),
                      sample_prior = TRUE,
                      data=FD_Envir_Domain,
                      family = gaussian(),
                      chains = 4, iter = 20000, warmup = 18000, cores = 4,
                      control = list(adapt_delta=0.98), seed = 123, thin = 1)

## 12.4. Functional divergence model --------------------------------------------
  mFDiv <- brm (bf(coef_fdiv~PC1+
                   PC2+
                   coef_ppt_scale+
                   elevation_scale +
                           (1|DomainName)),
                      prior = c(prior(normal(0, 10), 'b'), 
                                prior(normal(0, 50), 'Intercept'), 
                                prior(student_t(3, 0, 2.5), 'sd'), 
                                prior(student_t(3, 0, 2.5), 'sigma')),
                      sample_prior = TRUE,
                      data=FD_Envir_Domain,
                      family = gaussian(),
                      chains = 4, iter = 20000, warmup = 18000, cores = 4,
                      control = list(adapt_delta=0.98), seed = 123, thin = 1)

## 12.5. Functional originality model --------------------------------------------
  mFOri <- brm (bf(coef_fori~PC1+
                   PC2+
                   coef_ppt_scale+
                   elevation_scale +
                           (1|DomainName)),
                      prior = c(prior(normal(0, 10), 'b'), 
                                prior(normal(0, 50), 'Intercept'), 
                                prior(student_t(3, 0, 2.5), 'sd'), 
                                prior(student_t(3, 0, 2.5), 'sigma')),
                      sample_prior = TRUE,
                      data=FD_Envir_Domain,
                      family = gaussian(),
                      chains = 4, iter = 20000, warmup = 18000, cores = 4,
                      control = list(adapt_delta=0.98), seed = 123, thin = 1)

## 12.6. Trait Dimension 1 model --------------------------------------------
  mDim1  <- brm (bf(coef_Dim1~PC1+
                    PC2+
                    coef_ppt_scale+
                    Elevation +
                    (1|DomainName)),
               prior = c(prior(normal(0, 10), 'b'), 
                         prior(normal(0, 50), 'Intercept'), 
                         prior(student_t(3, 0, 2.5), 'sd'), 
                         prior(student_t(3, 0, 2.5), 'sigma')),
               sample_prior = TRUE,
               data=FD_Envir_Domain,
               family = gaussian(),
               chains = 4, iter = 20000, warmup = 18000, cores = 4,
               control = list(adapt_delta=0.98), seed = 123, thin = 1)


## 12.7. Trait Dimension 1 model --------------------------------------------

  mDim2  <- brm (bf(coef_Dim2~PC1+
                    PC2+
                    coef_ppt_scale+
                    Elevation +
                    (1|DomainName)),
               prior = c(prior(normal(0, 10), 'b'), 
                         prior(normal(0, 50), 'Intercept'), 
                         prior(student_t(3, 0, 2.5), 'sd'), 
                         prior(student_t(3, 0, 2.5), 'sigma')),
               sample_prior = TRUE,
               data=FD_Envir_Domain,
               family = gaussian(),
               chains = 4, iter = 20000, warmup = 18000, cores = 4,
               control = list(adapt_delta=0.98), seed = 123, thin = 1)
