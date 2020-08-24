#######################################################

### Exploring site covs for Nepal camera trap sites ###

#######################################################

### obtaining Land cover at each camera trap buffer ###
#adapted from: 
#https://gis.stackexchange.com/questions/249878/extract-raster-pixel-values-from-shapefile-in-qgis

setwd("Z:/biome_health_project_files/country_files/nepal/our_papers/management_effect_2019CTsurvey")
library(raster)
library(rgdal)
#### 500m buffer ####
#reading raster with Land Cover in Nepal in 2010 (6 categories of cover; Uddin et al 2010) and 500m buffer around CT sites
# raster from: https://rds.icimod.org/Home/DataDetail?metadataId=9224&searchlist=True
raster <- stack("C:/Users/Guilherme/Dropbox/Geoprocessamento/Base de dados vetoriais/Nepal/LandCoverNepal_2010/data/np_lc_2010_v2f.tif")
#### legend for raster layer
# 1-Forest; 2-Shrubland; 3-Grassland; 4-Agriculture; 5-Barren land/dry riverbed; 6-water

poly <- shapefile("Z:/biome_health_project_files/country_files/nepal/GIS_files/2019_camera_forest_buffers/500m_CT_buffer.shp")
dim(poly)
val <- extract(raster,poly) # result is a list
summary(val)
# just checking values, each element of the list is a buffer around a camera trap site
head(val[[1]],10)
head(val[[38]],10)
head(val[[138]],10)
tail(val[[138]],10)
tail(val[[1]],25)
tail(val[[60]],25)
names(poly)
# using a field from my vector to create an identifier (use unique values)
listnames <- poly$Real_name

# create a empty list to save data frames to export
valList <- list()

# create as many data frames as features used to extract
for(i in 1:length(val)){
  valList[[i]] <- data.frame(ID=listnames[[i]],Value = val[[i]][,1])
}

# join all data frames and save it to an csv file
write.csv(do.call(rbind.data.frame,valList),"LandCover_test.csv")
#read csv with pixel value (land cover) within each CT buffer
LC_data <- read.csv("LandCover_test.csv")
head(LC_data)
LC_CTbuffer <- as.table(table(LC_data$ID,LC_data$Value)) #number of pixels in each land cover category within each CT buffer
write.table(LC_CTbuffer,"LandCover_CT500mbuffer.csv", sep=",") 
#this csv was compared with layers in QGIS and matched perfectly for all buffers assessed (n=28)

#### 1000m buffer ####
# similar things for 1000m buffer
poly2 <- shapefile("Z:/biome_health_project_files/country_files/nepal/GIS_files/2019_camera_forest_buffers/1000m_CT_buffer.shp")
val2 <- extract(raster,poly2) # same raster

summary(val2)
# checking values, each element of the list is a buffer around a camera trap site

names(poly2)
# I will use a field from my vector to create an idintifier (use unique values)
listnames2 <- poly2$Real_name

# create a empty list to save data frames to export
valList2 <- list()

# create as many data frames as features used to extract
for(i in 1:length(val2)){
  valList2[[i]] <- data.frame(ID=listnames2[[i]],Value = val2[[i]][,1])
}

# join all data frames and save it to an csv file
write.csv(do.call(rbind.data.frame,valList2),"LandCover_1000mBuffer.csv")
#read csv with pixel value (land cover) within each CT buffer
LC_data1000 <- read.csv("LandCover_1000mBuffer.csv")
head(LC_data1000)
LC_CT1000buffer <- as.table(table(LC_data1000$ID,LC_data1000$Value)) #number of pixels in each land cover category within each CT buffer
write.table(LC_CT1000buffer,"LandCover_CT1000mbuffer.csv", sep=",") 
#this csv was compared with layers in QGIS and matched perfectly for all buffers assessed (n=14)

#### checking correlation between buffers of distinct sizes ####
# reading data again to change colnames and check correlation between 500m and 1k buffer
buffer1000 <- read.csv("LandCover_CT1000mbuffer.csv")
head(buffer1000)
colnames(buffer1000) <-  c("CT_site", "Forest", "Shrubland", 
                         "Grassland", "Agriculture", "Barren", "Water", "passed_visual_check")

buffer500 <- read.csv("LandCover_CT500mbuffer.csv")
head(buffer500)
colnames(buffer500) <-  c("CT_site", "Forest", "Shrubland", 
                           "Grassland", "Agriculture", "Barren", "Water", "passed_visual_check")

Forest <- cor(buffer500$Forest,buffer1000$Forest)
Shrubland <- cor(buffer500$Shrubland,buffer1000$Shrubland)
Grassland <- cor(buffer500$Grassland,buffer1000$Grassland)
Agriculture <- cor(buffer500$Agriculture,buffer1000$Agriculture)
Barren <- cor(buffer500$Barren,buffer1000$Barren)
Water <- cor(buffer500$Water,buffer1000$Water)
cor <- cbind(Forest,Shrubland, Grassland, Agriculture, Barren, Water)
head(cor)
write.csv(cor, "LandCover_correlation_500vs1000mBuffer.csv")

# lowest cor was 0.76 for shrubland; >0.95 cor for  the 2 most common land covers (forest & agric)
# thus buffer size unlikely to influence occupancy modelling
# will use 500m buffer for initial analysis

#### Calculating prop of land cover classes ####
head(buffer500)
numeric_buffer500 <- buffer500[,2:7] # object with numeric cols from buffer 500
buffer500$sumpix <- rowSums(numeric_buffer500) #getting total number of pixels (needed to get prop in the next step)
head(buffer500)
buffer500$propForest <- buffer500$Forest/buffer500$sumpix   #proportion of forest
buffer500$propAgric <- buffer500$Agric/buffer500$sumpix   #prop Agriculture

# prop natural cover, anything that's not agriculture; thus includes water and barren
# barren was included because majority (if not all) barren land in Bardiya landscape is dry riverbed
buffer500$propNatural <- 1- buffer500$propAgric 

# prop natural vegetation, sum of forest, shrubland and grassland
buffer500$propVegNat <- (buffer500$Forest+buffer500$Shrubland+buffer500$Grassland)/buffer500$sumpix
head(buffer500,15)
write.csv(buffer500, "LandCover_CT500mbuffer.csv")

buffer500 <- read.csv("LandCover_CT500mbuffer.csv") # reading again after a crash
#### adding Land Cover from Uddin 2010 to site variable spreadsheet ####
library(readxl)
covs <- read_excel("siteVariable_148sites_Nepal2019.xlsx")
head(buffer500)
prop <- buffer500[,c(2,11:14)]
head(prop)

covs <- merge(covs,prop, by="CT_site")
head(covs)

#### exploring site variables ####
# correlation
covs_numeric <- covs[,c(5,12:15,20:23)]
which(is.na(covs_numeric))
library(psych)
pairs.panels(covs_numeric, scale=FALSE)

#influence of management
covs$Management <- factor(covs$Management , levels=c("NP", "BZ", "OBZ"))
par(mfrow=c(2,2))
boxplot(covs$DistVillage ~ covs$Management, main= "village")
boxplot(covs$DistRoad ~ covs$Management, main= "road")
boxplot(covs$DistRiver ~ covs$Management, main= "river")
boxplot(covs$propVegNat ~ covs$Management, main= "vegetation")

# PCA
# variables used in PCA (Dist_city Dist_road Dist_water Dist_humset HDens_2km NDVImean_500m cattle_freq)
names(covs)
management <- covs[,4]
covsModel <- covs[,c("DistVillage","DistRoads","DistRiver","propVegNat")]

#### ggplot PCA
library(ggfortify)
pca_res <- prcomp(covsModel, center=TRUE, scale. = TRUE)
summary(pca_res)
p <- autoplot(pca_res, data=covs, colour="Management", size=3, alpha=1/1.5,
         loadings=TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)
p+theme_bw()+
  theme(legend.title=element_blank(), legend.position=c(0.15,0.85))
