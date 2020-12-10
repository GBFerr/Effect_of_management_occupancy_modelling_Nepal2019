#########################################################################

# Effect of management on mammal occupancy in Terai Arc Landscape, Nepal
# Creating species x date occurrence matrix
# Data from 2019 camera trap survey

#########################################################################

{##### Part 1 - building Species presence/absence matrices #####
  
  
  # Detection/non-detection data for all spp photographed in the 2019 CT survey
  # Restricted to period between 15 March 2019 and 15 April 2019
  # Not all images obtained in this period were tagged 
  # we adopted 1-min interval between sequential images to select images to be tagged
  
  library(unmarked)
  library(here)
  library(dplyr)
  library(data.table)
  library(chron)
  
  #setwd("Z:/biome_health_project_files/country_files/nepal/our_papers/management_effect_2019CTsurvey")
  # server is extremely slow atm, will copy files locally
  
  setwd("C:/Users/Guilherme/Dropbox/!BHP_storage/Nepal2019_occuAnalysis")
  
  # load file with tagged images and metadata (date and time of the image)
  # must have cols w/ at least: species/label, CT site, date and time
  #data <- fread("Z:/biome_health_project_files/country_files/nepal/processed_data/2019_master_tags&meta.csv") # if loading from the server
  data <- fread("2019_master_tags&meta.csv")
  
  #checking consistency in labels used
  head(data)
  sort(unique(data$label)) # ok
  
  # renaming some cols in Data to match the functions below
  data <- rename(data, Site = camera_loc)
  data <- rename(data, Species = label)
  data <- rename(data, DateTime = datetime_true)
  str(data)
  data$DateTime <- as.Date(data$DateTime,"%d/%m/%Y") # making sure date col is Date
  
  # typo in Nilgai
  sort(unique(data$Species))
  data$Species <- sub(pattern="Nilagi", replace="Nilgai",data$Species) # correcting typo
  
  # loading matrix indicating whether camera trap site was working (1) or not (NA)
  all_cams <- read.csv("SurveyEffort_Mar-Apr_CTs_2019.csv") 
  head(all_cams)
  str(all_cams)
  all_cams$X <- as.Date(all_cams$X,"%d/%m/%Y") # making sure date is Date
  
  startDate <- as.Date("15/03/2019", "%d/%m/%Y") # 1st survey date
  endDate <- as.Date("15/04/2019", "%d/%m/%Y") # last survey date
  
  ####  calcOcc function - this function 'reads' the csv with labels/species, sites, and dates to create a 1/0 matrix
  calcOcc <-
    function(species, # species name - in dataframe - that the function is to be run for
             d = d, # dataframe with species, site, and each date it was seen at that site - must have a columns called Species, Site and DateTime
             all_cams = all_cams, # matrix with all the survey dates, 1s for dates when a camera was working/deployed and NAs for when it wasn't
             startDate = startDate,#start date in date format
             endDate = endDate) {
      # Make a vector of breaks
      brks <-seq(startDate, endDate, by = "day")   #makes a sequence of all the days from start to end
      
      # Create an empty matrix of dim sites x time periods
      occ <-matrix(0, ncol = length(unique(d$Site)), nrow = length(brks))
      colnames(occ) <- sort(unique(d$Site))
      rownames(occ) <- strftime(brks, format = "%Y-%m-%d")
      
      for (s in unique(d$Site)) {
        #this loops through each site and inserts 1's on days which there were captures
        seen <- NA
        captures <-na.omit(d$DateTime[d$Species == species & d$Site == s])
        # Were animals seen at the site
        seen <- which(brks %in% captures)
        # If the species was seen, occ = 1
        col_i <- which(colnames(occ) == s)
        occ[seen, col_i] <- 1
      }
      
      occ <- occ * all_cams[, 2:ncol(all_cams)]
      print(paste0(species, " done!"))
      species_name <- gsub(" ", "", species)
      row.names(occ) <- brks
      write.csv(occ, paste0("1d_matrix_", species_name, ".csv"))
      return(occ)
      
      
    }
  # applying function to each species (label)
  lapply(
    X = unique(data$Species),
    FUN = calcOcc,
    d = data,
    all_cams=all_cams,
    startDate = startDate,
    endDate = endDate) # this will save CSVs of spp matrices in the working directory
  
  # each csv generated represents the detection history of a species recorded 
  # 1: spp detected on that day and site (presence)
  # 0: spp not detected on that day and site (absence)
  # NA: camera trap was not functional on that day and site
  
  #### here ends creation of site x date matrix for each species ####
  
  # code below will take these site x data matrices and collapse them in longer time intervals
  # i.e. it will aggregate each data (1 survey day) into survey occasions of > 1 day 
  # below I used a 5-day survey occasion; so the output will be site x 5-day spp occurrence matrices
  
  ####  aggregating days into longer survey occasions ####
  ## timestepper - creates matrices of a given timestep (i.e. number of days aggregated), 
  # can choose to include or exclude NAs
  # but in reality must always exclude NAs otherwise it will consider days when CT wasn't working as actual survey days
  
  timestepper <- function(occ_in, timestep, na_mode = "include") {
    if (na_mode == "include") {
      occ_in[is.na(occ_in)] <-0   #replacing NAs with 0s if we want to include them in analysis.
    }
    
    if (timestep > nrow(occ_in) / 2) {
      print(paste(
        "Time step is too large! Please reduce to",
        nrow(occ_in) / 2 ,
        "or less."
      ))
    } else {
      start <- seq(1, nrow(occ_in), by = timestep)
      end <- seq(timestep, nrow(occ_in), by = timestep)
      
      if (length(start) > length(end)) {
        start <- start[-length(start)]
      }
      
      timesteps <- matrix(nrow = length(start), ncol = ncol(occ_in))
      colnames(timesteps) <- colnames(occ_in)
      rownames(timesteps) <-
        paste(rownames(occ_in)[start], rownames(occ_in)[end], sep = ":")
      
      for (i in 1:length(start)) {
        timestep_out <- colSums(occ_in[start[i]:end[i], ])
        timesteps[i, ] <- timestep_out
        timesteps[timesteps > 0] <- 1
      }
      
      timesteps <- t(timesteps)
      return(timesteps)
      
    }
    
  }
  
  # reading in all CSVs with spp presence absence; these are the original site x date (1 day) matrices
  filenames <- list.files("C:/Users/Guilherme/Dropbox/!BHP_storage/Nepal2019_occuAnalysis/matrices_out", pattern="*.csv", full.names=TRUE)
  ldf <- lapply(filenames, read.csv) # reading all CSVs into a list
  head(ldf)
  ldf[[1]] # just checking spp 1
  
  # creating object with spp names
  label <- basename(filenames)
  label <- sub(pattern=".csv", replace="",label)
  label <- sub(pattern="1d_matrix_", replace="",label)
  names(ldf) <- label
  
  # eliminating col "X" with dates before applying timestepper func
  new_list <- lapply(ldf, function(x) x%>% select(-X)) 
  new_list[[1]]
  
  #applying timestpper func to creat 5-day matrices for each spp
  matrix_5d <- lapply(X = new_list,
                      FUN = timestepper,
                      timestep = 5, # this is defining the survey occasions as 5 days; can be modified to have a longer or shorter occasion
                      na_mode = "exclude")  # must be set as exclude otherwise will consider days when CT wasn't working as actual survey days
  
  matrix_5d[[1]] # matric for spp 1; sites on rows and survey occasions on cols
  names(matrix_5d) <- label # adding species names
  
  # keeping only native mammals - it generated matrices for people, peacock, dog, etc
  matrix_native_5d <- matrix_5d[c(1,5,7,12,14,16,17:19,21:23,25:30,33,35,37:40,43)]
  
  # just checking
  names(matrix_native_5d) # must have only native mammal spp
  barkingdeer <- matrix_native_5d[[1]]
  chital <- matrix_native_5d[[2]]
  sum(barkingdeer, na.rm = T) # number of detections for barking deer considering the 5-day survey occasion (i.e. max 1 detection per 5 days)
  sum(chital, na.rm = T) # same as above for chital
  
## here ends the process of creating spp occurrence matrices for the detected species ##
  
  # code below will create all-0 matrices for hypothetical spp never detected during the survey (augmented species)
  # this is part of the data-augmentation approach to estimate spp richness via multi-species occupancy modelling (Bayesian approach)
  # if running single-species occupancy models there is no need to create augmented spp
  
  Aug01 <- matrix_native_5d[[1]] # getting matrix for spp1 as base
  Aug01[Aug01==1] <-0 # making all detection 0
  sum(Aug01, na.rm=T )
  Aug02 <- Aug01
  sum(Aug02, na.rm=T)
  Aug03 <- Aug01
  Aug04 <- Aug01
  Aug05 <- Aug01
  Aug06 <- Aug01
  Aug07 <- Aug01
  Aug08 <- Aug01
  Aug09 <- Aug01
  Aug10 <- Aug01
  Aug11 <- Aug01
  Aug12 <- Aug01
  Aug13 <- Aug01
  Aug14 <- Aug01
  Aug15 <- Aug01
  sum(Aug10, na.rm=T) # checking
  
  # adding augmented spp to the list of native spp matrices
  matrix_native_5d$Aug01 <- Aug01
  matrix_native_5d$Aug02 <- Aug02
  matrix_native_5d$Aug03 <- Aug03
  matrix_native_5d$Aug04 <- Aug04
  matrix_native_5d$Aug05 <- Aug05
  matrix_native_5d$Aug06 <- Aug06
  matrix_native_5d$Aug07 <- Aug07
  matrix_native_5d$Aug08 <- Aug08
  matrix_native_5d$Aug09 <- Aug09
  matrix_native_5d$Aug10 <- Aug10
  matrix_native_5d$Aug11 <- Aug11
  matrix_native_5d$Aug12 <- Aug12
  matrix_native_5d$Aug13 <- Aug13
  matrix_native_5d$Aug14 <- Aug14
  matrix_native_5d$Aug15 <- Aug15
  
  names(matrix_native_5d)
  sum(matrix_native_5d[[28]], na.rm=T) # sum =0; ok as it is an aug spp
  
  # transforming list of matrices into a 3D array: site x survey occasion x species
  library(abind)
  X <- abind(matrix_native_5d, along=3) # this will be used in as the spp input of multi-species occupancy model
  str(X)
  X[19:142,,28] # just checking the matrix for spp 28
  
  
}