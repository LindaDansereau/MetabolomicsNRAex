# Convert peak area to concentrations for graphing and other analyses
# Load libraries and other functions from Metabolomics_01.R
# TODO: remove extra comment columns from standard curve table

# check current directory
here()

# transfer data files to project folder
# this is set up as a project so copy files into subfolders in the project folder
# check Standard curve info

# Expt specific info
# Enter information here
source("NM-08/NM-08_info.R") # Separate R script for experiment specific information
# copy over standard curve filter info to line 40

Time <- Sys.time() # use to set a timestamp in filename

# load files
data <- read_csv(paste(Exp,"data",datafile, sep = "/"), na = "N/A") # MQ data
StdData <- read_csv(paste(Exp,"data",stddatafile, sep = "/"), na = "N/A") # standard curve info
SampleData <- read_csv(paste(Exp,"data",sampleinfofile, sep = "/"), na = "N/A") # sample info

# Data preparation

# fix column headings
data <- data %>%
  clean_names()
StdData <- StdData %>%
  clean_names()
SampleData <- SampleData %>%
  clean_names()

# additional changes before merging data files
#keep_cols <- c("sample_name", "component_name", "expected_rt", "integration_type", "area", "retention_time")	
keep_cols <- c("sample_name", "component_name", "area", "retention_time")	
data <- data %>%
  select(all_of(keep_cols))

StdData <- StdData %>%
  filter (standard_curve == "NADmetA" | standard_curve == "NADmetB" | standard_curve == "NADmetC" | standard_curve == "NR_NAR")
#filter (standard_curve == "NADP" | standard_curve == "NADPH" | standard_curve == "NADmet")
#filter (standard_curve == "NAMN_standards" | standard_curve == "met_standards")

# merge data files
SampleInfoData <- inner_join(data, SampleData, by = "sample_name")
AllData <- inner_join(SampleInfoData, StdData, by = "component_name")

#calculations and other adjustments
data_2 <- AllData %>%
  mutate (id_number = factor(id_number)) %>%
  mutate (area = as.numeric(area)) %>%
  mutate (volume = as.numeric(volume)) %>%
  mutate (wet_weight = as.numeric(wet_weight)) %>%
  mutate (area_mg = area*volume/10^3/wet_weight) %>% 
  mutate (calc_conc = conc_c(area,slope,intercept)) %>%
  mutate (pmol_mg = calc_conc*volume/10^3/wet_weight)

# save file

filename <- paste0(base_filename, "compiled", format(Time, "_%Y%m%d_%H%M%S"), ".csv") 
write.csv(data_2, file = paste(Exp,"data",filename, sep = "/"), na = "")