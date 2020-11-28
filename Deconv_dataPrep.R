# Expt specific info
# Enter information here
source("NM-08/NM-08_info.R") # Separate R script for experiment specific information

Time <- Sys.time() # use to set a timestamp in filename

# load files
deConvData <- read_csv(paste(Exp,"data",deconvdatafile, sep = "/"), na = "N/A") # Data from Matlab output
SampleData <- read_csv(paste(Exp,"data",sampleinfofile, sep = "/"), na = "N/A") # sample info
StdData <- read_csv(paste(Exp,"data",stddatafile, sep = "/"), na = "N/A") # standard curve info


# Data preparation

# Matlab output

# rename first two columns as metabolite and component_name
# NAD_NAAD_test <- NAD_NAAD_test %>% 
#   rename(metabolite = X1) %>% 
#   rename(component_name = X2)

# merge first two columns and pivot
data <- deConvData %>% 
  unite(met_comp, c(metabolite, component_name), remove = FALSE, na.rm = TRUE) %>% 
  pivot_longer(- c(met_comp, metabolite, component_name), names_to = "fileName", values_to = "Area")


# replace file names to match sample_name format
data <- data %>% 
  add_column(sample_name = str_replace_all(data$fileName, ".mzML", ""), .after = "met_comp") 

data <- data %>% 
  mutate(sample_name = str_remove_all(sample_name, "[']")) %>% 
  mutate(fileName = NULL)


# fix column headings
data <- data %>%
  clean_names()
SampleData <- SampleData %>%
  clean_names()
StdData <- StdData %>%
  clean_names()

# additional changes before merging data files
StdData <- StdData %>%
  filter (standard_curve == "NADmetA" | standard_curve == "NADmetB" | standard_curve == "NADmetC")

# merge data files
SampleInfoData <- inner_join(data, SampleData, by = "sample_name")
AllData <- inner_join(SampleInfoData, StdData, by = "component_name") %>% 
  mutate(metabolite.y = NULL) %>% 
  rename(metabolite = metabolite.x)

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
filename <- paste0(base_filename, "Deconv_compiled", format(Time, "_%Y%m%d_%H%M%S"), ".csv") 
write.csv(data_2, file = paste(Exp,"data",filename, sep = "/"), na = "")
