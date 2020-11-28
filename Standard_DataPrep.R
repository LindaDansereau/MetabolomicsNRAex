# Prepare peak area from Standard curve runs for regression analysis
# Currently set for copy-paste into Prism
# TODO: Analysis within R 
# Load libraries and other functions from Metabolomics_01.R


# check current directory
here()

# transfer data files to project folder
# this is set up as a project so copy files into subfolders in the project folder

# Expt specific info
# Enter information here
source("NM-08/NM-08_info.R") # Separate R script for experiment specific information

Time <- Sys.time() # use to set a timestamp in filename

# load files
# data <- read_csv(paste(Exp,"data", datafile, sep = "/"), na = "N/A") # MQ data
StandardData <- read_csv(paste(Exp,"data",stdinfofile, sep = "/"), na = "N/A") # sample info
StandardConc <- read_csv(paste(Exp,"data",stdconcfile, sep = "/"), na = "N/A") # concentrations

# preparing deconvolved data
deConvData <- read_csv(paste(Exp,"data",deconvdatafile, sep = "/"), na = "N/A") # Data from Matlab output
data <- deConvData %>% 
  unite(met_comp, c(metabolite, component_name), remove = FALSE, na.rm = TRUE) %>% 
  pivot_longer(- c(met_comp, metabolite, component_name), names_to = "sample_name", values_to = "Area")


# Data preparation

# fix column headings
data <- data %>%  clean_names()
StandardData <- StandardData %>%  clean_names()
StandardConc <- StandardConc %>%  clean_names()

# additional changes before merging data files
#keep_cols <- c("sample_name", "component_name", "expected_rt", "integration_type", "area", "retention_time")	
keep_cols <- c("sample_name", "component_name", "area", "retention_time")	
data <- data %>%
  # select(all_of(keep_cols)) %>% 
  # filter(sample_name != "NH4OAc_blank")
  filter(str_detect(sample_name,  "NADmetC"))

# merge data files
AllData <- inner_join(data, StandardData, by = "sample_name")

#calculations and other adjustments
data_2 <- AllData %>%
  mutate (id_number = factor(id_number)) %>%
  mutate (area = as.numeric(area)) %>%
  mutate (metabolite = component_name)

# save file
filename <- paste0(base_filename, "NADmetC_standard_compiled", format(Time, "_%Y%m%d_%H%M%S"), ".csv") 
write.csv(data_2, file = paste(Exp,"data",filename, sep = "/"), na = "")

# #export for Prism - all metabolites
keep_cols2 <- c("sample_name2", "area", "id_number", "metabolite", "met_comp")	
data_print <- data_2 %>%
  select(all_of(keep_cols2)) %>% 
  spread(met_comp, area)

# rearrange for easier copy-paste, sort as needed
data_print <- data_print %>%
  pivot_longer(-c(id_number,sample_name2, metabolite)) %>%
  pivot_wider(names_from = id_number, values_from = value)

data_print <- left_join(data_print, StandardConc,by = "sample_name2")

filename <- paste0(base_filename, "NADmetC_standard_results", format(Time, "_%Y%m%d_%H%M%S"),".csv") 
write.csv(file = paste(Exp,"results", filename, sep = "/"), data_print)

