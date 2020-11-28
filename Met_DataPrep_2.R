# Apply recovery correction to Met_Data output
# Load libraries and other functions from Metabolomics_01.R
# Run Met_DataPrep.R, then SpikeRecovery.R

# check current directory
here()

# transfer data files to project folder
# this is set up as a project so copy files into subfolders in the project folder
# check Standard curve info

# Expt specific info
# Enter information here
# source("NM-08/NM-08_info.R") # Separate R script for experiment specific information
# copy over standard curve filter info to line 40

Time <- Sys.time() # use to set a timestamp in filename

# load files

# data output from Met_DataPrep.R
suffix <- "_20201127_053709" # select specific timestamped file, otherwise set as NULL
if(deconv == "Y") {
  comp_datafile <- paste0(base_filename, "Deconv_compiled", suffix, ".csv") # output from Deconv_DataPrep.R
} else {
  comp_datafile <- paste0(base_filename, "compiled", suffix, ".csv") # Output from Met_DataPrep.R
}

data <- read_csv(paste(Exp,"data",comp_datafile, sep = "/"), na = "N/A")

# recovered spike info from SpikeRecovery.R
recov_suffix <- "_20201127_054047" # select specific timestamped file, otherwise set as NULL
recov_datafile <- paste0(base_filename, "spk_recovery", recov_suffix, ".csv") # output from SpikeRecovery.R
recovInfo <-  read_csv(paste(Exp,"results", recov_datafile, sep = "/"), na = "N/A")

# data preparation
data <- data %>% clean_names()
recovInfo <- recovInfo %>% clean_names()


# Comparison with standard control mix
# add "met_comp" when using deconv
expect_std <- data %>% 
  filter (sample_name2 == "std_high" | sample_name2 == "std_low") %>% 
  group_by_at (vars(all_of(c("metabolite", "spike")))) %>%
  summarise (exp_avg = mean(calc_conc), exp_err = se(calc_conc)) %>% 
  mutate (level = spike) %>% 
  mutate (act_conc = case_when(level == "H" ~ 5000,level == "L" ~ 500)) %>% 
  mutate ('exp-act' = round((exp_avg - act_conc),3))

recov_prop <- left_join(recovInfo, expect_std, by = c("metabolite","level")) %>% 
  mutate (recov_div_exp = (recov_avg/exp_avg)) %>% 
  mutate (recov_div_act = (recov_avg/act_conc))

filename <- paste0(base_filename, "recovery", format(Time, "_%Y%m%d_%H%M%S"), ".csv") 
write.csv (recov_prop, file = paste(Exp,"results",filename, sep = "/"), na = "")

# Applying correction factor
# add "met_comp" when using deconv
recov_prop_select <- recov_prop %>% # not entirely whether to use the data from H, L or both
  filter (level == "H") %>% 
  select (all_of(c("metabolite","recov_div_exp","recov_div_act")))

corrected_data <- left_join(data,recov_prop_select, by = c("metabolite")) %>% 
  mutate (corr_conc = case_when(category == "Sample" ~ calc_conc/recov_div_exp)) %>% 
  mutate (corr_pmol_mg = corr_conc*volume/10^3/wet_weight)

filename <- paste0(base_filename, "compiled_corr", format(Time, "_%Y%m%d_%H%M%S"), ".csv") 
write.csv(corrected_data, file = paste(Exp,"data",filename, sep = "/"), na = "")