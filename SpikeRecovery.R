# data preparation for checking level of spike recovery in metabolome experiments
# graph output
# Load libraries and other functions from Metabolomics_01.R
# TODO: 

# start with compiled data from Met_DataPrep.R

# Expt specific info
# Enter information here
source("NM-08/NM-08_info.R") # Separate R script for experiment specific information

# Expt specific info
suffix <- "_20201127_053709" # select specific timestamped file, otherwise set as NULL
Note <- "NR_NaR only, Curve regression: Area, linear through 0, weighting = none" # comment for PDF

if(deconv == "Y") {
  comp_datafile <- paste0(base_filename, "Deconv_compiled", suffix, ".csv") # output from Deconv_DataPrep.R
} else {
  comp_datafile <- paste0(base_filename, "compiled", suffix, ".csv") # Output from Met_DataPrep.R
}
descrip <- "NR Aex Recovered Spike"


# PDF output info, change as needed
Graph.title <- "Recovered Spike Amounts"
Expt.title <- paste(Exp, descrip)
File.name <- paste(Exp, "_", descrip,".pdf",sep = "")
File.description <- "Metabolome analysis"
Today <- Sys.Date() # Today's date format:yyyy-mm-dd
Time <- Sys.time() 

# Filter and grouping lists, change as needed
# add "met_comp" to grouping spike for deconv results
# ids_interest <- c("52","1816","1819","258288","258290") # untested
grouping_spike <- c("id_number", "component_name", "metabolite", "slope", "spike", "intercept")
X_label <- "Spike Level"
Y_label <- "Concentration (nM)"
# reorder_list <- c("L","H")

# load data file
data <- read_csv(paste(Exp,"data",comp_datafile, sep = "/"), na = "N/A")
MetLists <- read_csv(paste(Exp,"data",met_datafile, sep = "/"), na = "N/A")


# additional changes to data
data <- data %>% 
  mutate (id_number = factor(id_number)) #%>% # id_number not always in expected format
  #mutate (replicate = factor(replicate))

MetLists <- MetLists %>% clean_names()

# keeping metabolites based on positive vs negative values

# NAD Metabolites to analyse
usedMet <- MetLists %>% filter(used_met == "Y") 
usedMet <- usedMet[['component_name']]

# Metabolites separated out due to different amounts in samples
# Lists change depending on experiment conditions (eg: tissue type)
# Check results with usedMet and then adjust lists - enter data into met_datafile and reload

topMet <-  MetLists %>% filter(used_met_level == "top") 
topMet <- topMet[['component_name']]

highMet <-  MetLists %>% filter(used_met_level == "high") 
highMet <- highMet[['component_name']]

medMet <-  MetLists %>% filter(used_met_level == "med") 
medMet <- medMet[['component_name']]

lowMet <-  MetLists %>% filter(used_met_level == "low") 
lowMet <- lowMet[['component_name']]

vlowMet <-  MetLists %>% filter(used_met_level == "vlow") 
vlowMet <- vlowMet[['component_name']]

keepMet <- MetLists %>% filter(keep_spike == "Y") 
keepMet <- keepMet[['component_name']]

# Determine recovered spike amount
Recov_data <- data %>% # average over individual samples
  filter (category == "Recov_cont") %>%
  mutate (spike = paste0("spike",spike)) %>% # change spike labels to make next steps easier
  group_by_at(vars(all_of(grouping_spike))) %>%
  summarise (RA = mean(area)) %>% # calculate recovery area - no actual effect of calculating mean when only one injection
  pivot_wider (names_from = spike, values_from = RA) %>% # spread spike info into columns
  mutate (across(c(starts_with("spike")), list( A = ~ . - spikeN))) %>% # subtract no spike area from each spike
  select (-spikeN_A) %>%  # remove no spike-no spike column
  mutate (across(c(contains("_A")), list( conc = ~ conc_c(. ,slope,intercept)))) %>% # determine concentration 
  pivot_longer(cols = contains("_conc"), names_to = "Level", values_to = "Conc") %>% # return data to long format
  select (- contains("spike")) %>% # remove unneeded columns
  mutate (Level = str_extract(Level, "[A-Z]{1,6}")) # return Level column to starting spike column status
  
# export data
# add "met_comp" to grouping spike for deconv results
Printdata <- Recov_data %>% 
  filter (!is.na(Conc)) %>% # remove any invalid values
  group_by_at(vars(all_of(c("metabolite","Level")))) %>% 
  summarise (Recov_Avg = mean(Conc), Recov_Err = se(Conc))

filename <- paste0(base_filename, "spk_recovery", format(Time, "_%Y%m%d_%H%M%S"), ".csv")
write.csv(file = paste(Exp,"results",filename, sep = "/"), Printdata)

# Reorder data if needed, create Avg and Name column for graph scripts
if(exists("reorder_list")) {
  Graphdata <- DataReorder(Recov_data, Level, reorder_list) %>%
    group_by_at(vars(all_of(c("id_number","Name","metabolite","component_name","Level")))) %>% 
    summarise (Avg = mean(Conc), Err = se(Conc))
  } else {
     Graphdata <- Recov_data %>% 
       mutate (Name = Level) %>% 
       group_by_at(vars(all_of(c("id_number","Name","metabolite","component_name","Level")))) %>% 
       summarise (Avg = mean(Conc), Err = se(Conc))
  }


# Colour selections for graphs
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
myPalette = getPalette(8) #not used?
colourCount = length(unique(Recov_data$id_number))
myPalette = getPalette(colourCount)

# Dot plots with spike levels, filter for Met levels
# Cycling through Metabolite levels assuming Met_analysis has been run
# check appropriate binwidth before making final graphs

All <- FilteredMetGraph(.ind_data = Graphdata, binwidth = 100)
Top <- FilteredMetGraph(.ind_data = Graphdata, .metlevel = topMet, binwidth = 100)
High <- FilteredMetGraph(.ind_data = Graphdata, .metlevel = highMet, binwidth = 500)
Med <- FilteredMetGraph(.ind_data = Graphdata, .metlevel = medMet, binwidth = 40)
Low <- FilteredMetGraph(.ind_data = Graphdata, .metlevel = lowMet, binwidth = 200)
VLow <- FilteredMetGraph(.ind_data = Graphdata, .metlevel = vlowMet, binwidth = 30)


# Filter for spike level and/or valid metabolites in keepMet 
# different graph info

X_label <- "Metabolite"
Y_label <-  "Concentration (nM)"
Graph.title <- "Recovered Spike Amounts"

# will generate graphs for each spike level
Spike <- FilteredMetGraphS(Graphdata)
Spike_keep <- FilteredMetGraphS(Graphdata, metlevel = keepMet)

# Generate Summary PDF

File.name <- paste0(base_filename,"spike recovery", format(Time, "_%Y%m%d_%H%M%S"), ".pdf")

pdf(paste(Exp,"results",File.name, sep = "/"))
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,File.description, pos = 4)
text(1,3.5, Expt.title, pos = 4)
text(1,3, "Updated:", pos = 4)
text(2,3, Today, pos = 4)
text(1,2.5, "Comments:", pos = 4)
text(1,2.3, "Expected spike conc: 0.5, 5 uM", pos = 4)
# text(1,2.1, "samples C1-3 only", pos = 4)
text(1,1.9, "Spike recovery only", pos = 4)
All
Top
High
Med
Low
VLow
FilteredMetGraphS(Graphdata)
FilteredMetGraphS(Graphdata, metlevel = keepMet)
dev.off()
