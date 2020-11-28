# Graphs metabolites with area data only
# Load libraries and other functions from Metabolomics_01.R
# 

# check current directory
here()

# Use OtherMet_DataPrep.R to set up data for analysis
# This is for graphing metabolites without standard curve info

# transfer data files to project folder
# this is set up as a project so copy files into subfolders in the project folder

# Expt specific info
source("NM-08/NM-08_info.R") # Separate R script for experiment specific information

suffix <- "_20201117_011337" # select specific timestamped file, otherwise set as NULL

comp_datafile <- paste0(base_filename, "noQuant_compiled", suffix,".csv") # Output from OtherMet_DataPrep.R
descrip <- "NR Aex relative changes"

Graph.title <- "Un-quantified Metabolites"

# PDF output info, change as needed
Expt.title <- paste(Exp, descrip)

File.description <- "Metabolome analysis"
Today <- Sys.Date() # Today's date format:yyyy-mm-dd
Time <- Sys.time() # use to set a timestamp in filename

# Filter and grouping lists, change as needed
# ids_interest <- c("52","1816","1819","258288","258290") # untested
Y_label <- "Peak area/mg"
grouping_base <-  c("sample_name2", "component_name", "metabolite")
Indiv_grouping <- c("id_number", "component_name", "metabolite")
#reorder_list <- c(7,13,18,22,33)
#reorder_list <- c("N","L","M","H")

# load data file
data <- read_csv(paste(Exp,"data",comp_datafile, sep = "/"), na = "N/A")
MetLists <- read_csv(paste(Exp,"data",met_datafile, sep = "/"), na = "N/A")

# additional changes to data
data <- data %>% mutate (id_number = factor(id_number)) # needed to fix issues with NAs after importing data

if("age" %in% colnames(data)) {
  data <- data %>% mutate (age = factor(age))
}

if("use" %in% colnames(data)) {
  data <- data %>% filter(use == "Y")
}

if(comb_vars == "Y") {
  comb_name <- X_var
  data <- data %>% unite (!!comb_name, c(var1, var2), na.rm = TRUE)  # merging some sample info
}

MetLists <- MetLists %>%
  clean_names()

# NAD Metabolites to analyse
# usedMet <- MetLists %>% 
#   filter(used_met == "Y") # use this for metabolites with standard curves 
# usedMet <- usedMet[['component_name']]

usedMet <- MetLists %>% filter(other_met == "Y") # use this for metabolites with out standard curves
usedMet <- usedMet[['component_name']]

# Metabolites separated out due to different amounts in samples
# Lists change depending on experiment conditions (eg: tissue type)
# Check results with usedMet and then adjust lists - enter data into met_datafile and reload
# 

topMet <-  MetLists %>% filter(other_met_level == "top") 
topMet <- topMet[['component_name']]

highMet <-  MetLists %>% filter(other_met_level == "high") 
highMet <- highMet[['component_name']]

medMet <-  MetLists %>% filter(other_met_level == "med") 
medMet <- medMet[['component_name']]

lowMet <-  MetLists %>% filter(other_met_level == "low") 
lowMet <- lowMet[['component_name']]

vlowMet <-  MetLists %>% filter(other_met_level == "vlow") 
vlowMet <- vlowMet[['component_name']]


# Data prep for export and graphs
# Experiment specific options to change as needed

# # average over a group using uncorrected concentrations ie: don't have wet weight data
# 
# Group_raw_data <- data %>% 
#   filter (category == "Sample" ) %>% 
#   group_by_at(vars(all_of(grouping_base),!!X_var)) %>%
#   summarise (Avg = mean(calc_conc), Err = se(calc_conc)) 
# 

# average over individual samples
# used for export for Prism and for dot plot graphs

Indiv_data <- data %>% 
  filter (category == "Sample") %>%
  filter (area_mg > 0) %>% # corrected and valid concentrations
  group_by_at(vars(all_of(Indiv_grouping), !!X_var)) %>%
  summarise (Avg = mean(area_mg), Err = se(area_mg))

Indiv_data <- Indiv_data %>%
  filter (component_name %in% usedMet)

# average over a group using uncorrected concentrations ie: don't have wet weight data

Group_data <- data %>% 
  filter (category == "Sample" ) %>%
  filter (area_mg > 0) %>% # corrected and valid concentrations
  group_by_at(vars(all_of(grouping_base),!!X_var)) %>%
  summarise (Avg = mean(area_mg), Err = se(area_mg)) 

Group_data <- Group_data %>%
  filter (component_name %in% usedMet)

#export for Prism - all metabolites
data_print <- Indiv_data %>%
  ungroup %>%
  select(-component_name) %>%
  mutate(Err = NULL) %>%
  spread(metabolite, Avg)

# rearrange for easier copy-paste, sort as needed
data_print <- data_print %>% 
  pivot_longer(-c(id_number,!!X_var)) %>% 
  pivot_wider(names_from = id_number, values_from = value)

Prismfile <- paste0(base_filename, "noQuant_results", format(Time, "_%Y%m%d_%H%M%S"),".csv")
write.csv(file = paste(Exp,"results",Prismfile, sep = "/"), data_print)

  
# Graph to check on metabolite levels and if data needs to be reordered

GraphCheck <- Group_data %>%
  ggplot(aes(x = !!X_var, y = Avg, ymin = Avg - Err, ymax = Avg + Err)) + 
  facet_wrap(~ component_name) + 
  geom_col() +
  geom_errorbar(width = 0.5) +
  labs(x = X_label,     
       y = Y_label,
       title = Graph.title)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # adjusts orientation of x axis labels

plot(GraphCheck)

# Reorder data if needed, create Name column if not
if(exists("reorder_list")) {
  Indiv_graphdata <- DataReorder(Indiv_data, !!X_var, reorder_list)
  Group_graphdata <- DataReorder(Group_data, !!X_var, reorder_list)
} else {
  Indiv_graphdata <- Indiv_data %>%
    mutate(Name = !!X_var)
  Group_graphdata <- Group_data %>%
    mutate(Name = !!X_var)
}

# Colour selections for graphs
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
myPalette = getPalette(8) # not used?
colourCount = length(unique(Indiv_graphdata$id_number))
myPalette = getPalette(colourCount)

BarAll <- FilteredMetGraph(.gr_data = Group_graphdata)

DotAll <- FilteredMetGraph(.ind_data = Indiv_graphdata, binwidth = 10000)

# Cycling through Metabolite levels

# # Bar graphs only
# TopBar <- FilteredMetGraph(.gr_data = Group_graphdata,.metlevel = topMet)
# 
# HighBar <- FilteredMetGraph(.gr_data = Group_graphdata,.metlevel = highMet)
# 
# MedBar <- FilteredMetGraph(.gr_data = Group_graphdata,.metlevel = medMet)
# 
# LowBar <- FilteredMetGraph(.gr_data = Group_graphdata,.metlevel = lowMet)
# 
# VLowBar <- FilteredMetGraph(.gr_data = Group_graphdata,.metlevel = vlowMet)

# Combined graphs 
Top <- FilteredMetGraph(.ind_data = Indiv_graphdata, .gr_data = Group_graphdata,
                        .metlevel = topMet, binwidth = 70000)
High <- FilteredMetGraph(.ind_data = Indiv_graphdata, .gr_data = Group_graphdata,
                         .metlevel = highMet, binwidth = 6000)
Med <- FilteredMetGraph(.ind_data = Indiv_graphdata, .gr_data = Group_graphdata,
                        .metlevel = medMet, binwidth = 10000)
Low <- FilteredMetGraph(.ind_data = Indiv_graphdata, .gr_data = Group_graphdata,
                        .metlevel = lowMet, binwidth = 400)
VLow <- FilteredMetGraph(.ind_data = Indiv_graphdata, .gr_data = Group_graphdata,
                         .metlevel = vlowMet, binwidth = 20)


# Generate Summary PDF

File.name <- paste(Exp, "_", descrip, format(Time, "_%Y%m%d_%H%M%S"), ".pdf",sep = "")

pdf(paste(Exp,"results",File.name, sep = "/"))
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,File.description, pos = 4)
text(1,3.5, Expt.title, pos = 4)
text(1,3, "Updated:", pos = 4)
text(2,3, Today, pos = 4)
text(1,2.5, "Comments:", pos = 4)
text(1,2.3, "Metabolite amounts based on peak area", pos = 4)
#text(1,2.1, "Using one injection/sample only", pos = 4)
BarAll
DotAll
# TopBar
# HighBar
# MedBar
# LowBar
# VLowBar
Top
High
Med
Low
VLow
dev.off()
