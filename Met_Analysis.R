# Graph quantified metabolites
# Load libraries and other functions from Metabolomics_01.R

# check current directory
here()

# Use Met_DataPrep.R to set up data for analysis
# If recovery info available use data from Met_DataPrep_2.R

# transfer data files to project folder
# this is set up as a project so copy files into subfolders in the project folder

# Expt specific info
source("NM-08/NM-08_info.R") # Separate R script for experiment specific information

suffix <- "_20201127_054153" # select specific timestamped file, otherwise set as NULL
Note <- "NR_NaR only, Curve regression: Area, linear through 0, weighting = none" # comment for PDF
Note2 <-  "Recov corr using H std spikes" # comment for PDF

recovery <- "Y" # Y or N if recovery info is avail
deconv <- "N" # Y or N if using deconvoluted data

if(deconv == "Y") {
  comp_datafile <- paste0(base_filename, "Deconv_compiled", suffix, ".csv") # output from Deconv_DataPrep.R
} else {
if(recovery == "Y"){
  comp_datafile <- paste0(base_filename, "compiled_corr", suffix, ".csv") # output from Met_DataPrep_2.R
} else {
comp_datafile <- paste0(base_filename, "compiled", suffix, ".csv") # Output from Met_DataPrep.R
}
}


# PDF output info, change as needed
descrip <- "NR Aex NAD metabolome"
Graph.title <- "NAD metabolites (selected transitions)"
Expt.title <- paste(Exp, descrip)
File.name <- paste(Exp, "_", descrip,".pdf",sep = "")
File.description <- "Metabolome analysis"
Today <- Sys.Date() # Today's date format:yyyy-mm-dd
Time <- Sys.time() 

# Filter and grouping lists, change as needed
# ids_interest <- c("52","1816","1819","258288","258290") # untested
grouping_base <-  c("sample_name2", "component_name", "metabolite")
Indiv_grouping <- c("id_number", "component_name", "metabolite")
grouping <- c("component_name", "metabolite")
Y_label <- "Concentration (pmol/mg)"
# reorder_list <- c(7,13,18,22,33)
#reorder_list <- c("N","L","M","H")

# load data files
data <- read_csv(paste(Exp,"data",comp_datafile, sep = "/"), na = "N/A")
MetLists <- read_csv(paste(Exp,"data",met_datafile, sep = "/"), na = "N/A")



# additional changes to data
data <- data %>% mutate (id_number = factor(id_number)) # id_number not always in expected format

if("age" %in% colnames(data)) {data <- data %>% mutate (age = factor(age))}

if("use" %in% colnames(data)) {data <- data %>% filter(use == "Y")}

# set in expt info.R script
if(comb_vars == "Y") {
  comb_name <- X_var
  data <- data %>% unite (!!comb_name, c(var1, var2), na.rm = TRUE)  # merging some sample info
}

MetLists <- MetLists %>%
  clean_names()

# NAD Metabolites to analyse

usedMet <- MetLists %>% filter(used_met == "Y") 
usedMet <- usedMet[['component_name']]

# Metabolites separated out due to different amounts in samples
# Lists change depending on experiment conditions (eg: tissue type)
# Check results with usedMet and then adjust lists - enter data into met_datafile and reload

if(recovery == "Y") {
topMet <-  MetLists %>% filter(used_met_rec_level == "top") 
topMet <- topMet[['component_name']]

highMet <-  MetLists %>% filter(used_met_rec_level == "high") 
highMet <- highMet[['component_name']]

medMet <-  MetLists %>% filter(used_met_rec_level == "med") 
medMet <- medMet[['component_name']]

lowMet <-  MetLists %>% filter(used_met_rec_level == "low") 
lowMet <- lowMet[['component_name']]

vlowMet <-  MetLists %>% filter(used_met_rec_level == "vlow") 
vlowMet <- vlowMet[['component_name']]

} else {

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

}


# Standard mix check
Std_data <- data %>% 
  filter (category == "Standard") %>%
  filter (component_name %in% usedMet) %>% 
  group_by_at(vars(all_of(grouping_base))) %>%
  summarise (Avg = mean(calc_conc), Err = se(calc_conc))

if(deconv == "Y") {
Std_BarGraph <- Std_data %>%
  ggplot(aes(x = sample_name2, y = Avg, fill = metabolite)) + 
  facet_wrap(~ component_name) + 
  geom_col(position = "dodge") +
  #geom_errorbar(width = 0.5) +
  labs(x = "Standard Mix",     
       y = "Concentration (nM)",
       title = "NAD metabolites: Std controls" )+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
} else {
Std_BarGraph <- Std_data %>%
  ggplot(aes(x = sample_name2, y = Avg)) + 
  facet_wrap(~ component_name) + 
  geom_col() +
  #geom_errorbar(width = 0.5) +
  labs(x = "Standard Mix",     
       y = "Concentration (nM)",
       title = "NAD metabolites: Std controls" )+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # adjusts orientation of x axis labels
}
plot(Std_BarGraph)

# Data prep for export and graphs
# Experiment specific options to change as needed

# average over individual samples
# used for export for Prism and for dot plot graphs

if(recovery == "Y") {
Indiv_data <- data %>% 
  filter (category == "Sample") %>%
  filter (pmol_mg > 0) %>% # corrected and valid concentrations
  group_by_at(vars(all_of(Indiv_grouping), !!X_var)) %>%
  summarise (Avg = mean(corr_pmol_mg), Err = se(corr_pmol_mg))

# average over grouped samples
# used for column graphs

Group_data <- data %>% 
  filter (category == "Sample") %>%
  filter (pmol_mg > 0) %>% # corrected and valid concentrations
  group_by_at(vars(all_of(grouping),!!X_var)) %>%
  summarise (Avg = mean(corr_pmol_mg), Err = se(corr_pmol_mg))

} else {
  
  Indiv_data <- data %>% 
    filter (category == "Sample") %>%
    filter (pmol_mg > 0) %>% # corrected and valid concentrations
    group_by_at(vars(all_of(Indiv_grouping), !!X_var)) %>%
    summarise (Avg = mean(pmol_mg), Err = se(pmol_mg))
  
  # average over grouped samples
  # used for column graphs
  
  Group_data <- data %>% 
    filter (category == "Sample") %>%
    filter (pmol_mg > 0) %>% # corrected and valid concentrations
    group_by_at(vars(all_of(grouping),!!X_var)) %>%
    summarise (Avg = mean(pmol_mg), Err = se(pmol_mg))
}

#export for Prism - all metabolites
if(deconv == "Y") {
Indiv_data <- Indiv_data %>%
  filter (component_name %in% usedMet) %>% 
  unite(met_comp, c(metabolite, component_name), remove = TRUE, na.rm = TRUE)

data_print <- Indiv_data %>%
  ungroup %>%
  mutate(Err = NULL) %>%
  spread(met_comp, Avg)

} else {
  
Indiv_data <- Indiv_data %>%
  filter (component_name %in% usedMet)

data_print <- Indiv_data %>%
  ungroup %>%
  select(-component_name) %>%
  mutate(Err = NULL) %>%
  spread(metabolite, Avg)
}

# rearrange for easier copy-paste, sort as needed
data_print <- data_print %>% 
  pivot_longer(-c(id_number,!!X_var)) %>% 
  pivot_wider(names_from = id_number, values_from = value)

if(deconv == "Y") {
  Prismfile <- paste0(base_filename, "results_deconv", format(Time, "_%Y%m%d_%H%M%S"), ".csv")
} else {
if(recovery == "Y"){
Prismfile <- paste0(base_filename, "results_corr", format(Time, "_%Y%m%d_%H%M%S"), ".csv")
} else {
  Prismfile <- paste0(base_filename, "results", format(Time, "_%Y%m%d_%H%M%S"), ".csv")
}
}

write.csv(file = paste(Exp,"results",Prismfile, sep = "/"), data_print)

# graphs

# Graph to check on metabolite levels and if data needs to be reordered

GraphCheck <- Group_data %>%
  ggplot(aes(x = !!X_var, y = Avg, ymin = Avg - Err, ymax = Avg + Err, fill = metabolite)) + 
  facet_wrap(~ component_name) + 
  geom_col(position = "dodge") +
  geom_errorbar(width = 0.5, position = "dodge") +
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

DotAll <- FilteredMetGraph(.ind_data = Indiv_graphdata, binwidth = 0.1)

# Cycling through Metabolite levels
# check appropriate binwidth before making final graphs

if (length(topMet) > 0) {Top <- FilteredMetGraph(.ind_data = Indiv_graphdata, .gr_data = Group_graphdata,
                        .metlevel = topMet, binwidth = 2)}
if (length(highMet) > 0) {High <- FilteredMetGraph(.ind_data = Indiv_graphdata, .gr_data = Group_graphdata,
                        .metlevel = highMet, binwidth = 1)}
if (length(medMet) > 0) {Med <- FilteredMetGraph(.ind_data = Indiv_graphdata, .gr_data = Group_graphdata,
                         .metlevel = medMet, binwidth = 0.2)}
if (length(lowMet) > 0) {Low <- FilteredMetGraph(.ind_data = Indiv_graphdata, .gr_data = Group_graphdata,
                        .metlevel = lowMet, binwidth = 0.05)}
if (length(vlowMet) > 0) {VLow <- FilteredMetGraph(.ind_data = Indiv_graphdata, .gr_data = Group_graphdata,
                        .metlevel = vlowMet, binwidth = 0.007)}

# Generate Summary PDF


if(recovery == "Y"){
  File.name <- paste0(base_filename, "quant changes_corr", format(Time, "_%Y%m%d_%H%M%S"), ".pdf")
} else {
  File.name <- paste0(base_filename,"quant changes", format(Time, "_%Y%m%d_%H%M%S"), ".pdf")
}

pdf(paste(Exp,"results",File.name, sep = "/"))
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,File.description, pos = 4)
text(1,3.5, Expt.title, pos = 4)
text(1,3, "Updated:", pos = 4)
text(2,3, Today, pos = 4)
text(1,2.5, "Comments:", pos = 4)
text(1,2.3, Note, pos = 4)
text(1,2.1, Note2, pos = 4)
plot(Std_BarGraph)
BarAll
DotAll
if(exists("Top")) {Top}
if(exists("High")) {High}
if(exists("Med")) {Med}
if(exists("Low")) {Low}
if(exists("VLow")) {VLow}
dev.off()