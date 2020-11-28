

Exp <- "NM-08"
datafile <- "NM-08_NR-NaR_final.csv" # Multiquant data export
base_filename <- "NM-08_NRAex_NR_NaR_" # prefix label for all output files

stdinfofile <- "NM-08_StandardInfo.csv" # Standard info for Standard_DataPrep.R
stdconcfile <- "NM-08_StandardConc.csv" # Standard concentrations for Standard_DataPrep.R

stddatafile <- "NM-08_StandardCurves.csv" # Standard curve info if avail
#stddatafile <- "NM-08_StandardCurves_highMeNAM.csv" # Standard curve file with different MeNAM curve
sampleinfofile <- "NM-08_SampleInfo.csv" # Sample info if avail
met_datafile <- "NM-08_MetaboliteLists.csv" # lists for filtering metabolites

deconvdatafile <- "NaAD_NAD_deconv_NM-08_samples.csv" # MATLAB output file


# info for main analysis
X_var <- sym("tr_tp")
X_label <- "Treatment timepoint"
reorder_list <- c("Placebo_O","Placebo_pre","Placebo_post","Placebo_3h post","NR_O","NR_pre","NR_post","NR_3h post")

# if combining info from multiple columns for X-var
comb_vars <- "Y" # Y or N 
var1 <- sym("treatment") # NULL or sym("var")
var2 <- sym("timepoint") # NULL or sym("var")


# Standard Curve info: copy-paste into Met_DataPrep line 40
# filter (standard_curve == "NADmetA" | standard_curve == "NADmetB" | standard_curve == "NADmetC")

# Spike recovery: 
# reorder_list <- c("L","H") # copy-paste into line 32
