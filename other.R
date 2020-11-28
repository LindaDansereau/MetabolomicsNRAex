# Print Standard controls only
File.name <- paste0(base_filename,"Standard controls", format(Time, "_%Y%m%d_%H%M%S"), ".pdf")

pdf(paste(Exp,"results",File.name, sep = "/"))
plot(Std_BarGraph)
dev.off()

# Print Spike recovery All graph only
File.name <- paste0(base_filename,"spike recovery", format(Time, "_%Y%m%d_%H%M%S"), ".pdf")

pdf(paste(Exp,"results",File.name, sep = "/"))
All
dev.off()

# manual range check

x <-data_2 %>% 
  filter(category == "Sample") %>% 
  #filter(sample_name != "NM-08_NR_Aex_34") %>% 
  filter(component_name == "MeNAM1")

hist(x$area)

max(x$area)

length(which(x$area > 38044827))

# TODO:  
# Warning: some scripts changed for deconvolution data and probably won't work for regular data
# make graphs easier to adapt to different projects

# - regressions for standard curves: run Standard_DataPrep.R, copy data into Prism, add info to StandardCurves.csv
# - run Deconv_DataPrep.R with new StandardCurves
# - goto Met_Analysis.R: run with deconv = Y, run to StdGraph and then .csv export
# - goto SpikeRecovery.R: run to .csv export 
# - goto Met_DataPrep_2.R: run to end 
# - goto Met_Analysis.R: run with recov = Y, run to StdGraph and then use deconv .csv export
# - copy results into Prism

# - NMN and NaAD reading higher - why?

# - new NR and NaR standard curves; update CSV standard curve file - done
# - Met_DataPrep.R using recovery data and first part of sample data
# - Met_Analysis.R
# - SpikeRecovery.R - confirm valid recovery info
# - pass on interim results to Andy/Ben
# - Met_DataPrep_2.R - 
# - Met_Analysis.R - recov = Y
# - send Excel file to Andy


