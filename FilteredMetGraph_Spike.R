# output:
# Bar graphs of spike levels in individual samples
# input: 
# Data
# Spike level
# list of metabolites as component name: usedMet (default), keepMet
# 
# global environment variables: 
# myPalette
# X_label
# Y_label
# Graph.title
# contents of met levels with component names
# reorder_list

FilteredMetGraphS <- function(ind_data, metlevel = usedMet){
  metlevel <- enquo(metlevel)
  
  Indiv_filterdata <- ind_data %>%
    filter(.data$component_name %in% (!!metlevel))

  if(exists("reorder_list")) {
  Recov_graph <- Indiv_filterdata %>%
    ggplot(aes(x = metabolite, y = Avg, fill = id_number)) + 
    geom_col(position = 'dodge') +
    scale_fill_manual(values = myPalette) +
    facet_wrap( ~ (Level = factor(Level, levels = reorder_list))) +
    labs(x = X_label,     
         y = Y_label,
         title = Graph.title)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) # adjusts orientation of x axis labels
  
  plot(Recov_graph)
  
  for(i in reorder_list) {
    New_filterdata <- Indiv_filterdata %>%
      filter (Level == i)
    
    Graphname <- paste("Recov_graph", i)
    New.Title <- paste(Graph.title, i)
    
    Graphname <- New_filterdata %>%
      ggplot(aes(x = metabolite, y = Avg, fill = id_number)) + 
      geom_col(position = 'dodge') +
      scale_fill_manual(values = myPalette) +
      labs(x = X_label,     
           y = Y_label,
           title = New.Title)+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) # adjusts orientation of x axis labels
    
    plot(Graphname)
  }
  } else {
    Recov_graph <- Indiv_filterdata %>%
      ggplot(aes(x = metabolite, y = Avg, fill = id_number)) + 
      geom_col(position = 'dodge') +
      scale_fill_manual(values = myPalette) +
      labs(x = X_label,     
           y = Y_label,
           title = Graph.title)+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) # adjusts orientation of x axis labels
    
    plot(Recov_graph)
  }
  
  
} # end of function