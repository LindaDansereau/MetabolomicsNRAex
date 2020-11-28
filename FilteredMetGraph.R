# output:
# various graphs
# 
# input: 
# data sources: .ind_data and/or .gr_data 
# list of metabolites: .metlevel (eg: topMet, default = usedMet)
# ... : binwidth; 
#
# global environment variables: 
# myPalette
# X_label
# Y_label
# Graph.title
# contents of met levels with component names


FilteredMetGraph <- function(.ind_data, .gr_data, ... , .metlevel = usedMet){
  
  ellipsis::check_dots_used()
  
  metlevel <- enquo(.metlevel)
  
  if(!missing(.ind_data)) {
    Indiv_filterdata <- .ind_data %>%
          filter(.data$component_name %in% (!!metlevel)) } # filter individual dataframe
  
  if(!missing(.gr_data)) {
    Group_filterdata <- .gr_data %>%
        filter(.data$component_name %in% (!!metlevel)) # filter group dataframe
    } 

  # make combo graph with both individual and group data
  if(!missing(.ind_data) & !missing(.gr_data)) {
    ComboGraph <- ggplot() +
      geom_col (data = Group_filterdata,
                aes(x = Name, y = Avg)) +
      geom_errorbar(data = Group_filterdata,
                    aes (x = Name, ymin = Avg-Err, ymax = Avg+Err), width = 0.5) +
      geom_dotplot (data = Indiv_filterdata,
                    aes (x = Name, y = Avg, fill = id_number), 
                    binwidth = ..., binaxis = "y", stackdir = "center") +
      scale_fill_manual(values = myPalette) +
      facet_wrap( ~ metabolite) + 
      labs(x = X_label,     
           y = Y_label,
           title = Graph.title )+
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) # adjusts orientation of x axis labels

    
    plot(ComboGraph) } else {
      
  # make column graph with group data, with option for no errorbars if Err not available
  if(!missing(.gr_data)) {
    if ("Err" %in% colnames(Group_filterdata)) {
      ColGraph <- ggplot() +
        geom_col (data = Group_filterdata,
                  aes(x = Name, y = Avg)) +
        geom_errorbar(data = Group_filterdata,
                      aes (x = Name, ymin = Avg-Err, ymax = Avg+Err), width = 0.5) +
        scale_fill_manual(values = myPalette) +
        facet_wrap( ~ metabolite) + 
        labs(x = X_label,     
             y = Y_label,
             title = Graph.title )+
        theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) # adjusts orientation of x axis labels

      plot(ColGraph)} else {
    
    ColGraph <- ggplot() +
    geom_col (data = Group_filterdata,
              aes(x = Name, y = Avg)) +
    scale_fill_manual(values = myPalette) +
    facet_wrap( ~ metabolite) +
    labs(x = X_label,
         y = Y_label,
         title = Graph.title )+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # adjusts orientation of x axis labels
    
    plot(ColGraph)  }

    } else {
  
  # make dot graph with individual data
  if(!missing(.ind_data)) {
    DotGraph <- Indiv_filterdata  %>%
      ggplot(aes(x = Name, y = Avg, fill = id_number)) + 
      facet_wrap(~ metabolite) + 
      geom_dotplot(binwidth = ..., binaxis = "y", stackdir = "center") +  # adjust binwidth to alter dot size
      scale_fill_manual(values = myPalette) +
      geom_hline(yintercept = 0, colour = "black", size = 0.5) +
      labs(x = X_label,     
           y = Y_label,
           title = Graph.title)+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) # adjusts orientation of x axis labels
    
    plot(DotGraph) }
    
  }
  }
  
} # end of function
  