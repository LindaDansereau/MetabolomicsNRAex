# output: dataframe with reordered levels
# input: data = starting data, col_name = column to reorder, new_order = levels
# give function call to variable to bring to global environment
# possible to just make new column without reordering?

DataReorder <- function(old_data,col_name,new_order){
  col_name = enquo(col_name)
  
  old_data %>%
    mutate (Name = factor((!!col_name), levels = new_order))

}
