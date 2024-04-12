#' handle_outliers
#'
#' @param normalyzed_data = this is the normalzed data.
#'
#' @return a data frame in the long format where outliers according to a dixon test have been handled. 
#' @export
#'
#' @examples
handle_outliers = function(normalyzed_data){
  
  # Putting into the long format, adding group information, based on sample_time_replicate naming convention
  
  long_data = normalyzed_data %>% 
    tidyr::pivot_longer(30:length(.)) %>% 
    tidyr::separate(name,sep = "_",into = c("sample","time","replicate")) %>% 
    dplyr::mutate(sample_id = paste0(sample,"_",time,"_",replicate),
           group = paste0(sample,".",time)) %>% 
    dplyr::filter(!is.na(value))
  
  
  # initalizing the data frame
  out = data.frame()
  
  # Looping through each group to check for outliers.
  for(i in unique(long_data$group)){
    
    #filtering data to only include the samples with the same condition and time
    f = dplyr::filter(long_data, group == i)
    
    # Looping though each peptide within each group
    for(j in unique(long_data$Peptide)){
      f2 = dplyr::filter(f,Peptide == j )
      # # only doing tests if there are between 3 and 30 values
       if(length(f2$value) > 3 & length(f2$value) < 30){
      #   
      #   #is there an outlier?
        stats = outliers::dixon.test(f2$value)$p.value
        
        #what would the outlier be if there was one?
        outlier = outliers::outlier(f2$value)
       
        # If there is a outlier, let's change the value to 0 for now, we will change it to NA after we impute. 
        if(stats < 0.05){
          f2[f2$value == outlier,"value"] = 0
        }
        # #adding the cv
        # f2 = f2 %>% 
        #     dplyr:: mutate(cv = sd(value,na.rm = T)/mean(value,na.rm = T))
       }
        # let's build our final data frame
        out = dplyr::bind_rows(out,f2)

      }
  }


  NAs_added_back = out %>% 
    #dplyr::filter(!is_outlier) %>% 
    dplyr::select(-sample,-group,-time,-replicate,-group) %>% 
    tidyr::pivot_wider(names_from = sample_id, values_from = value) %>% 
    tidyr::pivot_longer(30:length(.),names_to = "sample_id")
  
  
  return(NAs_added_back)
}

  
  
  
  
  
  
  