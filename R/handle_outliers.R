#' handle_outliers
#'
#' This function replaces outliers for MSP-MS data with NA based on the results of dixon tests as implemented in the outliers package
#'
#' @param normalyzed_data = this is the normalized data. Intended to be prepared using mspms::normalyze()
#' @param design_matrix = this is the design_matrix containing "sample","group","condition","time" columns. Intended to be read in using readr::read_csv()
#'
#' @return a data frame in the long format where outliers according to a dixon test have been handled.
#' @export
#' @examples
#'
#' # Here we handle the outliers from the test data
#'handled_outliers = handle_outliers(mspms::normalyzed_data,mspms::design_matrix)

handle_outliers = function(normalyzed_data,design_matrix){

  # dealing with no visible binding for global variable ‘.’ NOTE
  . = NULL

  # Putting into the long format, adding group information, based on sample_time_replicate naming convention
  index = which(names(normalyzed_data) == "Protein Accession")+1

  long_data = normalyzed_data %>%
    tidyr::pivot_longer(dplyr::all_of(index:length(.)),names_to = "sample") %>%
    dplyr::inner_join(design_matrix,by = "sample") %>%
    dplyr::filter(!is.na(.data$value))


  # initalizing the data frame
  out = data.frame()

  # Looping through each group to check for outliers.
  for(i in unique(long_data$group)){

    #filtering data to only include the samples with the same condition and time
    f = dplyr::filter(long_data, .data$group == i)

    # Looping though each peptide within each group
    for(j in unique(long_data$Peptide)){
      f2 = dplyr::filter(f,.data$Peptide == j )
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
    dplyr::select(-.data$group,-.data$condition,-.data$time) %>%
    tidyr::pivot_wider(names_from = .data$sample, values_from = .data$value) %>%
    tidyr::pivot_longer(index:length(.),names_to = "sample_id") %>%
    tibble::as_tibble()


  return(NAs_added_back)
}






