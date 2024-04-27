#' impute
#'
#' This function imputes data using according to the truncated normal distribution as implemented in truncnorm.
#'
#' @param outlier_handled_data = this is the data that has been normalized, outlier removed, and is ready for imputation.
#' @param noise = this is the level of "noise" that we want to use to create the distribution for imputation. Corresponds to the percentage of the data at the lower end of the distribution that we want to use for imputation.
#' @param nsd = this is the number of standard deviations used to calculate the max and min values for the truncnorm distribution.
#'
#' @return a data frame in the wide format where NAs have been imputed.
#' @export
#'
#' @examples
#' impute(mspms::outlier_handled_data,noise = 0.05,nsd = 1)
impute = function(outlier_handled_data,noise=0.05,nsd=1){
  # taking the lowest end of the values observed to use as a list of values to compose the distribution we will be imputing with.
  values_for_impute = outlier_handled_data %>%
    dplyr::filter(value != 0 ) %>%
    dplyr::arrange(-value) %>%
    dplyr::pull(value) %>%
    tail(noise * length(.))

  # Maximum-likelihood fitting of univariate distributions
  fit.norm = MASS::fitdistr(values_for_impute, "normal")

  mean = fit.norm$estimate[[1]]
  sd = fit.norm$estimate[[2]]

  # mean + sd = max
  max = fit.norm$estimate[[1]]+nsd*fit.norm$estimate[[2]]

  # mean - sd = min
  min = fit.norm$estimate[[1]]-nsd*fit.norm$estimate[[2]]


  #setting seed to make imputation reproducible.
  set.seed(2)


  #Need to get our outlier data in a matrix form, we could then pivot it back


  # Imputing values for NAs
  imputed = c()
  # Looping through each value
  for(i in outlier_handled_data$value){
    if(is.na(i)){
      i = round(truncnorm::rtruncnorm(1,
                                      mean = mean,
                                      sd = sd,
                                      a = min,
                                      b = max))
    }

    imputed = c(imputed,i)
  }

  # altering data frame. Putting in the imputed values
  imputed_data = dplyr::bind_cols(outlier_handled_data,imputed_values = imputed)

  #making our outliers, that are coded as 0s, NAs
  imputed_data$imputed_values[imputed_data$imputed_values == 0] = NA_real_


  # Putting in the wide format.
  output = imputed_data %>%
    dplyr::select(-value) %>%
    tidyr::pivot_wider(names_from = "sample_id",
                       values_from = imputed_values) %>%
    dplyr::mutate(Peptide = gsub("\\.","_",Peptide),
                  Peptide_no_cleavage = gsub("_","",Peptide),.after = Peptide) %>%
    tibble::as_tibble()

  return(output)

}
