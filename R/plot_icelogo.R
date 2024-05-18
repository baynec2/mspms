#' plot_icelogo
#'
#' This function plots the cleavage motifs that were enriched relative to background in a manner similar to IceLogo.
#' https://iomics.ugent.be/icelogoserver/resources/manual.pdf
#'
#' @param cleavage_seqs = these are the cleavage sequences that are observed in the experiment
#' @param background_universe = this is a list of all the possible clevage sequences in the peptide library used.
#' @param type = this is the type of visualization you would like to perform, can be either "percent_difference" or "fold_change".
#' @return a ggplot object that shows the motif of the cleavage sequences
#' @export
#'
#' @examples
#'cleavage_seqs = mspms::prepared_for_stats %>%
#' mspms::polish() %>%
#' dplyr::filter(condition == "DMSO",time == 240) %>%
#' dplyr::pull(cleavage_seq)
#'plot_icelogo(cleavage_seqs)
plot_icelogo = function(cleavage_seqs,
                        background_universe = mspms::all_possible_8mers_from_228_library,
                        type = "percent_difference") {
  ## First need to format this such that we have a matrix with the counts at each location of the sequence
  nchar_bg = nchar(background_universe[1])

  # do this for the background first.
  background_universe = tibble::tibble(sequences = background_universe)

  # separate the sequences into individual amino acids
  bg_seq = tidyr::separate(
    background_universe,
    col = 1,
    into = paste0("position", 1:(nchar_bg + 1)),
    sep = "",
    remove = TRUE
  ) %>%
    dplyr::select(-.data$position1)

  # Counting the number in each position
  bg_count = bg_seq %>%
    purrr::map_df(table) %>%
    t()

  #replacing NA with 0
  bg_count[is.na(bg_count)] = 0

  bg_prop = bg_count %>%
    as.data.frame() %>%
    tibble::rownames_to_column("AA") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), prop.table)) %>%
    tibble::column_to_rownames("AA") %>%
    as.matrix()

  # Now doing this for the experimental clevages
  nchar_cleav = nchar(cleavage_seqs[[1]])

  # done for the background first
  cleavage_seqs = tibble::tibble(sequences = cleavage_seqs)

  clev_seq = tidyr::separate(
    cleavage_seqs,
    col = 1,
    into = paste0("position", 1:(nchar_cleav + 1)),
    sep = "",
    remove = TRUE
  ) %>%
    dplyr::select(-.data$position1)

  # counting the number of time each AA appears at each position
  clev_count = clev_seq %>%
    purrr::map_df(table) %>%
    t()

  # replacing NA with 0
  clev_count[is.na(clev_count)] = 0

  clev_prop = clev_count %>%
    as.data.frame() %>%
    tibble::rownames_to_column("AA") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), prop.table)) %>%
    tibble::column_to_rownames("AA") %>%
    as.matrix()


  ### calculating the SD ###
  # defining function per Icelogo manual
  standard_deviation = function(proportion_mat, count_mat) {
    (proportion_mat / colSums(count_mat)) ^ (1 / 2)
  }

  # applying to bg data
  bg_sd = standard_deviation(bg_prop, clev_count)

  ### Calculating Z score ###

  # Defining function per icelogo manual
  zscore = function(cleav_prop, bg_prop, bg_sd) {
    (cleav_prop - bg_prop) / bg_sd
  }

  # applying to our data
  zscores = zscore(clev_prop, bg_prop, bg_sd) %>%
    as.data.frame()

  names(zscores) = paste0("P", c((nchar_cleav / 2):1,
                                 paste0(1:(nchar_cleav / 2), "'")))

  #finding the significant zscores
  sig_zscores = zscores %>%
    tibble::rownames_to_column("AA") %>%
    tidyr::pivot_longer(2:length(.), names_to = "position", values_to = "zscore") %>%
    dplyr::filter(abs(zscore) > 1.96) %>%
    dplyr::mutate(aa_position = paste0(.data$AA, ".", .data$position))

  # Calculating percent difference
  pd = (clev_prop * 100) - (bg_prop * 100) %>%
    as.data.frame()


  names(pd) = paste0("P", c((nchar_cleav / 2):1,
                            paste0(1:(nchar_cleav / 2), "'")))

  ### Determining the type of plot to create ###

  if (type == "percent_difference") {
    sig_final = pd %>%
      tibble::rownames_to_column("AA") %>%
      tidyr::pivot_longer(2:length(.), names_to = "position", values_to = "percent_difference") %>%
      dplyr::mutate(aa_position = paste0(.data$AA, ".", .data$position)) %>%
      dplyr::filter(.data$aa_position %in% sig_zscores$aa_position) %>%
      dplyr::select(-.data$aa_position) %>%
      tidyr::pivot_wider(
        names_from = .data$position,
        values_from = .data$percent_difference,
        names_sort = FALSE
      )

    # Dealing with the problem introduced by the pivot wider function
    missing_cols = 9 - ncol(sig_final)
    missing_data = as.data.frame(matrix(nrow = nrow(sig_final), ncol = missing_cols, NA))
    `%!in%` = Negate(`%in%`)
    missing_names = names(pd)[names(pd) %!in% names(sig_final)]
    names(missing_data) = missing_names

    final = sig_final %>%
      dplyr::bind_cols(missing_data) %>%
      dplyr::relocate(paste0("P", c((nchar_cleav / 2):1,
                                    paste0(
                                      1:(nchar_cleav / 2), "'"
                                    )))) %>%
      tibble::column_to_rownames("AA") %>%
      as.matrix()

    # Plotting the icelogo
    p1 = ggseqlogo::ggseqlogo(final,
                              font = "helvetica_light" ,
                              method = 'custom',
                              seq_type = "AA") +
      ggplot2::scale_x_continuous(labels = paste0("P", c((nchar_cleav / 2):1,
                                                         paste0(
                                                           1:(nchar_cleav / 2), "'"
                                                         ))),
                                  breaks = 1:nchar_cleav) +
      ggplot2::ylab("p < 0.05 Percent Difference")

    return(p1)


  } else if (type == "fold_change") {
    fold_change = (clev_prop * 100) / (bg_prop * 100)

    converted_fc = fold_change

    # Here we run into a problem with infinite values. If a value is 0 in the count, but greater than that in the background - we have a problem.
    # To deal with this, ice logo just shows an infinite value as taking up the entire scale.
    # We could mimic that here by setting the columns with infinite values to the max value of the fold change/ the number of things in the columnns

    max_sum = max(purrr::map_df(as.data.frame(converted_fc), sum), na.rm = T)
    final_converted_fc  = data.frame(row.names = rownames(converted_fc))

    for (i in 1:ncol(converted_fc)) {
      if (0 %in% converted_fc[, i]) {
        num_zero = sum(converted_fc[, i] == 0, na.rm = TRUE)

        col = replace(converted_fc[, i],
                      converted_fc[, i] == 0,
                      (max_sum / num_zero * -1))

        final_converted_fc =  cbind(final_converted_fc, col)
      } else{
        final_converted_fc = cbind(final_converted_fc, converted_fc[, i])
      }
    }
    names(final_converted_fc) = paste0("P", c((nchar_cleav / 2):1,
                                              paste0(1:(nchar_cleav / 2), "'")))
    sig_final = final_converted_fc %>%
      tibble::rownames_to_column("AA") %>%
      tidyr::pivot_longer(2:length(.), names_to = "position", values_to = "fold_change") %>%
      dplyr::mutate(aa_position = paste0(.data$AA, ".", .data$position)) %>%
      dplyr::filter(.data$aa_position %in% sig_zscores$aa_position) %>%
      dplyr::select(-.data$aa_position) %>%
      tidyr::pivot_wider(
        names_from = .data$position,
        values_from = .data$fold_change,
        names_sort = FALSE
      )

    # Dealing with the problem introduced by the pivot wider function
    missing_cols = 9 - ncol(sig_final)
    missing_data = as.data.frame(matrix(nrow = nrow(sig_final), ncol = missing_cols, NA))
    `%!in%` = Negate(`%in%`)
    missing_names = names(final_converted_fc)[names(final_converted_fc) %!in% names(sig_final)]

    names(missing_data) = missing_names
    final = sig_final %>%
      dplyr::bind_cols(missing_data) %>%
      dplyr::relocate(paste0("P", c((nchar_cleav / 2):1,
                                    paste0(
                                      1:(nchar_cleav / 2), "'"
                                    ))), ) %>%
      tibble::column_to_rownames("AA") %>%
      as.matrix()

    # Plotting the motif
    p1 = ggseqlogo::ggseqlogo(final,
                              font = "helvetica_light" ,
                              method = 'custom',
                              seq_type = "AA") +
      ggplot2::scale_x_continuous(labels = paste0("P", c((nchar_cleav / 2):1,
                                                         paste0(
                                                           1:(nchar_cleav / 2), "'"
                                                         ))),
                                  breaks = 1:nchar_cleav) +
      ggplot2::ylab("p < 0.05 fold change")

    return(p1)

  } else {
    stop("type must be either 'percent_difference' or 'fold_change'")
  }
}


