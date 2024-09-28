#' calc_AA_count_of_motif
#'
#' Calculate the counts of amino acids at each position of a motif for all the
#' sequences in a vector.
#'
#' @param cleavage_motif a vector of cleavage motifs
#' @return a matrix of counts
#' @keywords internal
calc_AA_count_of_motif <- function(cleavage_motif) {
    # calculating number of characters
    nchar <- nchar(cleavage_motif[1])

    # do this for the background first.
    sequences <- tibble::tibble(cleavage_motif = cleavage_motif)

    # separate the sequences into individual amino acids
    AAs <- tidyr::separate(
        sequences,
        col = 1,
        into = c("exclude", paste0("position", seq_len(nchar))),
        sep = "",
        remove = TRUE
    ) %>%
        dplyr::select(-"exclude")

    # Define the desired order.
    # This ensures that the amino acids will always be in the same order.
    desired_order <- c(
        "A", "D", "E", "F", "G", "H", "I", "K", "L", "n",
        "N", "P", "Q", "R", "S", "T", "V", "W", "X", "Y"
    )

    # Counting the number in each position
    count <- AAs %>%
        purrr::map_df(table)

    # Are any columns missing?
    `%!in%` <- Negate(`%in%`)
    missing_cols <- desired_order[desired_order %!in% names(count)]

    # building a df with missing columns
    missing_matrix <- matrix(nrow = nchar, ncol = length(missing_cols))
    colnames(missing_matrix) <- missing_cols
    missing_df <- tibble::as_tibble(missing_matrix)

    count_matrix <- dplyr::bind_cols(count, missing_df) %>%
        dplyr::relocate(dplyr::all_of(desired_order)) %>%
        t()

    # replacing NA with 0
    count_matrix[is.na(count_matrix)] <- 0

    colnames(count_matrix) <- paste0("P", c(
        seq(from = (nchar / 2), to = 1),
        paste0(seq_len(nchar / 2), "'")
    ))

    return(count_matrix)
}

#' calc_AA_prop_of_motif
#'
#' Calculate the proportion of amino acids at each position in a vector of
#' motifs.
#'
#' @param count_matrix this is a matrix of the counts of cleavage motifs
#'
#' @return a matrix with proportions of counts.
#' @keywords internal
calc_AA_prop_of_motif <- function(count_matrix) {
    prop <- count_matrix %>%
        as.data.frame() %>%
        tibble::rownames_to_column("AA") %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), prop.table)) %>%
        tibble::column_to_rownames("AA") %>%
        as.matrix()
    return(prop)
}

#' calc_AA_motif_zscore
#'
#' Calculate the Z score for the amino acids at each position
#'
#' @param background_count_matrix the count matrix from the
#' background sequences
#' @param background_prop_matrix the proportion matrix from the
#' background sequences
#' @param experimental_count_matrix the count matrix from the
#' experimental sequences
#' @param experimental_prop_matrix the proportion matrix from the
#' experimental sequences
#'
#' @return a data frame of Zscores for each amino acid at each position.
#' @keywords internal

calc_AA_motif_zscore <- function(background_count_matrix,
                                 background_prop_matrix,
                                 experimental_count_matrix,
                                 experimental_prop_matrix) {
    ### calculating the SD ###
    # defining function per Icelogo manual
    standard_deviation <- function(proportion_matrix, count_matrix) {
        (proportion_matrix / colSums(count_matrix))^(1 / 2)
    }

    # calculating the standard deviation
    bg_sd <- standard_deviation(
        background_prop_matrix,
        experimental_count_matrix
    )

    ### Calculating Z score ###

    # Defining function per icelogo manual
    zscore <- function(experimental_prop_matrix, background_prop_matrix, bg_sd) {
        (experimental_prop_matrix - background_prop_matrix) / bg_sd
    }

    # applying to our data
    zscores <- zscore(experimental_prop_matrix, background_prop_matrix, bg_sd) %>%
        as.data.frame()

    return(zscores)
}

#' calc_sig_zscores
#' Determine which Zscores are significant at the given alpha for a
#' matrix of scores
#' @param zscores = a data frame of zscores
#' @param pval = p value threshold for significance. Default is 0.05
#'
#' @return a tibble of significant zscores
#' @keywords internal

calc_sig_zscores <- function(zscores, pval = 0.05) {
    # no visible binding for global variable .
    . <- NULL
    # converting p value to zscore threshold. Divide by two since it is two tailed
    threshold <- stats::qnorm(p = pval / 2, lower.tail = FALSE)

    sig_zscores <- zscores %>%
        tibble::rownames_to_column("AA") %>%
        tidyr::pivot_longer(2:length(.),
            names_to = "position",
            values_to = "zscore"
        ) %>%
        dplyr::filter(abs(.data$zscore) > threshold) %>%
        dplyr::mutate(aa_position = paste0(.data$AA, ".", .data$position))

    return(sig_zscores)
}

#' calc_AA_percent_difference
#'
#' Calculate the percent difference between a matrix of background proportions
#' and a matrix of experimentally observed proportions.
#'
#' @param background_prop_matrix a proportion matrix of amino acids
#' per position from background cleavage sequences
#' @param experimental_prop_matrix a proportion matrix of amino acids per
#' position from experimental cleavage sequences
#'
#' @return a data frame of percent differences
#' @keywords internal
calc_AA_percent_difference <- function(background_prop_matrix,
                                       experimental_prop_matrix) {
    # Calculating percent difference
    pd <- (experimental_prop_matrix * 100) - (background_prop_matrix * 100) %>%
        as.data.frame()

    return(pd)
}

#' calc_AA_fc
#'
#' Calculate the fold change of each amino acid by position.
#'
#' @param experimental_prop_matrix a matrix of the experimental
#' proportions (from your vector of cleavage sequences) at each position.
#' @param background_prop_matrix a matrix of the background proportions
#'  of AAs at each position
#' @param sig_zscores a tibble of the significant zscores.
#' @return a matrix
#' @keywords internal
calc_AA_fc <- function(experimental_prop_matrix,
                       background_prop_matrix,
                       sig_zscores) {
    # calculating the fold change of each amino acid by position
    fold_change <- (experimental_prop_matrix * 100) /
        (background_prop_matrix * 100)
    # Filtering to only include significant z-scores
    prepared_fc <- prepare_fc(
        as.data.frame(fold_change),
        sig_zscores
    )
    converted_fc <- prepared_fc
    # Here we run into a problem with infinite values.
    # If a value is 0 in the count, but greater than that in the background
    # Icelogo shows infinite value as taking up the entire scale.
    max_sum <- max(purrr::map_df(as.data.frame(converted_fc), sum, na.rm = TRUE),
        na.rm = TRUE
    )
    final_converted_fc <- data.frame(row.names = rownames(converted_fc))
    for (i in seq_len(ncol(converted_fc))) {
        if (0 %in% converted_fc[, i]) {
            num_zero <- sum(converted_fc[, i] == 0, na.rm = TRUE)
            col <- replace(
                converted_fc[, i],
                converted_fc[, i] == 0,
                (max_sum / num_zero * -1)
            )
            final_converted_fc <- cbind(final_converted_fc, col)
        } else {
            final_converted_fc <- cbind(final_converted_fc, converted_fc[, i])
        }
    }
    colnames(final_converted_fc) <- colnames(converted_fc)
    # now converting values less than 1 to converted fold change.
    final_converted_fc2 <- final_converted_fc %>%
        dplyr::mutate(
            dplyr::across(
                dplyr::all_of(names(final_converted_fc)),
                ~ case_when(.x < 1 & .x > 0 ~ (1 / .x) * -1,
                    .default = .x
                )
            )
        )
    final_converted_fc2 <- as.matrix(final_converted_fc2)
    return(final_converted_fc2)
}
#' prepare_sig_p_dif
#'
#' Prepare significant percent difference data frame for iceLogo
#'
#' @param percent_difference  a data frame containing the percent differences
#' @param sig_zscores a matrix of significant amino acids at each position
#' based on z-scores
#'
#' @return a tibble
#' @keywords internal
prepare_sig_p_dif <- function(percent_difference, sig_zscores) {
    # No visible binding
    . <- NULL
    sig_final <- percent_difference %>%
        tibble::rownames_to_column("AA") %>%
        tidyr::pivot_longer(2:length(.),
            names_to = "position",
            values_to = "percent_difference"
        ) %>%
        dplyr::mutate(aa_position = paste0(.data$AA, ".", .data$position)) %>%
        dplyr::filter(.data$aa_position %in% sig_zscores$aa_position) %>%
        dplyr::select(-"aa_position") %>%
        tidyr::pivot_wider(
            names_from = "position",
            values_from = "percent_difference",
            names_sort = FALSE
        )
    # Dealing with the problem introduced by the pivot wider function
    missing_cols <- ncol(percent_difference) - (ncol(sig_final) - 1)
    missing_data <- as.data.frame(matrix(
        nrow = nrow(sig_final),
        ncol = missing_cols, NA
    ))
    `%!in%` <- Negate(`%in%`)
    missing_names <- names(percent_difference)[names(percent_difference) %!in%
        names(sig_final)]
    names(missing_data) <- missing_names
    final <- sig_final %>%
        dplyr::bind_cols(missing_data) %>%
        dplyr::relocate(paste0("P", c(
            (ncol(percent_difference) / 2):1,
            paste0(
                seq_len(ncol(percent_difference) / 2), "'"
            )
        ))) %>%
        tibble::column_to_rownames("AA") %>%
        as.matrix()
    return(final)
}

#' prepare_icelogo_data
#'
#' Prepare the final matrix containing iceLogo data for plotting.
#'
#' @param cleavage_seqs the cleavage sequences that are observed in the
#' experiment
#' @param background_universe a vector of the cleavage sequences to use as
#' the background.
#' @param pval the p-value threshold to consdier
#' @param type the type of iceLogo calculation to perform. Accepted values are
#' "percent_difference" or "fold_change".
#'
#' @return a matrix of enriched amino acids per position
#' @keywords internal
prepare_icelogo_data <- function(cleavage_seqs,
                                 background_universe =
                                     mspms::all_possible_8mers_from_228_library,
                                 pval = 0.05,
                                 type = "percent_difference") {
    # calculation proportions of background
    background_counts <- calc_AA_count_of_motif(background_universe)
    background_proportions <- calc_AA_prop_of_motif(background_counts)

    # calculate proportions of experimentally identified motifs
    experimental_counts <- calc_AA_count_of_motif(cleavage_seqs)
    experimental_proprotions <- calc_AA_prop_of_motif(experimental_counts)

    # calculate zscores
    zscores <- calc_AA_motif_zscore(
        background_count_matrix = background_counts,
        background_prop_matrix = background_proportions,
        experimental_count_matrix = experimental_counts,
        experimental_prop_matrix = experimental_proprotions
    )
    # calculate significant zscores
    sig_zscores <- calc_sig_zscores(zscores, pval)
    if (type == "percent_difference") {
        pd <- calc_AA_percent_difference(
            experimental_prop_matrix = experimental_proprotions,
            background_prop_matrix = background_proportions
        )
        final_pd <- prepare_sig_p_dif(pd, sig_zscores)
        return(final_pd)
    } else if (type == "fold_change") {
        fc <- calc_AA_fc(
            experimental_prop_matrix = experimental_proprotions,
            background_prop_matrix = background_proportions,
            sig_zscores = sig_zscores
        )
        return(fc)
    } else {
        stop("type must be either 'percent_difference' or 'fold_change'")
    }
}

#' prepare_fc
#'
#' Prepare fold changes of amino acids by position for Icelogo visualization.
#'
#' @param fold_change  a matrix of the fold changes of the AA by
#' position.
#' @param sig_zscores  a tibble of the significant zscores.
#'
#' @return a matrix of the fold changes of the significant AAs at each position.
#' @keywords internal

prepare_fc <- function(fold_change, sig_zscores) {
    # No visible binding for global variable .
    . <- NULL
    sig_final <- fold_change %>%
        tibble::rownames_to_column("AA") %>%
        tidyr::pivot_longer(2:length(.),
            names_to = "position",
            values_to = "fold_change"
        ) %>%
        dplyr::mutate(aa_position = paste0(.data$AA, ".", .data$position)) %>%
        dplyr::filter(.data$aa_position %in% sig_zscores$aa_position) %>%
        dplyr::select(-"aa_position") %>%
        tidyr::pivot_wider(
            names_from = "position",
            values_from = "fold_change",
            names_sort = FALSE
        )
    # Dealing with the problem introduced by the pivot wider function
    missing_cols <- ncol(fold_change) - (ncol(sig_final) - 1)
    missing_data <- as.data.frame(matrix(
        nrow = nrow(sig_final),
        ncol = missing_cols, NA
    ))
    `%!in%` <- Negate(`%in%`)
    missing_names <- names(fold_change)[names(fold_change) %!in%
        names(sig_final)]
    names(missing_data) <- missing_names
    final <- sig_final %>%
        dplyr::bind_cols(missing_data) %>%
        dplyr::relocate(paste0("P", c(
            (ncol(fold_change) / 2):1,
            paste0(
                seq_len(ncol(fold_change) / 2), "'"
            )
        ))) %>%
        tibble::column_to_rownames("AA") %>%
        as.matrix()
    return(final)
}
