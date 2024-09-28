#' plot_icelogo
#'
#' This function plots the cleavage motifs that were enriched relative to
#' background as implemented in the iceLogo method.
#' https://iomics.ugent.be/icelogoserver/resources/manual.pdf
#'
#' @param cleavage_seqs  these are the cleavage sequences of interest
#' @param background_universe  this is a list of cleavage sequences
#' to use as the background in building the iceLogo.
#' @param pval  this is the pvalue threshold (<=) to consider significant when
#' determining the significance of the sig_cleavages relative to the background
#' at each position of the iceLogo.
#' @param type this is the type of visualization you would like to perform,
#' accepted values are  either "percent_difference" or "fold_change".
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' # Determining significant cleavages for catA
#' catA_sig_cleavages <- mspms::log2fc_t_test_data %>%
#'     dplyr::filter(p.adj <= 0.05, log2fc > 3) %>%
#'     dplyr::filter(condition == "CatA") %>%
#'     dplyr::pull(cleavage_seq) %>%
#'     unique()
#'
#' # Plotting icelogo
#' plot_icelogo(catA_sig_cleavages,
#'     background_universe = all_possible_8mers_from_228_library
#' )
plot_icelogo <- function(cleavage_seqs,
                         background_universe =
                             mspms::all_possible_8mers_from_228_library,
                         pval = 0.05,
                         type = "percent_difference") {
    # calculate n
    n <- length(cleavage_seqs)
    # processing data to plot
    icelogo_data <- prepare_icelogo_data(
        cleavage_seqs,
        background_universe,
        pval,
        type
    )
    nchar_cleav <- ncol(icelogo_data)
    # Creating color scheme for mspms data
    col_scheme <- icelogo_col_scheme()
    # Plotting the motif
    p1 <- ggseqlogo::ggseqlogo(icelogo_data,
        font = "helvetica_light",
        method = "custom",
        seq_type = "AA",
        col_scheme = col_scheme
    ) +
        ggplot2::scale_x_continuous(
            labels = paste0("P", c(
                (nchar_cleav / 2):1,
                paste0(
                    seq_len(nchar_cleav / 2), "'"
                )
            )),
            breaks = seq_len(nchar_cleav)
        ) +
        ggplot2::labs(subtitle = paste0("n = ", n)) +
        if (type == "percent_difference") {
            ggplot2::ylab("percent difference")
        } else if (type == "fold_change") {
            ggplot2::ylab("fold change")
        } else {
            stop("type must be either percent_difference or fold_change")
        }
    return(p1)
}
