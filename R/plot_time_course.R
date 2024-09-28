#' plot_time_course
#'
#' Easily plot a time course of all peptides in a QFeatures object by
#' peptide.
#'
#'
#' @param mspms_tidy_data tidy mspms data (prepared from QFeatures object
#' by mspms_tidy())
#' @param value_colname the name of the column containing values.
#' @param summarize_by_mean whether to summarise by mean (TRUE- show error bars
#' +- 1 standard deviation) or not (FALSE)
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' # Determining peptide of interest
#' max_log2fc_pep <- mspms::log2fc_t_test_data %>%
#'     dplyr::filter(p.adj <= 0.05, log2fc > 3) %>%
#'     dplyr::filter(log2fc == max(log2fc)) %>%
#'     dplyr::pull(peptide)
#'
#' # Defining QFeatures filter
#' filtered <- mspms::mspms_tidy_data %>%
#'     dplyr::filter(peptide == max_log2fc_pep) %>%
#'     plot_time_course()
plot_time_course <- function(mspms_tidy_data,
                             value_colname = "peptides_norm",
                             summarize_by_mean = FALSE) {
    value_colname <- dplyr::sym(value_colname)
    if (isTRUE(summarize_by_mean)) {
        p1 <- mspms_tidy_data %>%
            dplyr::group_by(
                .data$peptide, .data$condition, .data$time,
                .data$cleavage_seq
            ) %>%
            dplyr::summarize(
                mean = mean(!!value_colname, na.rm = TRUE),
                sd = stats::sd(!!value_colname, na.rm = TRUE)
            ) %>%
            ggplot2::ggplot(ggplot2::aes(
                x = .data$time,
                y = .data$mean,
                color = .data$condition
            )) +
            ggplot2::geom_point() +
            ggplot2::geom_line() +
            ggplot2::geom_errorbar(ggplot2::aes(
                ymax = .data$mean + .data$sd,
                ymin = .data$mean - .data$sd
            ), width = 15) +
            ggplot2::facet_wrap(~ peptide + cleavage_seq, scales = "free_y")
    } else {
        p1 <- mspms_tidy_data %>%
            ggplot2::ggplot(ggplot2::aes(
                x = .data$time,
                y = !!value_colname,
                color = .data$condition
            )) +
            ggplot2::geom_point() +
            ggplot2::stat_summary(ggplot2::aes(group = .data$condition),
                fun = mean,
                geom = "line"
            ) +
            ggplot2::facet_wrap(~ peptide + cleavage_seq, scales = "free_y")
    }
    return(p1)
}
