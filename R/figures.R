# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' Creates and saves binned normalized Z-Score figure per region.
#'
#' @importFrom magrittr %>%
#' @param data A binned Z-scored data frame with columns \emph{focus}, \emph{start},
#'  \emph{PPDX_norm}, \emph{sample} and \emph{reference}.
#' @param sample_name sample name of interest
save_bin_plot <- function(samples, sample_name) {
  capt <-
    paste(
      "Supplementary bin figure complementing the final result",
      "figure. Each subregion depicts the bins from the region of interest",
      "and their deviation from the reference group."
    )
  
  ggplot2::ggsave(
    paste0(sample_name, ".bins", ".png"),
    ggplot2::ggplot(
      samples,
      ggplot2::aes(
        as.numeric(start),
        PPDX_norm,
        group = sample,
        color = reference
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ focus, scales = "free", ncol = 2) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::scale_x_continuous(
        labels = function(x)
          format(
            x,
            big.mark = " ",
            decimal.mark = ",",
            scientific = FALSE
          )
      ) +
      ggplot2::scale_color_manual(values = c("red", "grey")) +
      ggplot2::ylab("Normalized Z-score") +
      ggplot2::labs(
        title = paste(
          sample_name,
          "scientific",
          utils::packageName(),
          "bin report"
        ),
        subtitle = paste("Version: ", utils::packageVersion(utils::packageName())),
        caption = stringr::str_wrap(capt, width = 100)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(),
        strip.background = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
        plot.caption = ggplot2::element_text(
          hjust = 0,
          vjust = 0,
          face = "italic",
          color = "black"
        )
      )
    ,
    width = 10,
    height = max(1, round(length(
      unique(samples$focus)
    ) / 2)) * 2
  )
}

