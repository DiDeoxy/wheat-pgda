#' Loads gene info
#'
#' Loads the info on the selected gene alignments
#'
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>%
#' @importFrom readr read_csv
#'
#' @param csv the csv file containg the alignments
#' @param base used further on as the y axis position of the genes during
#' plotting
#'
#' @return a tibble contianign the parsed info
#'
#' @export

load_genes <- function(csv, base = 0) {
  import::from(dplyr, "mutate", "select")
  import::from(magrittr, "%>%")
  import::from(readr, "read_csv")

  read_csv(
    csv, col_names = c("id", "chrom", "pos", "sleng", "salign", "%id")
  ) %>%
    mutate(pos_mb = pos / 1e6) %>%
    select(id, chrom, pos_mb) %>%
    cbind(base = base)
}