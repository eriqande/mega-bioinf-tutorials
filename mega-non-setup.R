# this is a few functions to make it easy to pull data off the NextSeq
# and prep it for doing mega-non-model runs, and light prelim QC

require(tidyverse)



#' Read blocked CSV (like in the sample sheets off the NextSeq) into a named list of tibbles
#'
#' @param path  The path of the SampleSheet.csv file
slurp_sample_sheet <- function(path) {
  lvec <- read_lines(path)
  blocks_tib <- tibble(
    start_positions = which(str_detect(lvec, "^\\[[a-zA-Z_]+\\]"))
  ) %>%
    mutate(
      names = str_replace_all(lvec[start_positions], "[^a-zA-Z_]", ""),
      end_positions = lead(start_positions, default = length(lvec)),
      table = pmap(
        .l = list(s = start_positions, e = end_positions, n = names),
        .f = function(s, e, n) {
          lines <- lvec[(s + 1):(e - 1)]
          read_csv(paste(lines, collapse = "\n"), col_names = str_detect(n, "_Data"))
        }
      )
    )
  # return a list of tibbles
  ret <- blocks_tib$table
  names(ret) <- blocks_tib$names
  ret
}


