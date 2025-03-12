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
          read_csv(paste(lines, collapse = "\n"), col_names = str_detect(n, "_Data"), show_col_types = FALSE)
        }
      )
    )
  # return a list of tibbles
  ret <- blocks_tib$table
  names(ret) <- blocks_tib$names
  ret
}



#' get the fastq paths relative to the current working directory
#'
#' You should run this from the directory that will be your data_parent_dir
#' for any mega-non-model runs you are going to do.
#'
#' This returns the path to the file that the file listing is in.
#' @param fastq_dir the relative path to the folder that holds the fastqs named xxx.fastq.gz
get_fastq_paths <- function(fastq_dir) {
  tfile <- tempfile()
  CALL <- paste("ls -l ", fastq_dir, "/*.fastq.gz > ", tfile, collapse = "", sep = "")
  system(CALL)
  tfile  # return the tempfile path
}



#' Create the units.csv file for a run directory and write it back out to that directory
#'
#' @param fastq_dir the relative path to the folder that holds the fastqs named xxx.fastq.gz
#'
#' You should run this from the directory that will be your data_parent_dir
#' for any mega-non-model runs you are going to do.
#'
#' This returns the path to that units file
write_units_file <- function(fastq_dir) {

  # get the file listing and reduce it to paths and kb
  tf <- get_fastq_paths(fastq_dir)
  file_listing <- read_table(tf, col_names = FALSE, show_col_types = FALSE)
  f2 <- file_listing %>%
    rename(fq = X9) %>%
    mutate(kb = as.numeric(X5) / 1024) %>%
    select(fq, kb)

  # get the sample sheet.  This assumes it is named SampleSheet.csv and is stored
  # in the fastq_dir
  ssList <- slurp_sample_sheet(file.path(fastq_dir, "SampleSheet.csv"))

  # now we extract sample ids, etc. from the fastq names and join them and pivot, etc.
  f2_wide <- f2 %>%
    mutate(base = basename(fq)) %>%
    extract(
      base,
      into = c("sample", "read"),
      regex = "^(.*)_S[0-9]+_R([12])_001\\.fastq\\.gz",
      convert = TRUE,
      remove = FALSE
    ) %>%
    select(-base) %>%
    pivot_wider(names_from = read, values_from = c(kb, fq), names_sep = "") %>%
    mutate(
      lane = 1L,  # these are all assumed to be on a single lane (for now),
      sample_id = sample,
      platform = "ILLUMINA"
    )

  # and now we can join library and barcode information. We do a full join
  # so we can warn about paths missing from the SampleSheet (and vice versa).
  joined <- full_join(f2_wide, ssList$Cloud_Data, by = join_by(sample == Sample_ID)) %>%
    rename(
      flowcell = ProjectName,
      library = LibraryName
    ) %>%
    full_join(ssList$BCLConvert_Data, by = join_by(sample == Sample_ID)) %>%
    mutate(
      barcode = paste(Index, Index2, sep = "+"),
    )


  # now we do a little bit of error checking, catching situations where the SampleSheet
  # had a sample that we don't have a path for, and also where we have paths that
  # were not in the sample sheet.
  ready <- joined %>%
    filter(!is.na(fq1) & !is.na(Index)) %>%
    select("sample", "library", "flowcell", "platform",	"lane", "sample_id", "barcode", "fq1", "fq2", "kb1",	"kb2") %>%
    group_by(sample) %>%
    mutate(unit = 1:n(), .after = sample) %>%
    ungroup()

  missing_paths <-  joined %>%
    filter(is.na(fq1))
  missing_barcodes <- joined %>%
    filter(is.na(Index))

  if(nrow(ready) == 0) {
    stop("No fastq paths ready for writing to units file. Something has gone horribly awry.")
  }
  if(nrow(missing_paths) > 0) {
    warning("Not all samples in SampleSheet.csv have fastq paths.  Check return tibble $missing_paths for info.", immediate. = TRUE)
  }
  if(nrow(missing_barcodes) > 0) {
    warning("Not all fastq paths have corresponding entries in SampleSheet.csv.  Check return tibble $missing_barcodes for info.", immediate. = TRUE)
  }

  write_tsv(ready, file = file.path(fastq_dir, "units.tsv"))

  # return the different tibbles in a list
  list(
    ready = ready,
    missing_paths = missing_paths,
    missing_barcodes = missing_barcodes
  )
}
