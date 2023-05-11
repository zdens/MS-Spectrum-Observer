
parse_filename <- function(filenames) {
  fn_splitted <- stri_split_fixed(filenames, "/")
  dev_list <- sapply(fn_splitted, function(x){
    ifelse(x[[1]] == "1", "1", "2")
  })
  fn_list <- sapply(fn_splitted, function(x){
    x[[2]]
  })
  return(list(fi = stri_join(dev_list, fn_list, sep = ", "),
              protocol_name = sapply(fn_splitted, function(x){
                x[[1]]
              })))
}

get_files_from_db <- function(root_path) {
  select_spec_files <- data.frame(filename = list.files(root_path, recursive = TRUE, pattern = ".cdf"))
  select_spec_files$id <- seq_len(nrow(select_spec_files))
  file_prot_info <- parse_filename(select_spec_files$filename)
  select_spec_files$fi <- file_prot_info$fi
  select_spec_files$protocol <- file_prot_info$protocol_name
  return(select_spec_files)
}