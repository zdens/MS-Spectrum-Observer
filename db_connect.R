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
  select_spec_files <- data.frame(filename = list.files(root_path, recursive = TRUE, pattern = ".rds"))
  select_spec_files$id <- seq_len(nrow(select_spec_files))
  return(select_spec_files)
}