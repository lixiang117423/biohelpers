# 示例数据
test_data = iris
usethis::use_data(test_data, overwrite = TRUE)

# reorder2heatmap
readr::read_delim("raw.data/reorder2heatmap.txt") %>% 
  magrittr::set_names(c("sample", "meta", "value")) %>% 
  dplyr::mutate(sample = as.character(sample))-> df.reorder2heatmap

usethis::use_data(df.reorder2heatmap, overwrite = TRUE)
