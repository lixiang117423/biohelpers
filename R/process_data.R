# 示例数据
test_data = iris
usethis::use_data(test_data, overwrite = TRUE)

# reorder2heatmap
readr::read_delim("raw.data/reorder2heatmap.txt") %>% 
  magrittr::set_names(c("sample", "meta", "value")) %>% 
  dplyr::mutate(sample = as.character(sample))-> df.reorder2heatmap

usethis::use_data(df.reorder2heatmap, overwrite = TRUE)

# call_DEGs_DESeq2
readxl::read_excel("D://OneDrive/NAS/科研相关/PhData/data/03.生信挖掘/07.水稻PikPia/data/水稻转录组信息收集.xlsx", sheet = "final") %>%
  dplyr::filter(BioProject == "PRJNA661210") %>%
  dplyr::mutate(LibraryLayout = "SINGLE") %>%
  dplyr::filter(Treatment.last == "12h") -> df.sample

readr::read_delim("D://OneDrive/NAS/科研相关/PhData/data/03.生信挖掘/08.作物多效基因/data/all.expression.table.txt") %>%
  dplyr::filter(Run %in% df.sample$Run) -> all.expre

df.sample %>% 
  dplyr::select(Run, Treatment.minor) %>% 
  magrittr::set_names(c("sample", "group")) %>% 
  dplyr::mutate(group = factor(group, level = c("CK","Guy11"))) %>% 
  tibble::column_to_rownames(var = "sample") -> df.sample.rnaseq

all.expre %>% 
  dplyr::select(gene, FPKM, Run) %>% 
  dplyr::mutate(FPKM = round(FPKM, 0)) %>% 
  dplyr::distinct_all() %>% 
  dplyr::mutate(tmp = paste0(gene, Run)) %>% 
  dplyr::filter(!duplicated(tmp)) %>% 
  dplyr::select(-tmp) %>% 
  tidyr::pivot_wider(names_from = Run, values_from = FPKM) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  dplyr::select(rownames(df.sample.rnaseq))-> df.gene.rnaseq

usethis::use_data(df.gene.rnaseq, overwrite = TRUE)
usethis::use_data(df.sample.rnaseq, overwrite = TRUE)
