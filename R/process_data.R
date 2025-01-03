# # 示例数据
# test_data = iris
# usethis::use_data(test_data, overwrite = TRUE)

# # reorder2heatmap
# readr::read_delim("raw.data/reorder2heatmap.txt") %>%
#   magrittr::set_names(c("sample", "meta", "value")) %>%
#   dplyr::mutate(sample = as.character(sample))-> df.reorder2heatmap

# usethis::use_data(df.reorder2heatmap, overwrite = TRUE)

# # call_DEGs_DESeq2
# readxl::read_excel("D://OneDrive/NAS/科研相关/PhData/data/03.生信挖掘/07.水稻PikPia/data/水稻转录组信息收集.xlsx", sheet = "final") %>%
#   dplyr::filter(BioProject == "PRJNA661210") %>%
#   dplyr::mutate(LibraryLayout = "SINGLE") %>%
#   dplyr::filter(Treatment.last == "12h") -> df.sample

# readr::read_delim("D://OneDrive/NAS/科研相关/PhData/data/03.生信挖掘/08.作物多效基因/data/all.expression.table.txt") %>%
#   dplyr::filter(Run %in% df.sample$Run) -> all.expre

# df.sample %>%
#   dplyr::select(Run, Treatment.minor) %>%
#   magrittr::set_names(c("sample", "group")) %>%
#   dplyr::mutate(group = factor(group, level = c("CK","Guy11"))) %>%
#   tibble::column_to_rownames(var = "sample") -> df.rnaseq.sample

# all.expre %>%
#   dplyr::select(gene, FPKM, Run) %>%
#   dplyr::mutate(FPKM = round(FPKM, 0)) %>%
#   dplyr::distinct_all() %>%
#   dplyr::mutate(tmp = paste0(gene, Run)) %>%
#   dplyr::filter(!duplicated(tmp)) %>%
#   dplyr::select(-tmp) %>%
#   tidyr::pivot_wider(names_from = Run, values_from = FPKM) %>%
#   tibble::column_to_rownames(var = "gene") %>%
#   dplyr::select(rownames(df.rnaseq.sample))-> df.rnaseq.gene

# usethis::use_data(df.rnaseq.gene, overwrite = TRUE)
# usethis::use_data(df.rnaseq.sample, overwrite = TRUE)

# # enrich_GO
# load("D:/OneDrive/NAS/科研相关/PhData/data/02.元阳梯田/1569转录+代谢+激素文章/data/acuce.go.RData")

# df.go %>%
#   dplyr::select(1,3:5) %>%
#   magrittr::set_names(c("gene", "go.id", "go.term", "go.ontology")) -> df.rnaseq.go

# usethis::use_data(df.rnaseq.go, overwrite = TRUE)

# # 随机选择500个基因作为差异表达基因
# set.seed(123)
# df.rnaseq.go %>%
#   dplyr::select(gene) %>%
#   dplyr::sample_n(500, replace = FALSE) -> df.rnaseq.degs

# usethis::use_data(df.rnaseq.degs, overwrite = TRUE)

# # enrich_KEGG
# readxl::read_excel("D:/OneDrive/NAS/科研相关/PhData/data/99.其他项目/GOandKEGG/kegg.20220727.xlsx") %>%
#   dplyr::select(geneid, keggid, description) %>%
#   magrittr::set_names(c("gene", "kegg.id", "kegg.term")) -> df.rnaseq.kegg

# usethis::use_data(df.rnaseq.kegg, overwrite = TRUE)

# # plot_volcano
# library(dplyr)
# library(biohelpers)

# data(df.rnaseq.gene)
# data(df.rnaseq.sample)
# call_DEGs_DESeq2(data = df.rnaseq.gene,
#                  sample = df.rnaseq.sample,
#                  group = "group") -> df.rnaseq.plot_volcano

# usethis::use_data(df.rnaseq.plot_volcano, overwrite = TRUE)

# # PCoA
# readxl::read_excel("D:/OneDrive/NAS/科研相关/PhData/data/03.生信挖掘/10.松针腐解/data/微生物组数据.最终数据.xlsx", sheet = "otu") %>%
#   dplyr::filter(environment == "Cropland", soil != "non-continuous cropping soil", group == "bacteria") -> df.otu

# df.otu %>%
#   dplyr::select(treatment, replicates, OTU, value) %>%
#   dplyr::mutate(sample = paste0(treatment, "_", replicates)) %>%
#   dplyr::select(sample, OTU, value) %>%
#   tidyr::pivot_wider(names_from = OTU, values_from = value) %>%
#   tibble::column_to_rownames(var = "sample") -> df.pcoa.otu

# usethis::use_data(df.pcoa.otu, overwrite = TRUE)

# rownames(df.pcoa.otu) %>%
#   as.data.frame() %>%
#   magrittr::set_names("sample") %>%
#   dplyr::mutate(group = stringr::str_split(sample, "_") %>% sapply("[", 1)) -> df.pcoa.sample

# usethis::use_data(df.pcoa.sample, overwrite = TRUE)

# # top_10
# readxl::read_excel("D:/OneDrive/NAS/科研相关/PhData/data/03.生信挖掘/10.松针腐解/data/微生物组数据.最终数据.xlsx", sheet = "otu") -> df.otu

# df.otu %>%
#   dplyr::filter(environment == "Cropland", soil != "non-continuous cropping soil", group == "bacteria") %>%
#   dplyr::select(treatment, replicates, OTU, value, phylum) %>%
#   dplyr::mutate(sample = paste0(treatment, "_", replicates)) %>%
#   dplyr::select(sample, OTU, value) %>%
#   tidyr::pivot_wider(names_from = OTU, values_from = value) %>%
#   tibble::column_to_rownames(var = "sample") -> df.top10.otu

# usethis::use_data(df.top10.otu, overwrite = TRUE)

# df.top10.otu %>%
#   rownames() %>%
#   as.data.frame() %>%
#   magrittr::set_names("sample") %>%
#   dplyr::mutate(group = stringr::str_split(sample, "_") %>% sapply("[", 1)) -> df.top10.sample

# usethis::use_data(df.top10.sample, overwrite = TRUE)

# df.otu %>%
#   dplyr::filter(environment == "Cropland", soil != "non-continuous cropping soil", group == "bacteria") %>%
#   dplyr::select(OTU, phylum, class, order, family, genus, species) %>%
#   dplyr::distinct_all() -> df.top10.class

# usethis::use_data(df.top10.class, overwrite = TRUE)

# call_DAMs_LEfSe
# readxl::read_excel("D:/OneDrive/NAS/科研相关/PhData/data/03.生信挖掘/10.松针腐解/data/微生物组数据.最终数据.xlsx", sheet = "otu") %>%
#   dplyr::filter(environment == "Cropland", soil != "non-continuous cropping soil", group == "bacteria", treatment %in% c("CK", "Cover")) -> df.otu

# df.otu %>%
#   dplyr::select(treatment, OTU, value) %>%
#   dplyr::group_by(treatment, OTU) %>%
#   dplyr::summarise(sd = sd(value)) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(sd != 0) -> df.filtered.otu

# df.otu %>%
#   dplyr::select(treatment, replicates, OTU, value) %>%
#   dplyr::mutate(sample = paste0(treatment, replicates)) %>%
#   dplyr::select(sample, OTU, value) %>%
#   tidyr::pivot_wider(names_from = sample, values_from = value) %>%
#   dplyr::filter(OTU %in% df.filtered.otu$OTU) %>%
#   tibble::column_to_rownames(var = "OTU") -> df.call_DAMs_LEfSe.otu

# usethis::use_data(df.call_DAMs_LEfSe.otu, overwrite = TRUE)

# df.otu %>%
#   dplyr::select(treatment, replicates) %>%
#   dplyr::mutate(sample = paste0(treatment, replicates)) %>%
#   dplyr::select(sample, treatment) %>%
#   dplyr::rename(group = treatment) %>%
#   dplyr::distinct_all() %>%
#   tibble::column_to_rownames(var = "sample") -> df.call_DAMs_LEfSe.sample

# usethis::use_data(df.call_DAMs_LEfSe.sample, overwrite = TRUE)

# SummarizedExperiment(assays = list(counts = otu), colData = sample) %>%
#   lefser::lefser(groupCol = "group", kruskal.threshold = 0.05, wilcox.threshold = 0.05, lda.threshold = 1) -> lefse.result
