# 检查R包名称的可用性，外网可用
available::available("biohelpers")

# 添加MIT许可证信息
usethis::use_mit_license()

# 添加文献引用信息
# usethis::use_citation()

# 添加开发记录文件
# file.create("deve.log.R")

# 添加github链接和bug提交链接，然后手动编辑DESCRIPTION这个文件中对应的内容
usethis::use_github_links()

# 添加示例数据,直接存储为.rda格式
# file.create("./R/data.R") # 解释每个示例数据
# file.create("./R/process_data.R") # 数据处理记录脚本

# 将编译时用不到的文件添加到.Rbuildignore这个文件夹中，不会改变文件位置，只是在编译时忽略这些文件
usethis::use_build_ignore("README.md")
usethis::use_build_ignore("deve.log.R")
usethis::use_build_ignore("test_function.R")
usethis::use_build_ignore("biohelpers.Rproj")
usethis::use_build_ignore("R/process_data.R") # 编译时忽略这个脚本
usethis::use_git_ignore("raw.data/")
usethis::use_build_ignore("raw.data/")
usethis::use_build_ignore("R/process_data.R")  
usethis::use_git_ignore("R/process_data.R")

# 添加全局变量，所有的全局变量全部添加在这个脚本里面
# file.create("R/biohelpers-global.R")

# 添加函数,添加一次即可，避免覆盖
# file.create("R/pca_in_one.R")
# file.create("R/cor_and_plot.R")
# file.create("R/lm_and_plot.R")
# file.create("R/find_outliner.R")
# file.create("R/plot_theme.R")
# file.create("R/reorder2heatmap.R")
# file.create("R/show_minor_tricks.R")
# file.create("R/call_DEGs_DESeq2.R")
# file.create("R/enrich_GO.R")
# file.create("R/enrich_KEGG.R")
# file.create("R/plot_volcano.R")
# file.create("R/PCoA.R")
# file.create("R/top_10.R")
# file.create("R/call_DEMs_DESeq2.R")
# file.create("R/call_DAMs_LEfSe.R")
# file.create("R/Microbiome.RDA.R")
# file.create("R/CommonlyUsed.df_to_list.R")
# file.create("R/CommonlyUsed.anova_and_post_test.R")
# file.create("R/Microbiome.permanova.R")
# file.create("R/Metabolome.splsda.R")
# file.create("R/Metabolome.oplsda.R")
# file.create("R/Microbiome.rrarefy.R")
# file.create("R/PopGen.manhattan_plot.R")


# 修改版本
usethis::use_version("major") # 第一位数字，当你做了不向后兼容的 API 修改时，增加主版本号。
usethis::use_version("minor") # 第二位数字，当你在不破坏向后兼容性的前提下增加功能时，增加次版本号。
usethis::use_version("patch") # 第三位数字，当你做了向后兼容的问题修正时，增加修订号。
usethis::use_version("dev") # 第四位数字，通常用来表示开发中的版本，可能不稳定或不完整。

# 添加用到的R包和函数
usethis::use_import_from("dplyr","%>%")
usethis::use_import_from("dplyr","arrange")
usethis::use_import_from("dplyr","mutate")
usethis::use_import_from("dplyr","select")
usethis::use_import_from("dplyr","left_join")
usethis::use_import_from("FactoMineR","PCA")
usethis::use_import_from("factoextra","get_eigenvalue")
usethis::use_import_from("tibble","rownames_to_column")
usethis::use_import_from("ggplot2", "ggplot")
usethis::use_import_from("ggplot2", "aes")
usethis::use_import_from("ggplot2", "geom_vline")
usethis::use_import_from("ggplot2", "geom_hline")
usethis::use_import_from("ggplot2", "geom_point")
usethis::use_import_from("ggplot2", "labs")
usethis::use_import_from("magrittr", "set_names")
usethis::use_import_from("stringr", "str_replace")
usethis::use_import_from("WGCNA","corAndPvalue")
usethis::use_import_from("tidyr", "pivot_longer")
usethis::use_import_from("tidyr", "nest")
usethis::use_import_from("purrr", "map")
usethis::use_import_from("purrr", "map_dbl")
usethis::use_import_from("ggplot2", "geom_point")
usethis::use_import_from("ggplot2", "geom_smooth")
usethis::use_import_from("stats", "anova")
usethis::use_import_from("rlang", "sym")
usethis::use_import_from("stats", "quantile")
usethis::use_import_from("stats", "IQR")
usethis::use_import_from("ggthemes", "theme_foundation")
usethis::use_import_from("ggplot2", "theme")
usethis::use_import_from("ggplot2", "element_line")
usethis::use_import_from("ggplot2", "element_rect")
usethis::use_import_from("ggplot2", "element_text")
usethis::use_import_from("ggplot2", "element_blank")
usethis::use_import_from("ggplot2", "unit")
usethis::use_import_from("ggplot2", "rel")
usethis::use_import_from("dplyr","slice_max")
usethis::use_import_from("dplyr","rename")
usethis::use_import_from("dplyr","count")
usethis::use_import_from("dplyr","arrange")
usethis::use_import_from("dplyr","inner_join")
usethis::use_import_from("dplyr","n")
usethis::use_import_from("stats", "reorder")
usethis::use_import_from("ggprism", "guide_prism_offset_minor")
usethis::use_import_from("ggplot2", "guides")
usethis::use_import_from("DESeq2", "DESeqDataSetFromMatrix")
usethis::use_import_from("DESeq2", "DESeq")
usethis::use_import_from("DESeq2", "results")
usethis::use_import_from("tibble", "rownames_to_column")
usethis::use_import_from("stringr", "str_split")
usethis::use_import_from("clusterProfiler", "enricher")
usethis::use_import_from("ggplot2", "geom_hline")
usethis::use_import_from("ggplot2", "geom_vline")
usethis::use_import_from("ggsci", "scale_color_d3")
usethis::use_import_from("vegan", "vegdist")
usethis::use_import_from("ape", "pcoa")
usethis::use_import_from("rstatix", "wilcox_test")
usethis::use_import_from("scales", "percent")
usethis::use_import_from("ggplot2", "position_fill")
usethis::use_import_from("SummarizedExperiment", "SummarizedExperiment")
usethis::use_import_from("lefser", "lefser")
usethis::use_import_from("vegan", "decostand")
usethis::use_import_from("vegan", "rda")
usethis::use_import_from("vegan", "scores")
usethis::use_import_from("vegan", "envfit")
usethis::use_import_from("ggrepel", "geom_label_repel")
usethis::use_import_from("ggplot2", "arrow")
usethis::use_import_from("ggplot2", "unit")
usethis::use_import_from("ggplot2", "geom_segment")
usethis::use_import_from("multcomp","glht")
usethis::use_import_from("multcomp","mcp")
usethis::use_import_from("multcomp","cld")
usethis::use_import_from("broom","tidy")
usethis::use_import_from("agricolae","duncan.test")
usethis::use_import_from("stats", "aov")
usethis::use_import_from("vegan", "adonis")
usethis::use_import_from("mixOmics", "splsda")
usethis::use_import_from("ropls", "opls")
usethis::use_import_from("chemhelper", "get_scores")
usethis::use_import_from("dplyr", "lag")
usethis::use_import_from("janitor", "clean_names")
usethis::use_import_from("normentR", "simulateGWAS")
usethis::use_import_from("ggtext", "element_markdown")
usethis::use_import_from("forcats", "as_factor")


123

# 将某些文件格式化为tidyverse风格
dir("R") %>% 
  as.data.frame() %>% 
  magrittr::set_names("file") %>% 
  dplyr::mutate(path = paste0("R/", file)) %>% 
  dplyr::pull(path) %>%  
  purrr::walk(~ styler::style_file(.))

123

# 编译vignettes
# usethis::use_vignette(name = "pca_in_one") # 运行第二次会覆盖之前的
# usethis::use_vignette(name = "cor_and_plot")
# usethis::use_vignette(name = "lm_and_plot")
# usethis::use_vignette(name = "ind_outliner")
# usethis::use_vignette(name = "plot_theme")
# usethis::use_vignette(name = "reorder2heatmap")

# devtools::build_vignettes()

# 检查
file.remove("./NAMESPACE")
devtools::load_all()
devtools::document()
usethis::use_tidy_description()
devtools::check()

# 编译R包并安装
# devtools::build()
# devtools::check_built("../biohelpers_0.0.0.5.tar.gz")
devtools::install_local(force = TRUE)
# usethis::use_github_release(publish = TRUE)

# pkgdown添加网页,后续推送会自动更新
# usethis::use_pkgdown_github_pages()
# pkgdown::build_site(new_process = FALSE)
