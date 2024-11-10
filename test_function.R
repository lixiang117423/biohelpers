###################################################
############# pca_in_one ##########################
###################################################

library(tidyverse)

df = iris[,1:4]
sample = iris$Species %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  magrittr::set_names(c("sample", "species"))


FactoMineR::PCA(df, ncp = 5, graph = FALSE) -> pca.res

pca.res %>%
  factoextra::get_eigenvalue() %>%
  base::as.data.frame() %>%
  dplyr::mutate(pc = base::paste0("PC", 1:nrow(.))) %>%
  dplyr::arrange(-variance.percent) %>%
  dplyr::mutate(pc = base::factor(pc, levels = base::unique(pc))) %>%
  dplyr::mutate(variance.percent = base::round(variance.percent, 2)) %>%
  dplyr::select(pc, variance.percent) -> eig.val


pca.res[["ind"]][["coord"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample")  %>%
  dplyr::left_join(sample,  by = "sample") %>%
  ggplot(aes(Dim.1, Dim.2, color = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size = 3) +
  labs(
    x = paste0("PC1 (", eig.val$variance.percent[1], "%)"),
    y = paste0("PC2 (", eig.val$variance.percent[2], "%)")
  ) -> plot.pca


# 测试函数
devtools::load_all()

data <- iris[, 1:4]
sample <- iris$Species %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  magrittr::set_names(c("sample", "species"))

pca_in_one(data, sample, x = "pc1", y = "pc3") -> result.pca

result.pca$plot

# 相关性函数
df = iris[,1:4]

WGCNA::corAndPvalue(df) -> cor.res

cor.res$cor %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "from") %>% 
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "to", values_to = "cor") -> df.cor

cor.res$p %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "from") %>% 
  tidyr::pivot_longer(cols = 2:ncol(.), names_to = "to", values_to = "pvalue") -> df.pvalue

dplyr::left_join(df.cor, df.pvalue)












