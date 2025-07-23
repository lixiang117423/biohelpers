# biohelpers

<!-- badges: start -->
[![R-CMD-check](https://github.com/lixiang117423/biohelpers/workflows/R-CMD-check/badge.svg)](https://github.com/lixiang117423/biohelpers/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/biohelpers)](https://CRAN.R-project.org/package=biohelpers)
<!-- badges: end -->

## 简介

biohelpers 是一个为生物数据处理提供便利函数的R包。该包简化了常见的生物数据分析流程，提供标准化的分析和可视化输出。

## 安装

你可以从GitHub安装开发版本：

```r
# 安装开发版本
devtools::install_github("lixiang117423/biohelpers")
```

## 主要功能

### 通用分析函数
- `pca_analysis()`: 主成分分析 (PCA)
- `cor_analysis()`: 相关性分析和可视化
- `lm_analysis()`: 线性回归分析和可视化
- `find_outliers()`: 异常值检测
- `theme_bio()`: 生物数据可视化主题
- `reorder_heatmap()`: 热图数据重排序

### 代谢组学分析
- `opls_analysis()`: OPLS-DA分析
- `spls_analysis()`: sPLS-DA分析

### 微生物组学分析
- `pcoa_analysis()`: 主坐标分析 (PCoA)
- `rda_analysis()`: 冗余分析 (RDA)
- `find_dams_deseq2()`: DESeq2差异分析
- `find_dams_lefse()`: LEfSe差异分析
- `permanova_test()`: PERMANOVA置换检验
- `rarefy_table()`: 稀释抽样
- `top_taxa()`: 获取主要物种

### 转录组学分析
- `find_degs_deseq2()`: DESeq2差异表达分析
- `enrich_go()`: GO功能富集分析
- `enrich_kegg()`: KEGG通路富集分析
- `volcano_plot()`: 火山图绘制

### 群体遗传学分析
- `manhattan_plot()`: 曼哈顿图绘制

## 使用示例

### PCA分析
```r
library(biohelpers)

# 准备数据
data <- iris[,1:4]
sample_info <- data.frame(
  sample = paste0("sample", 1:150),
  species = iris$Species
)

# 进行PCA分析
pca_result <- pca_analysis(data = data, sample = sample_info)
```

### 相关性分析
```r
# 相关性分析
cor_result <- cor_analysis(data.1 = iris[,1:2], data.2 = iris[,3:4])
```

### 线性回归分析
```r
# 线性回归分析
lm_result <- lm_analysis(
  data = iris,
  x = "Sepal.Length",
  y = "Sepal.Width", 
  group = "Species",
  color = "Species"
)
```

## 贡献

欢迎提出问题和贡献代码！请在[GitHub Issues](https://github.com/lixiang117423/biohelpers/issues)中报告bug或建议新功能。

## 许可证

本项目采用MIT许可证 - 查看[LICENSE.md](LICENSE.md)文件了解详情。

## 作者

* **Xiang LI** - *项目维护者* - [lixiang117423](https://github.com/lixiang117423)

## 致谢

感谢所有为此项目做出贡献的开发者和用户。
