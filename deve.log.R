# check r package name availability
available::available("biohelpers")

# add MIT license
usethis::use_mit_license()

# add .Rbuildignore
usethis::use_build_ignore("README.md")
usethis::use_build_ignore("deve.log.R")
# usethis::use_build_ignore("LICENSE")
# usethis::use_build_ignore("LICENSE.md")
usethis::use_build_ignore("biohelpers.Rproj")

# add function PCAinALL
# usethis::use_r("PCAinALL")
usethis::rename_files("pca_in_all.R", "pca_in_one.R")
styler::style_file("R/pca_in_all.R")
# usethis::use_package("FactoMineR")
devtools::document()

# 添加用到的R包和其中的部分功能
# 会自动添加到R/biohelpers-package.R这个文件中
usethis::use_import_from("FactoMineR", "PCA")
usethis::use_import_from("factoextra", "get_eigenvalue")
usethis::use_import_from("dplyr", "mutate")
usethis::use_import_from("dplyr", "arrange")
usethis::use_import_from("base", "paste0")
usethis::use_import_from("base", "factor")
usethis::use_import_from("base", "unique")
usethis::use_import_from("base", "round")
usethis::use_import_from("dplyr", "select")
usethis::use_import_from("base", "as.data.frame")
# usethis::use_import_from("dpl")

# 格式化DESCIPTION
usethis::use_tidy_description()

# 添加引用
usethis::use_citation()

# 添加github链接和bug提交链接
usethis::use_github_links()

# 添加示例数据,直接存储为.rda格式
test_data = iris
usethis::use_data(test_data)

# 修改版本
usethis::use_version("major") # 第一位数字
usethis::use_version("minor") # 第二位数字
usethis::use_version("patch") # 第三位数字
usethis::use_version("dev") # 第四位数字

# check
devtools::load_all()
devtools::check()


