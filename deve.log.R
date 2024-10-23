# add MIT license
usethis::use_mit_license()

# check
devtools::check()

# add function PCAinALL
# usethis::use_r("PCAinALL")
styler::style_file("R/PCAinALL.R")
usethis::use_package("FactoMineR")
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
