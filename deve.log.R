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
usethis::rename_files("PCAinALL.R", "pca_in_all.R")
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

# check
devtools::load_all()
devtools::check()
