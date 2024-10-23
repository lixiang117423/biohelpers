# 

# add MIT license
usethis::use_mit_license()

# check
devtools::check()

# add function PCAinALL
usethis::use_r("PCAinALL")
styler::style_file("R/PCAinALL.R")
