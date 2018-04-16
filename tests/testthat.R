library("testthat")
library("assertthat")

library("GeneAccord")
# "When R CMD check runs tests, it sets R_TESTS. When the tests
# themeselves run R CMD xxxx, as is the case with the tests in
# devtools, having R_TESTS set causes errors because it confuses
# the R subprocesses. Unsetting it here avoids those problems.
#"R_TESTS" = "" "
Sys.setenv("R_TESTS" = "")

test_check("GeneAccord")

