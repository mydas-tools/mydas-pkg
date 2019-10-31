# #https://cran.r-project.org/web/packages/qicharts2/vignettes/qicharts2.html
# 
# x <- data.frame(
# `Number of useful observations`       = n,
# `Upper limit for longest run`         = round(log2(n) + 3),
# `Lower limit for number of crossings` = qbinom(0.05, n - 1, 0.5),
# check.names                           = FALSE)
# #knitr::kable(x)
