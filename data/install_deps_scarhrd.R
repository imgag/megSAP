install.packages("optparse", repos="https://cran.rstudio.com")
install.packages("devtools", repos="https://cran.rstudio.com")
install.packages("sequenza", repos="https://cran.rstudio.com")

library(devtools)
install_github('sztup/scarHRD',build_vignettes = FALSE)
install_github('aroneklund/copynumber')
