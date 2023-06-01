secret <- commandArgs(trailing = TRUE)
devtools::install_github('shahcompbio/signals', dependencies = TRUE, auth_token = secret)
