##                          RNA-seq shiny app                                 ##
## -------------------------------------------------------------------------- ##
##                                                                            ##
##                                                                            ##
## Running the code on this script will start the Shiny App on your computer, ##
## assuming all necessary data are in ./input_data!                           ##
##                                                                            ##
##                                                                            ##
##       To run, highlight all text and hit "ctrl+enter" or "cmnd+enter"      ##
##                                                                            ##
## -------------------------------------------------------------------------- ##

library(shiny, warn.conflicts = F, quietly = T)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## Optional config if your input_data file names are different than expected ##

annotation_file <<- "dm6_dHisC_chrHis.gtf"
counts_file <<- "featureCounts.txt"
samples_file <<- "sampleSheet.txt"
synonym_file <<- "fb_synonym_fb_2021_03.tsv"


runApp()
