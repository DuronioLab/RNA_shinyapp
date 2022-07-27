##                          RNA-seq shiny app                                 ##
## -------------------------------------------------------------------------- ##
##                                                                            ##
##                                                                            ##
## Running the code on this script will start the Shiny App on your computer, ##
## assuming all necessary data are in ./input_data!                           ##
##                                                                            ##
##                                                                            ##
## -------------------------------------------------------------------------- ##



# To run, highlight all text and hit "ctrl+enter" or "cmnd+enter"

library(shiny, warn.conflicts = F, quietly = T)

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

runApp()
