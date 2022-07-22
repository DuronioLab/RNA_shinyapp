## Running the code on this script will start the Shiny App locally on your computer, assuming all necessary data are here!

# To run, highlight all text and hit "ctrl+enter" or "cmnd+enter"

library(shiny, warn.conflicts = F, quietly = T)

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(wd)
setwd(wd)


runApp()
