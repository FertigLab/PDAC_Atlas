# Author: Jacob Mitchell

# Creation of initial file structure

# dir for R scripts
if (!dir.exists("R")) {
  dir.create("R")
}

# dir for the combined cds object
if (!dir.exists("data")) {
  dir.create("data")
}

# dir for data passed to CoGAPS and the params given to the algorithm
if (!dir.exists("CoGAPS")) {
  dir.create("CoGAPS")
}

# dir for output of the CoGAPS algorithm
if (!dir.exists("CoGAPS_output")) {
  dir.create("CoGAPS_output")
}

# dir for resulting cds object containing integrated CoGAPS results
if (!dir.exists("results")) {
  dir.create("results")
}



