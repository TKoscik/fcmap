# Example use of fcmap ---------------------------------------------------------
# You can edit this script and call it from the terminal using:
# >> Rscript /full/path/to/script.R

# Setup R Environment ----------------------------------------------------------
# install fcmap functions
if (!("devtools" %in% rownames(installed.packages()))) {
  install.packages("devtools")
}
if (!("fcmap" %in% rownames(installed.packages()))) {
  devtools::install_github("TKoscik/fcmap")
}

# Clear R environment
rm(list=ls())
gc()

# Load fcmp functions
library("fcmap")

# Method 1 ---------------------------------------------------------------------
epi.nii <- "Fill in appropriate input"
network_cortex <- "Fill in appropriate input"
network_ROI <- "Fill in appropriate input"
save_file_name <- "Fill in appropriate input"

fcmap.m1(epi.nii, network_cortex, network_ROI, save_file_name)

# Method 2 ---------------------------------------------------------------------
epi.nii <- "Fill in appropriate input"
cortex <- "Fill in appropriate input"
mask_ROI <- "Fill in appropriate input"
save_file_name <- "Fill in appropriate input"
cluster <- NULL # default input
cluster_limit <- NULL # default input

fcmap.m2(epi.nii, cortex, mask_ROI, save_file_name, cluster, cluster_limit)
