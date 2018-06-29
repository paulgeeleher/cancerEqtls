#' Set the root directory where the data will be stored. NB: this directory needs to be set / created based on your own system!! The file "theRootDir.R" is sourced by most of these scripts and should include the desired root directory used by your entire analysis.
theRootDir <- "/mnt/data_scratch/prediXcanProj/"

#' Note, I'm creating a new operator to concatenate strings:
'%&%' <- function(x, y)paste(x,y, sep= "")
