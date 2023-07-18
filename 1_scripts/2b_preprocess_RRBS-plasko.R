# Dependencies
library(tidyverse)
library(readxl)
library(methylKit)

clinic_plasko <- read_excel("0_data/clinic_revised_full_rrbs_wts.xlsx",
                            sheet = "plasko_komplet"
)

# Custom function for fetching paths
fetchFiles <- function(path, pattern) {
  filePaths <- list.files(path, full.names = TRUE, pattern = pattern)
  fileNames <- list.files(path, full.names = FALSE, pattern = pattern)
  return(as.list(filePaths[gsub("^([^^_]*).*", "\\1", fileNames) %in% clinic_plasko$pacjent]))
}

# Custom function for getting sample id
getSampleId <- function(filePaths) {
  sapply(filePaths, function(x) {
    splitPath <- strsplit(x, "/")[[1]]
    splitFileName <- strsplit(splitPath[length(splitPath)], "_")[[1]]
    paste0(splitFileName[1], "_", splitFileName[2])
  })
}

# Fetch files with the custom function
Tn_plasko_f <- fetchFiles("0_data/bisCov_profiles/plasko/Tn", "*Tn*")
Tc_plasko_f <- fetchFiles("0_data/bisCov_profiles/plasko/Tc", "*Tc*")
Tp_plasko_f <- fetchFiles("0_data/bisCov_profiles/plasko/Tp", "*Tp*")

# Create Tabix Methyl DB
Tc_grczuo_DB <- methRead(
  c(Tc_plasko_f, Tn_plasko_f),
  sample.id = c(getSampleId(Tc_plasko_f), getSampleId(Tn_plasko_f)),
  assembly = "GRCh38",
  treatment = c(rep(1, length(Tc_plasko_f)), rep(0, length(Tn_plasko_f))),
  context = "CpG",
  pipeline = "bismarkCoverage",
  header = FALSE,
  dbtype = "tabix",
  dbdir = "0_data/methylKit_DBs/methylDB_plasko_Tc_vs_Tn",
  mincov = 10 # default
)

Tp_grczuo_DB <- methRead(
  c(Tp_plasko_f, Tn_plasko_f),
  sample.id = c(getSampleId(Tp_plasko_f), getSampleId(Tn_plasko_f)),
  assembly = "GRCh38",
  treatment = c(rep(1, length(Tc_plasko_f)), rep(0, length(Tn_plasko_f))),
  context = "CpG",
  pipeline = "bismarkCoverage",
  header = FALSE,
  dbtype = "tabix",
  dbdir = "0_data/methylKit_DBs/methylDB_plasko_Tp_vs_Tn",
  mincov = 10 # default
)
