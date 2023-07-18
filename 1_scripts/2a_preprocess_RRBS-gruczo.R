# Dependencies
library(tidyverse)
library(readxl)
library(methylKit)

clinic_gruczo <- read_excel("0_data/clinic_revised_full_rrbs_wts.xlsx",
                            sheet = "gruczolak_komplet"
)

# Custom function for fetching paths
fetchFiles <- function(path, pattern) {
  filePaths <- list.files(path, full.names = TRUE, pattern = pattern)
  fileNames <- list.files(path, full.names = FALSE, pattern = pattern)
  return(as.list(filePaths[gsub("^([^^_]*).*", "\\1", fileNames) %in% clinic_gruczo$pacjent]))
}

# Custom function for getting sample id
getSampleId <- function(filePaths) {
  sapply(filePaths, function(x) {
    splitPath <- strsplit(x, "/")[[1]]
    splitFileName <- strsplit(splitPath[length(splitPath)], "_")[[1]]
    as.list(paste0(splitFileName[1], "_", splitFileName[2]))
  })
}

# Fetch files with the custom function
Tn_gruczo_f <- fetchFiles("0_data/bisCov_profiles/gruczolak/Tn", "*Tn*")
Tc_gruczo_f <- fetchFiles("0_data/bisCov_profiles/gruczolak/Tc", "*Tc*")
Tp_gruczo_f <- fetchFiles("0_data/bisCov_profiles/gruczolak/Tp", "*Tp*")

# Create Tabix Methyl DB
# Create a named vector for treatment
treatment <- c(rep(1, length(Tc_gruczo_f)), rep(0, length(Tn_gruczo_f)))
names(treatment) <- c(getSampleId(Tc_gruczo_f), getSampleId(Tn_gruczo_f))

Tc_grczuo_DB <- methRead(
  location = c(Tc_gruczo_f, Tn_gruczo_f),
  sample.id = c(getSampleId(Tc_gruczo_f), getSampleId(Tn_gruczo_f)),
  assembly = "GRCh38",
  treatment = treatment, 
  context = "CpG",
  pipeline = "bismarkCoverage",
  header = FALSE,
  dbtype = "tabix",
  dbdir = "0_data/methylKit_DBs/methylDB_gruczo_Tc_vs_Tn",
  mincov = 10 # default
)

Tp_grczuo_DB <- methRead(
  c(Tp_gruczo_f, Tn_gruczo_f),
  sample.id = c(getSampleId(Tp_gruczo_f), getSampleId(Tn_gruczo_f)),
  assembly = "GRCh38",
  treatment = c(rep(1, length(Tc_gruczo_f)), rep(0, length(Tn_gruczo_f))),
  context = "CpG",
  pipeline = "bismarkCoverage",
  header = FALSE,
  dbtype = "tabix",
  dbdir = "0_data/methylKit_DBs/methylDB_gruczo_Tp_vs_Tn",
  mincov = 10 # default
)



