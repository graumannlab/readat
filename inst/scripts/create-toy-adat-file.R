library(readat)
library(data.table)
library(digest)
library(dplyr)
library(magrittr)
library(pdapmain)
library(stringr)

# Metadata
metadata <- list(
  Version = "1.2",
  AssayType = "PharmaServices",
  AssayRobot = "Beckman BioMek Fx",
  AssayVersion = "v3.1",
  CreatedBy = "GLP-Like Assay Experiment",
  CreatedDate = format(Sys.Date()),
  DerivedFrom = "/some/file/path/some_file.adat",
  EnteredBy = "Richie",
  ExpDate = toString(Sys.Date() - c(7, 5, 2)),
  ExpIds = "",
  GeneratedBy = "",
  MasterMixVersion = "V3.2 Plasma",
  PlateMedianCal_Plate_A = 1, # Fix these later
  PlateMedianCal_Plate_B = 1,
  PlateMedianCal_Plate_C = 1,
  PlateMedianCal_Plate_D = 1,
  PlateMedianCal_Plate_E = 1,
  ProcessSteps = "/path/to/methodology/file.txt",
  ProteinEffectiveDate = format(Sys.Date()),
  ReportType = "Prepared for Release",
  StudyMatrix = "EDTA Plasma",
  StudyOrganism = "Human",
  Title = ""
)

#################################################

# Dims
nAptamers <- 1129
nSamples <- 1000

#################################################

# Sequence data
sequenceData <- ids[
  j = list(
    SeqId, SomaId, Target, TargetFullName,
    UniProt = vapply(UniProtId, paste0, character(1), collapse = " "),
    EntrezGeneID = vapply(EntrezGeneId, paste0, character(1), collapse = " "),
    EntrezGeneSymbol = vapply(EntrezGeneSymbol, paste0, character(1), collapse = " "),
    Type)
]
sequenceData <- sequenceData[
  j = `:=`(
    Organism = "Human",
    ColCheck = ifelse(runif(nAptamers) < 0.95, "PASS", "FAIL"),
    Units = "RFU"
  )
]

#################################################

# Sample data
sampleData <- data.table(
  PlateId   = gl(5, nSamples / 5, labels = paste("Plate", LETTERS[1:5])),
  SlideId   = sample(gl(10, nSamples / 10, labels = paste("Slide", 1:10))),
  SampleId  = factor(paste("Sample", seq_len(nSamples))),
  ExtIdentifier = factor(paste0("EID", 123455 + seq_len(nSamples))),
  SubjectId = gl(nSamples / 5, 1, nSamples, labels = paste("Subject", seq_len(nSamples / 5))),
  SiteId    = factor("Site 1"),
  NormScale = rnorm(nSamples, 1, 0.05)
)
sampleData$RowCheck  = with(
  sampleData,
  ifelse(NormScale >= 0.9 && NormScale <= 1.1, "PASS", "FAIL")
)

#################################################

# Add in calibrations
medianCalibrationByPlate <- sampleData[
  j = list(MedianCalibration = median(NormScale)),
  by = list(PlateId)
]
medianCalibrationByPlate <- medianCalibrationByPlate[
  j = PlateName := paste0("PlateMedianCal_", str_replace_all(PlateId, " ", "_"))
]
medianCalibrationByPlateList <- with(
  medianCalibrationByPlate,
  setNames(as.list(MedianCalibration), PlateName)
)

metadata_names <- names(metadata)
metadata <- merge(medianCalibrationByPlateList, metadata, warn_on_dupes = FALSE)[metadata_names]

#################################################

# Intensities
intensityMeans <- rlnorm(nAptamers, 7, 1.5)
sampleSds <- rf(nAptamers, nAptamers, nAptamers)
intensities <- Map(
  function(m, s) rnorm(nSamples, m, sqrt(m) * s),
  intensityMeans,
  sampleSds
) %>%
  #setNames(paste0("SeqId.", sequenceData$SeqId)) %>%
  setNames(rep_len(" ", nAptamers)) %>%
  as.data.table

#################################################

# Add in calibration references
# This is not how SomaLogic calculated these values, but it's mostly quite close
# (I think)
sequenceCalibrationReferences <- intensities[j = lapply(.SD, mean)] %>%
  unlist(use.names = FALSE)

tmp <- intensities[
  j = Map(function(x, y) mean(x) / y, .SD, sequenceCalibrationReferences),
  by = sampleData$PlateId
] %>%
  t
sequenceCalibrations <-tmp[-1,] %>%
  as.data.table %>%
  setNames(paste0("Cal_", str_replace_all(levels(sampleData$PlateId), " ", "_")))

sequenceData <- bind_cols(
  sequenceData,
  sequenceCalibrations,
  data.table(
    CalReference = sequenceCalibrationReferences,
    Dilution     = sample(c(0.005, 1, 40), nAptamers, replace = TRUE)
  )
)

#################################################

# Combine things
adatData <- bind_cols(
  sampleData,
  data.frame(" " = character(nSamples), check.names = FALSE),
  intensities
)

checksum <- digest(adatData, "sha1")

#################################################

# Export
adatFile <- "readat/inst/extdata/sample.adat"

# setattr(adatData, "SequenceData", sequenceData)
# setattr(adatData, "Metadata", metadata)
# setattr(adatData, "class", c("WideSomaLogicData", "data.table", "data.frame"))

getDataType <- function(x)
{
  # lists will be unlisted when written to file.
  class_x <- class(unlist(x))
  # Wildly guessing that ADAT uses Visual Basic data types.  The only value
  # in our test files is 'String'; the rest are guesses.
  # https://msdn.microsoft.com/en-us/library/47zceaw7.aspx
  switch(
    class_x,
    character = ,
    factor    = "String",
    numeric   = "Double",
    integer   = "Integer",
    logical   = "Boolean",
    stop("Data type ", class_x, "not yet implemented.")
  )
}

# Now write to file!
writeToAdat <- function()
{
  fcon <- file(adatFile, "w+")
  on.exit(close(fcon))
  cat(paste0("!Checksum\t", checksum, "\n"), file = fcon)
  cat("^HEADER\n", file = fcon)
  cat(
    paste0(
      paste0("!", names(metadata), "\t", metadata),
      collapse = "\n"
    ),
    file = fcon
  )
  cat("\n^COL_DATA\n", file = fcon)
  cat(
    paste0(c("!Name", colnames(sequenceData)), collapse = "\t"),
    fill = TRUE,
    file = fcon
  )
  cat(
    paste0(
      c("!Type", vapply(sequenceData, getDataType, character(1))),
      collapse = "\t"
    ),
    fill = TRUE,
    file = fcon
  )
  cat("^ROW_DATA\n", file = fcon)
  cat(
    paste0(c("!Name", colnames(sampleData)), collapse = "\t"),
    fill = TRUE,
    file = fcon
  )
  cat(
    paste0(
      c("!Type", vapply(sampleData, getDataType, character(1))),
      collapse = "\t"
    ),
    fill = TRUE,
    file = fcon
  )
  cat("^TABLE_BEGIN\n", file = fcon)
  blank_cols <- replicate(
    ncol(sampleData),
    character(ncol(sequenceData)),
    simplify = FALSE
  ) %>%
    as.data.frame()
  pdapmain::write_table(
    bind_cols(
      blank_cols,
      data.frame(" " = colnames(sequenceData), check.names = FALSE),
      as.data.frame(t(sequenceData))
    ),
    fcon,
    sep       = "\t",
    append    = TRUE,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  pdapmain::write_table(
    adatData,
    fcon,
    sep       = "\t",
    append    = TRUE,
    quote     = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}
writeToAdat()

