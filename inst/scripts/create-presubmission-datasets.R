library(magrittr)
library(stringi)
library(assertive.base)

generateNumberStrings <- function(nNumbers, nCharsPerNumber)
{
  # See http://stackoverflow.com/q/37298911/134830
  stri_rand_strings(nNumbers, nCharsPerNumber, "[0-9]")
}

write.csv0 <- function(x)
{
  fname <- paste0(get_name_in_parent(x), ".csv")
  utils::write.table(
    x,
    file      = file.path("inst/extdata", fname),
    sep       = ",",
    row.names = FALSE,
    col.names = FALSE,
    qmethod   = "double"
  )
}

nSamplesPerPlate <- 96

# Generate Example Comments -----------------------------------------------

nComments <- 20
comments <- data.frame(
  PlatePositions = sample(readat:::PLATE_POSITIONS$PlatePosition, nComments),
  SampleNotes    = sample(
    c("", "red", "yellow", "turbid", "red, turbid", "yellow, turbid"),
    nComments,
    replace = TRUE
  ),
  AssayNotes     = sample(
    c("", "centrifuged at 14000g", "subarray leakage"),
    nComments,
    replace = TRUE,
    prob = c(10, 2, 1)
  )
)



# Generate Example Controls -----------------------------------------------

nControls <- 12
nCharsInBarCode <- 6
barCodes <- generateNumberStrings(nControls, nCharsInBarCode) %>%
  paste0("I", .)
# Controls appear in order A1, A2, A3, ... rather than A1, B1, C1, ... as
# seen in the submission form.
pp <- readat:::PLATE_POSITIONS[
  with(readat:::PLATE_POSITIONS, order(Subarray, Slide)),
  "PlatePosition"
]
controls <- data.frame(
  PlatePositions = pp,
  BarCode        = `[<-`(
    rep.int(NA_character_, nSamplesPerPlate),
    sample(nSamplesPerPlate, nControls),
    barCodes
  )
)


# Generate Example Samples ------------------------------------------------

nCharsInSampleId <- 9
nActualSamples <- nSamplesPerPlate - nControls
sampleIds <- generateNumberStrings(nActualSamples, nCharsInSampleId)
samples <- data.frame(
  PlatePositions = pp,
  SampleId        = `[<-`(
    rep.int("NO READ", nSamplesPerPlate), is.na(controls$BarCode), sampleIds)
)

# Generate Example Slides -------------------------------------------------

nSlides <- 12
nCharsIdSlideId <- 12
slides <- data.frame(
  SlideId = generateNumberStrings(nSlides, nCharsIdSlideId)
)

write.csv0(comments)
write.csv0(controls)
write.csv0(samples)
write.csv0(slides)


