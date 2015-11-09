context("Test setup")

zipFiles <- dir(
  system.file("extdata", package = "readat"),
  pattern = "\\.adat\\.zip$",
  full.names = TRUE
)
tmp <- tempfile("readat-tests")
adatFiles <- vapply(
  zipFiles,
  function(x) unzip(x, exdir = tmp),
  character(1)
)

context("Test readAdat")

test_that(
  "readAdat can read all files in the extdata dir",
  {
    expectedNSamples <- rep_len(20L, length(adatFiles))
    expectedNSeqs <- ifelse(grepl("1.3k", adatFiles, fixed = TRUE), 1310L, 1129L)
    actual <- lapply(adatFiles, readAdat, keepOnlyPasses = FALSE)
    for(i in seq_along(actual))
    {
      expect_is(actual[[i]], "WideSomaLogicData")
      expect_identical(nrow(actual[[i]]), expectedNSamples[i])
      expect_identical(nrow(getSequenceData(actual[[i]])), expectedNSeqs[i])
    }
  }
)

context("Test teardown")

unlink(adatFiles, force = TRUE)





