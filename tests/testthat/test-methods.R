# setup -------------------------------------------------------------------

library(reshape2)
library(stringi)

wide <- readAdat(extractSampleData(), keepOnlyPasses = FALSE)
long <- melt(wide)

# accessors ----------------------------------------------------------------

# getChecksum -------------------------------------------------------------

test_that(
  "getChecksum for WideSomaLogicData objects returns a string",
  {
    actual <- getChecksum(wide)
    expect_is(actual, "character")
    expect_identical(length(actual), 1L)

    expected <- attr(wide, "Checksum")
    expect_identical(actual, expected)

    # Check that x not changed by reference
    expect_is(wide, c("WideSomaLogicData", "data.table", "data.frame"))
  }
)

test_that(
  "getChecksum for LongSomaLogicData objects returns a string",
  {
    actual <- getChecksum(long)
    expect_is(actual, "character")
    expect_identical(length(actual), 1L)

    expected <- attr(long, "Checksum")
    expect_identical(actual, expected)

    # Check that x not changed by reference
    expect_is(long, c("LongSomaLogicData", "data.table", "data.frame"))
  }
)

# getMetadata -------------------------------------------------------------

test_that(
  "getChecksum for WideSomaLogicData objects returns a list",
  {
    actual <- getMetadata(wide)
    expect_is(actual, "list")

    expected <- attr(wide, "Metadata")
    expect_identical(actual, expected)

    # Check that x not changed by reference
    expect_is(wide, c("WideSomaLogicData", "data.table", "data.frame"))
  }
)

test_that(
  "getMetadata for LongSomaLogicData objects returns a list",
  {
    actual <- getMetadata(long)
    expect_is(actual, "list")

    expected <- attr(long, "Metadata")
    expect_identical(actual, expected)

    # Check that x not changed by reference
    expect_is(long, c("LongSomaLogicData", "data.table", "data.frame"))
  }
)

# getSequences ------------------------------------------------------------

test_that(
  "getSequences for WideSomaLogicData objects returns a data table",
  {
    actual <- getSequenceData(wide)
    expect_is(actual, c("data.table", "data.frame"))

    expected <- attr(wide, "SequenceData")
    expect_identical(actual, expected)

    # Check that x not changed by reference
    expect_is(wide, c("WideSomaLogicData", "data.table", "data.frame"))
  }
)

test_that(
  "getSequences for LongSomaLogicData objects returns a data table",
  {
    actual <- getSequenceData(long)
    expect_is(actual, c("data.table", "data.frame"))

    expected <- attr(long, "SequenceData")
    expect_identical(actual, expected)

    # Check that x not changed by reference
    expect_is(long, c("LongSomaLogicData", "data.table", "data.frame"))
  }
)

# getIntensities ----------------------------------------------------------

test_that(
  "getIntensities for WideSomaLogicData objects returns a matrix",
  {
    actual <- getIntensities(wide)
    expect_is(actual, "matrix")

    #Check values the same
#     is_intensity_column <- stri_detect_regex(colnames(wide), "^SeqId\\.")
#     expected_values <- unlist(
#       wide[, is_intensity_column, with = FALSE],
#       use.names =  FALSE
#     )
#     actual_values <- as.numeric(actual)
#     expect_equal(actual_values, expected_values)

    # Check that x not changed by reference
    expect_is(wide, c("WideSomaLogicData", "data.table", "data.frame"))
  }
)

test_that(
  "getIntensities for LongSomaLogicData objects returns a data table",
  {
    actual <- getIntensities(long)
    expect_is(actual, c("data.table", "data.frame"))
    # Check that x not changed by reference
    expect_is(long, c("LongSomaLogicData", "data.table", "data.frame"))
  }
)


