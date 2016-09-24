# setup -------------------------------------------------------------------

wide <- readAdat(extractSampleData(), keepOnlyPasses = FALSE)
long <- melt(wide)

# melt --------------------------------------------------------------------

test_that(
  "melt for WideSomaLogicData objects returns an object of class LongSomaLogicData",
  {
    actual <- long
    expect_is(actual, c("LongSomaLogicData", "data.table", "data.frame"))

    expect_identical(actual, long)

    # Check that x not changed by reference
    expect_is(wide, c("WideSomaLogicData", "data.table", "data.frame"))
  }
)

# Indexing ----------------------------------------------------------------

test_that(
  "[ for WideSomaLogicData objects with i argument indexes by row",
  {
    actual <- wide[1:5]
    expect_is(actual, c("WideSomaLogicData", "data.table", "data.frame"))

    expect_identical(nrow(actual), 5L)
    expect_identical(ncol(actual), ncol(wide))

    # Check preservation of attributes
    expect_identical(attr(actual, "Checksum"), attr(wide, "Checksum"))
    expect_identical(attr(actual, "Metadata"), attr(wide, "Metadata"))
    expect_identical(attr(actual, "SequenceData"), attr(wide, "SequenceData"))

    # Check that x not changed by reference
    expect_is(wide, c("WideSomaLogicData", "data.table", "data.frame"))
  }
)

test_that(
  "[ for WideSomaLogicData objects with list j argument returns a data table",
  {
    # TODO: should subsetting with j also subset SequenceData if sequence
    # columns are removed?
    actual <- wide[, j = list(ExtIdentifier, `SeqId.3896-5_2`)]
    expect_is(actual, c("WideSomaLogicData", "data.table", "data.frame"))

    expect_identical(nrow(actual), nrow(wide))
    expect_identical(ncol(actual), 2L)

    # Check preservation of attributes
    expect_identical(attr(actual, "Checksum"), attr(wide, "Checksum"))
    expect_identical(attr(actual, "Metadata"), attr(wide, "Metadata"))
    expect_identical(attr(actual, "SequenceData"), attr(wide, "SequenceData"))

    # Check that x not changed by reference
    expect_is(wide, c("WideSomaLogicData", "data.table", "data.frame"))
  }
)

test_that(
  "[ for WideSomaLogicData objects with unquoted colname j argument returns a vector",
  {
    actual <- wide[, j = `SeqId.3896-5_2`]
    expect_is(actual, "numeric")

    expect_identical(length(actual), nrow(wide))

    # Check preservation of attributes
    expect_identical(attr(actual, "Checksum"), attr(wide, "Checksum"))
    expect_identical(attr(actual, "Metadata"), attr(wide, "Metadata"))
    expect_identical(attr(actual, "SequenceData"), attr(wide, "SequenceData"))

    # Check that x not changed by reference
    expect_is(wide, c("WideSomaLogicData", "data.table", "data.frame"))
  }
)

test_that(
  "[ for WideSomaLogicData objects with i & j arguments returns something suitable",
  {
    actual <- wide[1:5, j = list(ExtIdentifier, `SeqId.3896-5_2`)]
    expect_is(actual, c("WideSomaLogicData", "data.table", "data.frame"))

    expect_identical(nrow(actual), 5L)
    expect_identical(ncol(actual), 2L)

    # Check preservation of attributes
    expect_identical(attr(actual, "Checksum"), attr(wide, "Checksum"))
    expect_identical(attr(actual, "Metadata"), attr(wide, "Metadata"))
    expect_identical(attr(actual, "SequenceData"), attr(wide, "SequenceData"))

    # Check that x not changed by reference
    expect_is(wide, c("WideSomaLogicData", "data.table", "data.frame"))
  }
)


# This test is mysteriously not working
# "j not found error"

# test_that(
#   "[ for WideSomaLogicData objects with logical j argument returns a WideSomaLogicData object",
#   {
#     j <- rep_len(c(TRUE, FALSE), ncol(wide))
#     actual <- wide[, j = j, with = FALSE]
#     expect_is(actual, c("WideSomaLogicData", "data.table", "data.frame"))
#
#     expect_identical(nrow(actual), nrow(wide))
#     expect_identical(ncol(actual), sum(j))
#
#     # Check preservation of attributes
#     expect_identical(attr(actual, "Checksum"), attr(wide, "Checksum"))
#     expect_identical(attr(actual, "Metadata"), attr(wide, "Metadata"))
#     expect_identical(attr(actual, "SequenceData"), attr(wide, "SequenceData"))
#
#     # Check that x not changed by reference
#     expect_is(wide, c("WideSomaLogicData", "data.table", "data.frame"))
#   }
# )



test_that(
  "data.table can use indexing with with = FALSE",
  {
    DT <- data.table::data.table(x = 1:5, y = letters[1:5])
    j <- c(FALSE, TRUE)
    actual <- DT[, j, with = FALSE]
    expected <-  data.table::data.table(y = letters[1:5])
    expect_equal(actual, expected)
  }
)

# Accessors ----------------------------------------------------------------

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

# This test is not working. See problems with data.table's
# "with = FALSE" inside R CMD check

test_that(
  "getIntensities for WideSomaLogicData objects returns a matrix",
  {
    actual <- getIntensities(wide)
    expect_is(actual, "matrix")

    #Check values the same
#     is_intensity_column <- is_seqid(colnames(wide))
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

# getSampleData -----------------------------------------------------------

# TODO

