library(testthat)
library(GenomicLLM)
library(GenomicRanges)
library(data.table)

test_that("GRangesData class can be instantiated", {
  gr <- GRanges(
    seqnames = Rle(factor(c("chr1", "chr2"))),
    ranges = IRanges(start = c(10000, 20000), end = c(11000, 21000)),
    strand = Rle(c("+", "-"))
  )
  metadata <- data.table(name = c("gene1", "gene2"), type = c("type1", "type2"), feature_level = c(10, 20))
  gr_data <- new("GRangesData", granges = gr, metadata = metadata)
  expect_s4_class(gr_data, "GRangesData")
  expect_equal(as.character(seqnames(gr_data@granges)), c("chr1", "chr2"))
})

test_that("ComparisonData class can be instantiated", {
  gr1 <- GRanges(
    seqnames = Rle(factor(c("chr1", "chr2"))),
    ranges = IRanges(start = c(10000, 20000), end = c(11000, 21000)),
    strand = Rle(c("+", "-"))
  )
  metadata1 <- data.table(name = c("gene1", "gene2"), type = c("type1", "type2"), feature_level = c(10, 20))
  gr_data1 <- new("GRangesData", granges = gr1, metadata = metadata1)

  gr2 <- GRanges(
    seqnames = Rle(factor(c("chr1", "chr2"))),
    ranges = IRanges(start = c(10500, 20500), end = c(11500, 21500)),  # Ensure overlap
    strand = Rle(c("+", "-"))
  )
  metadata2 <- data.table(name = c("gene3", "gene4"), type = c("type3", "type4"), feature_level = c(30, 40))
  gr_data2 <- new("GRangesData", granges = gr2, metadata = metadata2)

  comparison <- new("ComparisonData", condition1 = gr_data1, condition2 = gr_data2, features = c("methylation"), ranges = "upstream", bp_window = 1000, output_consequence = "isoform_switch_consequence", system_prompt = "Analyze the following sequences:")
  expect_s4_class(comparison, "ComparisonData")
  expect_equal(comparison@features, c("methylation"))
})

test_that("readFile function reads BED files correctly", {
  temp_bed <- tempfile(fileext = ".bed")
  writeLines("chr1\t9999\t10999\tgene1\t0\t+\t9999\t10999\t0\t1\t1000\t0\ttype1\t10\nchr2\t19999\t20999\tgene2\t0\t-\t19999\t20999\t0\t1\t1000\t0\ttype2\t20", temp_bed)
  gr_data <- readFile(temp_bed, format = "BED")
  expect_s4_class(gr_data, "GRangesData")
  expect_equal(as.character(seqnames(gr_data@granges)), c("chr1", "chr2"))
})

test_that("aggregateFeature function aggregates feature data correctly", {
  gr1 <- GRanges(
    seqnames = Rle(factor(c("chr1", "chr2"))),
    ranges = IRanges(start = c(10000, 20000), end = c(11000, 21000)),
    strand = Rle(c("+", "-"))
  )
  metadata1 <- data.table(name = c("gene1", "gene2"), type = c("type1", "type2"), feature_level = c(10, 20))
  gr_data1 <- new("GRangesData", granges = gr1, metadata = metadata1)

  gr2 <- GRanges(
    seqnames = Rle(factor(c("chr1", "chr2"))),
    ranges = IRanges(start = c(10500, 20500), end = c(11500, 21500)),  # Ensure overlap
    strand = Rle(c("+", "-"))
  )
  metadata2 <- data.table(name = c("gene3", "gene4"), type = c("type3", "type4"), feature_level = c(30, 40))
  gr_data2 <- new("GRangesData", granges = gr2, metadata = metadata2)

  comparison <- new("ComparisonData", condition1 = gr_data1, condition2 = gr_data2, features = c("methylation"), ranges = "upstream", bp_window = 1000, output_consequence = "isoform_switch_consequence", system_prompt = "Analyze the following sequences:")
  comparison <- aggregateFeature(comparison)
  expect_s4_class(comparison, "ComparisonData")
  expect_equal(names(comparison@results), c("feature", "avg_feature1", "avg_feature2"))
})

test_that("normalizeFeatures function normalizes feature data correctly", {
  gr1 <- GRanges(
    seqnames = Rle(factor(c("chr1", "chr2"))),
    ranges = IRanges(start = c(10000, 20000), end = c(11000, 21000)),
    strand = Rle(c("+", "-"))
  )
  metadata1 <- data.table(name = c("gene1", "gene2"), type = c("type1", "type2"), feature_level = c(10, 20))
  gr_data1 <- new("GRangesData", granges = gr1, metadata = metadata1)

  gr2 <- GRanges(
    seqnames = Rle(factor(c("chr1", "chr2"))),
    ranges = IRanges(start = c(10500, 20500), end = c(11500, 21500)),  # Ensure overlap
    strand = Rle(c("+", "-"))
  )
  metadata2 <- data.table(name = c("gene3", "gene4"), type = c("type3", "type4"), feature_level = c(30, 40))
  gr_data2 <- new("GRangesData", granges = gr2, metadata = metadata2)

  comparison <- new("ComparisonData", condition1 = gr_data1, condition2 = gr_data2, features = c("methylation"), ranges = "upstream", bp_window = 1000, output_consequence = "isoform_switch_consequence", system_prompt = "Analyze the following sequences:")
  comparison <- normalizeFeatures(comparison, method = "z-score")

  expect_equal(mean(comparison@condition1@metadata$feature_level, na.rm = TRUE), 0)
  expect_equal(sd(comparison@condition1@metadata$feature_level, na.rm = TRUE), 1)
  expect_equal(mean(comparison@condition2@metadata$feature_level, na.rm = TRUE), 0)
  expect_equal(sd(comparison@condition2@metadata$feature_level, na.rm = TRUE), 1)
})

test_that("toLLMFormat function converts data to JSON format correctly", {
  gr1 <- GRanges(
    seqnames = Rle(factor(c("chr1", "chr2"))),
    ranges = IRanges(start = c(10000, 20000), end = c(11000, 21000)),
    strand = Rle(c("+", "-"))
  )
  metadata1 <- data.table(name = c("gene1", "gene2"), type = c("type1", "type2"), feature_level = c(10, 20))
  gr_data1 <- new("GRangesData", granges = gr1, metadata = metadata1)

  gr2 <- GRanges(
    seqnames = Rle(factor(c("chr1", "chr2"))),
    ranges = IRanges(start = c(10500, 20500), end = c(11500, 21500)),  # Ensure overlap
    strand = Rle(c("+", "-"))
  )
  metadata2 <- data.table(name = c("gene3", "gene4"), type = c("type3", "type4"), feature_level = c(30, 40))
  gr_data2 <- new("GRangesData", granges = gr2, metadata = metadata2)

  comparison <- new("ComparisonData", condition1 = gr_data1, condition2 = gr_data2, features = c("methylation"), ranges = "upstream", bp_window = 1000, output_consequence = "isoform_switch_consequence", system_prompt = "Analyze the following sequences:")
  comparison <- aggregateFeature(comparison)
  llm_data <- toLLMFormat(comparison, split_ratio = 0.8, style = "Bracket-Spacing")
  expect_s4_class(llm_data, "LLMFormatData")
  expect_true(length(llm_data@training_data) > 0)
  expect_true(length(llm_data@test_data) > 0)
})

test_that("createFormattedInput handles masking correctly", {
  results <- data.frame(
    avg_feature1 = c(10, 20, 30, 40),
    avg_feature2 = c(15, 25, 35, 45),
    gap = c(5, 10, 15, 20)
  )
  mask_ranges <- list(c(1, 2))

  input_bracket_spacing <- createFormattedInput(results, style = "Bracket-Spacing", mask_ranges = mask_ranges)
  expect_equal(input_bracket_spacing, "[NA, NA] {5} [NA, NA] {10} [30, 35] {15} [40, 45]")

  input_letter_digit <- createFormattedInput(results, style = "Letter-Digit", mask_ranges = mask_ranges)
  expect_equal(input_letter_digit, "A(NA,NA) S(5) A(NA,NA) S(10) A(30,35) S(15) A(40,45)")

  input_verbose <- createFormattedInput(results, style = "Verbose", mask_ranges = mask_ranges)
  expect_equal(input_verbose, "Feature1: NA, Feature2: NA, Gap: 5Feature1: NA, Feature2: NA, Gap: 10Feature1: 30, Feature2: 35, Gap: 15Feature1: 40, Feature2: 45, Gap: 20")
})
