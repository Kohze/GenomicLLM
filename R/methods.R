# Define the classes
setClass(
  Class = "GRangesData",
  slots = c(
    granges = "GRanges",
    metadata = "data.table"
  )
)

setClass(
  Class = "ComparisonData",
  slots = c(
    condition1 = "GRangesData",
    condition2 = "GRangesData",
    features = "character",
    ranges = "character",
    bp_window = "numeric",
    output_consequence = "character",
    system_prompt = "character",
    results = "data.frame"
  )
)

setClass(
  Class = "LLMFormatData",
  slots = c(
    training_data = "character",
    test_data = "character"
  )
)

#' Read a BED file and create a GRangesData object
#'
#' @param file Path to the BED file
#' @param format Format of the file (default is "BED")
#' @return A GRangesData object
#' @export
#' @importFrom data.table fread
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom methods new
setGeneric("readFile", function(file, ...) standardGeneric("readFile"))

#' @describeIn readFile Method to read BED files
setMethod("readFile", "character",
          function(file, format = "BED") {
            if (format == "BED") {
              bed_data <- fread(file)
              colnames(bed_data) <- c("seqnames", "start", "end", "name", "score", "strand",
                                      "thickStart", "thickEnd", "itemRgb", "blockCount",
                                      "blockSizes", "blockStarts", "type", "feature_level")
              bed_data$start <- bed_data$start + 1
              gr <- makeGRangesFromDataFrame(bed_data, keep.extra.columns = TRUE)
              metadata <- bed_data[, .(name, type, feature_level)]
              return(new("GRangesData", granges = gr, metadata = metadata))
            } else {
              stop("Unsupported file format")
            }
          })

#' Aggregate Feature Data for Comparison
#'
#' @param comparison A ComparisonData object
#' @param upstream_bp Number of base pairs upstream for aggregation (default is 1000)
#' @param mask_ranges A list of numeric vectors specifying ranges to mask (default is NULL)
#' @return A ComparisonData object with aggregated feature data
#' @export
#' @importFrom GenomicRanges resize
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
setGeneric("aggregateFeature", function(comparison, ...) standardGeneric("aggregateFeature"))

#' @describeIn aggregateFeature Method to aggregate feature data
setMethod("aggregateFeature", "ComparisonData",
          function(comparison, upstream_bp = 1000, mask_ranges = NULL) {

            extended_gr1 <- resize(comparison@condition1@granges, fix = "start", width = width(comparison@condition1@granges) + upstream_bp)
            extended_gr2 <- resize(comparison@condition2@granges, fix = "start", width = width(comparison@condition2@granges) + upstream_bp)

            if (!is.null(mask_ranges)) {
              for (mask_range in mask_ranges) {
                mask_start <- mask_range[1]
                mask_end <- mask_range[2]
                extended_gr1 <- extended_gr1[!(start(extended_gr1) >= mask_start & end(extended_gr1) <= mask_end)]
                extended_gr2 <- extended_gr2[!(start(extended_gr2) >= mask_start & end(extended_gr2) <= mask_end)]
              }
            }

            overlaps1 <- findOverlaps(extended_gr1, comparison@condition2@granges)
            overlaps2 <- findOverlaps(extended_gr2, comparison@condition1@granges)

            feature_values1 <- as.numeric(comparison@condition2@granges$feature_level[subjectHits(overlaps1)])
            feature_values2 <- as.numeric(comparison@condition1@granges$feature_level[subjectHits(overlaps2)])

            if (length(feature_values1) == 0 || length(feature_values2) == 0) {
              avg_feature1 <- numeric(0)
              avg_feature2 <- numeric(0)
            } else {
              avg_feature1 <- sapply(split(feature_values1, as.vector(queryHits(overlaps1))), mean, na.rm = TRUE)
              avg_feature2 <- sapply(split(feature_values2, as.vector(queryHits(overlaps2))), mean, na.rm = TRUE)
            }

            len <- min(length(avg_feature1), length(avg_feature2))
            avg_feature1 <- avg_feature1[seq_len(len)]
            avg_feature2 <- avg_feature2[seq_len(len)]

            if (len == 0) {
              avg_feature1 <- NA
              avg_feature2 <- NA
            }

            results <- data.frame(
              feature = comparison@features,
              avg_feature1 = avg_feature1,
              avg_feature2 = avg_feature2
            )
            comparison@results <- results

            return(comparison)
          })

#' Convert Comparison Data to JSON Format
#'
#' @param comparison A ComparisonData object
#' @param split_ratio Ratio for splitting data into training and test sets (default is NULL)
#' @param style The style of the output format (Bracket-Spacing, Letter-Digit, Verbose)
#' @param system_prompt A character string for the system prompt in the LLM JSON output (default is "Analyze the following sequences:")
#' @return A list containing JSON formatted data for training and test sets or a single JSON object
#' @export
#' @importFrom jsonlite toJSON
#' @importFrom methods new
#' @importFrom stats start
setGeneric("toLLMFormat", function(comparison, split_ratio = NULL, style = "Bracket-Spacing", system_prompt = "Analyze the following sequences:") standardGeneric("toLLMFormat"))

#' @describeIn toLLMFormat Method to convert comparison data to JSON format
setMethod("toLLMFormat", "ComparisonData",
          function(comparison, split_ratio = NULL, style = "Bracket-Spacing", system_prompt = "Analyze the following sequences:") {
            if (length(comparison@features) > 10) {
              stop("The number of features exceeds the maximum allowed (10).")
            }

            formatted_data <- data.frame(
              instruction = system_prompt,
              input = paste0("[", comparison@results$avg_feature1, ", ", comparison@results$avg_feature2, "]"),
              output = comparison@output_consequence
            )

            json_output <- toJSON(formatted_data, pretty = TRUE, auto_unbox = TRUE)

            if (!is.null(split_ratio)) {
              n <- nrow(formatted_data)
              train_indices <- sample(1:n, size = floor(split_ratio * n))
              test_indices <- setdiff(1:n, train_indices)

              training_data <- toJSON(formatted_data[train_indices, ], pretty = TRUE, auto_unbox = TRUE)
              test_data <- toJSON(formatted_data[test_indices, ], pretty = TRUE, auto_unbox = TRUE)

              combined_json <- list(training_data = as.character(training_data), test_data = as.character(test_data))
              return(new("LLMFormatData", training_data = combined_json$training_data, test_data = combined_json$test_data))
            } else {
              return(json_output)
            }
          })

#' Create Formatted Input
#'
#' @param results A data frame containing the comparison results
#' @param style The style of the formatted input (Bracket-Spacing, Letter-Digit, Verbose)
#' @param mask_ranges A list of numeric vectors specifying ranges to mask (default is NULL)
#' @return A formatted input string
#' @export
setGeneric("createFormattedInput", function(results, style, mask_ranges = list()) standardGeneric("createFormattedInput"))

setMethod("createFormattedInput", "data.frame",
          function(results, style, mask_ranges = list()) {
            input <- ""

            # Apply masking to results
            for (mask_range in mask_ranges) {
              results$avg_feature1[mask_range] <- NA
              results$avg_feature2[mask_range] <- NA
            }

            if (style == "Bracket-Spacing") {
              for (i in 1:nrow(results)) {
                input <- paste0(input, "[", results$avg_feature1[i], ", ", results$avg_feature2[i], "]")
                if (i < nrow(results)) {
                  input <- paste0(input, " {", results$gap[i], "} ")
                }
              }
            } else if (style == "Letter-Digit") {
              for (i in 1:nrow(results)) {
                input <- paste0(input, "A(", results$avg_feature1[i], ",", results$avg_feature2[i], ")")
                if (i < nrow(results)) {
                  input <- paste0(input, " S(", results$gap[i], ") ")
                }
              }
            } else if (style == "Verbose") {
              for (i in 1:nrow(results)) {
                input <- paste0(input, "Feature1: ", results$avg_feature1[i], ", Feature2: ", results$avg_feature2[i], ", Gap: ", results$gap[i])
              }
            }
            return(input)
          })

#' Normalize Feature Data
#'
#' @param comparison A ComparisonData object
#' @param method The normalization method, either "z-score" or "min-max" (default is "z-score")
#' @return A ComparisonData object with normalized feature data
#' @export
setGeneric("normalizeFeatures", function(comparison, method = "z-score") standardGeneric("normalizeFeatures"))

#' @describeIn normalizeFeatures Method to normalize feature data
setMethod("normalizeFeatures", "ComparisonData",
          function(comparison, method = "z-score") {
            normalize <- function(x, method) {
              if (method == "z-score") {
                return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
              } else if (method == "min-max") {
                return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
              } else {
                stop("Unsupported normalization method")
              }
            }

            comparison@condition1@metadata$feature_level <- normalize(comparison@condition1@metadata$feature_level, method)
            comparison@condition2@metadata$feature_level <- normalize(comparison@condition2@metadata$feature_level, method)

            return(comparison)
          })
