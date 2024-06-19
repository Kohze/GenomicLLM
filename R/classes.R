#' GRangesData Class
#'
#' This class represents a GRanges object along with its metadata.
#'
#' @slot granges A GRanges object containing genomic ranges.
#' @slot metadata A data frame containing metadata for the GRanges object.
#' @export
setClass("GRangesData",
         slots = list(
           granges = "GRanges",
           metadata = "data.frame"
         ))

#' ComparisonData Class
#'
#' This class represents the comparison data between two conditions, including the features, ranges, base pair window size, and results.
#'
#' @slot condition1 A GRangesData object for the first condition.
#' @slot condition2 A GRangesData object for the second condition.
#' @slot features A character vector of features to be compared.
#' @slot ranges A character string indicating the range type (upstream, downstream, etc.).
#' @slot bp_window An integer indicating the base pair window size.
#' @slot results A data frame containing comparison results.
#' @slot output_consequence A character string indicating the consequence of the comparison.
#' @slot system_prompt A character string for the system prompt in the LLM JSON output.
#' @export
setClass("ComparisonData",
         slots = list(
           condition1 = "GRangesData",
           condition2 = "GRangesData",
           features = "character",
           ranges = "character",
           bp_window = "numeric",
           results = "data.frame",
           output_consequence = "character",
           system_prompt = "character"
         ))

#' LLMFormatData Class
#'
#' This class represents the LLM formatted data, including training and test datasets.
#'
#' @slot training_data JSON formatted training data.
#' @slot test_data JSON formatted test data.
#' @export
setClass("LLMFormatData",
         slots = list(
           training_data = "character",
           test_data = "character"
         ))
