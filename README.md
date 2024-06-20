# GenomicLLM: A Comprehensive Package for Preparing Genomic Data for LLM Training

[![codecov](https://codecov.io/gh/username/GenomicLLM/branch/main/graph/badge.svg)](https://codecov.io/gh/username/GenomicLLM)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction

Welcome to GenomicLLM, a powerful R package designed to transform genomic comparison data into formats suitable for training large language models (LLMs). This package facilitates the handling of genomic feature data stored in `GRanges` objects, normalizes features, performs sophisticated data masking, and converts the data into various structured formats for LLM training.

## Features

- **Read and process BED files**: Convert BED files into `GRangesData` objects.
- **Aggregate feature data**: Combine genomic features from multiple conditions.
- **Normalize features**: Apply various normalization techniques to feature data.
- **Mask specific data ranges**: Perform sophisticated data masking.
- **Convert data for LLM training**: Format data into JSON and other structures for training LLMs.
- **Support multiple output styles**: Bracket-Spacing, Letter-Digit, Verbose.

## Installation

To install GenomicLLM, you need to have R installed on your system. You can install the development version from GitHub:

```R
# Install devtools if you don't have it
install.packages("devtools")

# Install GenomicLLM from GitHub
devtools::install_github("username/GenomicLLM")
```

## Usage

### Reading and Processing BED Files

```R
library(GenomicLLM)

# Read a BED file
bed_file <- "path/to/your/file.bed"
gr_data <- readFile(bed_file, format = "BED")
print(gr_data)
```

### Aggregating Feature Data

```R
# Create example GRangesData objects
gr1 <- GRanges(seqnames = Rle(factor(c("chr1", "chr2"))),
               ranges = IRanges(start = c(10000, 20000), end = c(11000, 21000)),
               strand = Rle(c("+", "-")))
metadata1 <- data.table(name = c("gene1", "gene2"), type = c("type1", "type2"), feature_level = c(10, 20))
gr_data1 <- new("GRangesData", granges = gr1, metadata = metadata1)

gr2 <- GRanges(seqnames = Rle(factor(c("chr1", "chr2"))),
               ranges = IRanges(start = c(10500, 20500), end = c(11500, 21500)),
               strand = Rle(c("+", "-")))
metadata2 <- data.table(name = c("gene3", "gene4"), type = c("type3", "type4"), feature_level = c(30, 40))
gr_data2 <- new("GRangesData", granges = gr2, metadata = metadata2)

# Create ComparisonData object
comparison <- new("ComparisonData", condition1 = gr_data1, condition2 = gr_data2,
                  features = c("methylation"), ranges = "upstream", bp_window = 1000,
                  output_consequence = "isoform_switch_consequence", system_prompt = "Analyze the following sequences:")

# Aggregate feature data
comparison <- aggregateFeature(comparison)
print(comparison@results)
```

### Normalizing Feature Data

```R
# Normalize feature data using z-score normalization
comparison <- normalizeFeatures(comparison, method = "z-score")
print(comparison@condition1@metadata)
print(comparison@condition2@metadata)
```

### Converting Data for LLM Training

```R
# Convert data to LLM format with a split ratio
llm_data <- toLLMFormat(comparison, split_ratio = 0.8, style = "Bracket-Spacing")
print(llm_data@training_data)
print(llm_data@test_data)
```

### Data Masking

```R
# Mask specific ranges in the data
mask_ranges <- list(c(1, 2))  # Mask the first two feature pairs
input_bracket_spacing <- createFormattedInput(comparison@results, style = "Bracket-Spacing", mask_ranges = mask_ranges)
print(input_bracket_spacing)
```

### Output Styles

GenomicLLM supports multiple output styles for formatting data. Below are examples for each style.

| Style             | Example Output                                                                                     |
|-------------------|----------------------------------------------------------------------------------------------------|
| Bracket-Spacing   | `[10, 15] {5} [20, 25] {10} [30, 35] {15} [40, 45]`                                                |
| Letter-Digit      | `A(10,15) S(5) A(20,25) S(10) A(30,35) S(15) A(40,45)`                                             |
| Verbose           | `Feature1: 10, Feature2: 15, Gap: 5 Feature1: 20, Feature2: 25, Gap: 10 Feature1: 30, Feature2: 35, Gap: 15 Feature1: 40, Feature2: 45, Gap: 20` |

```R
{ 
  [
    system: 'identify the isoform consequence of this sequence. Options: longer, shorter, same',
    input: '[10, 15] {5} [20, 25] {10} [30, 35] {15} [40, 45]',
    output: 'longer'
  ],
  ...
}
```

### LLM Training Follow UP

- Pick a LLM training tool of choice: We recommend https://github.com/hiyouga/LLaMA-Factory
- Upload the output .json to your huggingface repository
- Train your model (e.g. Llama 3)
- Utilize the masking functionality to discover feature attention and causal relations

## Contributing

We welcome contributions from the community. If you find a bug or have a feature request, please open an issue on GitHub. If you would like to contribute code, please fork the repository and submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgements

We would like to extend our gratitude to the Bioconductor community for their invaluable resources and support in developing this package. Additionally, we would like to thank the Cambridge Department of Genetics and the AF Smith group for their essential contributions and support, which have been instrumental in enabling this research.

## Contact

For further information or questions, please contact robin (a t) gounder (dot) com .

## Citation

We are not yet published. Once available, a citation quote will be added.

---

Happy coding and happy genomics!
