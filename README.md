# NGS_Data_Analyser
this project is about genomic data analysing.

# Genetic Data Analysis Tool

This Python script is a Genetic Data Analysis Tool that allows you to analyze genetic data from different family configurations and generate informative Excel reports based on the analysis results.

## Introduction

The Genetic Data Analysis Tool reads genetic data in CSV format for various family configurations, including father-mother, mother-child, and father-mother-child relationships. It performs filtering and analysis on the data and produces Excel reports containing valuable insights.

## Features

- Supports analysis for different family configurations: `father_mother`, `mother_child`, and `father_mother_child`.
- Filters genetic data based on specified criteria and generates informative Excel reports.
- Allows the choice of keeping intronic genes or excluding them from the analysis.

## Requirements

- Python 3.x
- Required Python libraries: `pandas`, `numpy`, `xlsxwriter`

## Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/genetic-analysis-tool.git
   ```

2. Navigate to the repository directory:

   ```bash
   cd genetic-analysis-tool
   ```

3. Run the script with appropriate command-line arguments to perform the desired analysis:

   ```bash
   python script.py <mode> [options] <input_files>
   ```

   Available modes:
   - `father_mother`: Analyze father-mother genetic data.
   - `mother_child`: Analyze mother-child genetic data.
   - `father_mother_child`: Analyze father-mother-child genetic data.

   Example usage:
   ```bash
   python script.py father_mother mother_genes.csv father_genes.csv mother_pathogen.csv father_pathogen.csv omim.txt
   ```

   Use the `--keep-intronic` flag to include intronic genes in the analysis.

## Output

The script generates Excel reports with different sections for each type of analysis, such as shared genes, compound genes, dangerous genes, and more. The reports provide insights into the genetic data for the specified family configuration.

## Configuration

The script provides various methods for filtering and analyzing genetic data, each with specific criteria. You can modify these methods to adjust the filtering criteria according to your analysis needs.








# Gene Filtering Script

This script is designed to filter and analyze genetic data files for specific patterns and relationships among genes. It utilizes multithreading for enhanced performance when processing large datasets.

## Prerequisites

Before using the script, ensure you have the following prerequisites in place:

- Python (3.6+)
- The `xlsxwriter` library (Install using `pip install XlsxWriter`)

## Usage

Run the script with the desired mode and options:

```bash
python script_name.py [options] <mode> [mode-specific arguments]
```

### Options

- `--keep-intronic`: Include intronic genes in the analysis.
- `--no-keep-intronic`: Exclude intronic genes from the analysis.

### Modes

1. `mother_child`: Filter the mother-child dataset.
   ```bash
   python script_name.py mother_child <path_to_mother_file> <path_to_child_file> [--keep-intronic]
   ```

2. `father_mother_child`: Filter the father-mother-child dataset.
   ```bash
   python script_name.py father_mother_child <path_to_father_file> <path_to_mother_file> <path_to_child_file> [--keep-intronic]
   ```

## Output

The script generates an Excel spreadsheet (.xlsx) containing the filtered and analyzed gene data. The resulting file will be named based on the child gene data file.

## Sample Run

To filter the mother-child dataset while excluding intronic genes:

```bash
python script_name.py mother_child mother_gene.csv child_gene.csv --no-keep-intronic
```

## Performance

The script takes advantage of multithreading to concurrently process data, which can significantly enhance performance for larger datasets.

## Execution Time

The script will print the execution time at the end of its run.




# Genetic Data Analysis Tool

This script is designed to analyze genetic data files for specific genetic variations (snips) and their associations with phenotypes. It utilizes multithreading for efficient data processing and generates frequency-based plots for visualizing the results.

## Prerequisites

Before using the script, ensure you have the following prerequisites in place:

- Python (3.6+)
- The `matplotlib` library (Install using `pip install matplotlib`)

## Usage

Run the script with the required arguments:

```bash
python script_name.py <metadata_file> <chromosomes_file> <chr> <pos> <ref> <alt> <phenotypes>
```

### Arguments

- `metadata_file`: Path to the metadata file.
- `chromosomes_file`: Path to the chromosomes file.
- `chr`: Snip chromosome type.
- `pos`: Snip position.
- `ref`: Snip reference type.
- `alt`: Snip alternative type.
- `phenotypes`: Phenotypes to search, separated by commas.

## Output

The script generates a PNG image that visualizes the frequency distribution of specified phenotypes associated with the given genetic variation. The image is saved using the format: `<chr>_<pos>_<ref>_<alt>.png`.

## Sample Run

To analyze genetic data with the provided metadata and chromosomes files, search for a specific snip, and visualize the frequency distribution of phenotypes:

```bash
python script_name.py TestData_metaData.txt chr11Data_test.txt.gz chr 12345 A T phenotype1,phenotype2
```

## Visualization

The script generates a frequency-based plot using the `matplotlib` library. The X-axis represents different phenotypes, while the Y-axis represents the frequency of the snip's association with those phenotypes.

## Notes

- The script utilizes multithreading for efficient data processing.
- It processes the provided metadata and chromosomes files to analyze genetic associations.



# Genetic Data Analysis GUI

This graphical user interface (GUI) application is designed to facilitate the analysis of genetic data files, particularly focusing on specific genetic variations (SNPs) and their associations with phenotypes. The application provides options to filter and visualize the data. It utilizes multithreading for efficient data processing and generates interactive plots for visualizing the results.

## Prerequisites

Before using the GUI, ensure you have the following prerequisites in place:

- Python (3.6+)
- The `tkinter` library (usually included with Python)
- The `matplotlib` library (Install using `pip install matplotlib`)

## Features

- Open and process meta data and chromosome files.
- Enter the chromosome, position, reference, alternative, and phenotypes for filtering.
- Visualize the frequency distribution of SNP associations with phenotypes.
- Interactive plot that allows hovering and clicking for detailed information.
- Save the generated plot as a PNG image.

## Usage

1. Run the script using the following command:

   ```bash
   python script_name.py
   ```

   Replace `script_name.py` with the actual name of your script.

2. The GUI window will appear, providing the following options:

   - Enter chromosome, position, reference, alternative, and phenotypes in the left panel.
   - Open meta data and chromosome files using the "Open meta data file" and "Open chromosome file" buttons.
   - Click the "Search" button to start the analysis. Progress will be displayed using a progress bar.
   - The generated plot will be displayed in the right panel, allowing interaction.

3. After processing, you can save the generated plot as a PNG image using the "Save plot as png" button.

## Interactive Plot

The generated plot allows interaction:

- Hovering over bars displays frequency information for each phenotype.
- Clicking on a bar shows detailed information about associated IDs and phenotypes.

## Notes

- The script utilizes multithreading for efficient data processing.
- It provides a user-friendly GUI for analyzing genetic data and generating interactive plots.

