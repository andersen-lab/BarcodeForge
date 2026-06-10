# BarcodeForge

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15654220.svg)](https://doi.org/10.5281/zenodo.15654220)

A CLI tool for generating pathogen-specific barcodes for Freyja.

## Installation

BarcodeForge requires **Python >= 3.10** (tested on 3.10–3.13).

### Conda (recommended)

BarcodeForge is available on [Bioconda](https://anaconda.org/bioconda/barcodeforge):

```
conda install bioconda::barcodeforge
```

### From source (GitHub)

Install the latest version directly from the GitHub repository:

```
pip install git+https://github.com/gp201/BarcodeForge.git
```

Or clone the repository and install it locally:

```
git clone https://github.com/gp201/BarcodeForge.git
cd BarcodeForge
pip install .
```

For development (includes test and formatting tools):

```
pip install -e ".[dev-requirements,test-requirements]"
```

### Additional dependencies

The following tools are required and must be installed separately via Conda
(from the `bioconda` and `conda-forge` channels):

```
conda install -c conda-forge -c bioconda usher ucsc-fatovcf
```

## Usage

Once installed, run the tool with:

```
barcodeforge --help
```

For example, to generate a barcode:

```
barcodeforge barcode reference.fasta aligned.fasta tree.nwk lineages.tsv --threads 4 --prefix RSVa
```
