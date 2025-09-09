# `bigbugdata`

## Setup

```
$ git clone https://github.com/tombch/bigbugdata.git
$ cd bigbugdata/
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install .
```

## Usage

```
$ bigbugdata --help
usage: bigbugdata [-h] -r REPORTS [REPORTS ...] [-o OUTPUT] [-n SAMPLE GROUP] [-R RANK] [-t TOPHITS] [-v]

options:
  -h, --help            show this help message and exit
  -r REPORTS [REPORTS ...], --reports REPORTS [REPORTS ...]
                        Input KrakenUniq report files
  -o OUTPUT, --output OUTPUT
                        Directory to store the output files (default: results)
  -n SAMPLE GROUP, --nc-group SAMPLE GROUP
                        Provide REGEX patterns to match a negative control and its group of samples
  -R RANK, --rank RANK  Taxonomic rank to filter the reports by (default: species)
  -t TOPHITS, --tophits TOPHITS
                        Number of top hits to include in the tophits output (default: 15)
  -v, --version         Output the version of bigbugdata
```

### Example Usage

Lets consider an example naming system, where samples are split into four groups: `CF_DNA`, `CF_RNA`, `CP_DNA` and `CP_RNA`. There is a negative control for each group - these are found in the samples starting with the prefixes `CF_DNA_Negative`, `CF_RNA_Negative`, `CP_DNA_Negative` and `CP_RNA_Negative`.

Then an example `bigbugdata` command would be:

```
$ bigbugdata --reports /path/to/reports/dir --output /path/to/output/dir --nc-group CF_DNA_Negative CF_DNA --nc-group CF_RNA_Negative CF_RNA --nc-group CP_DNA_Negative CP_DNA --nc-group CP_RNA_Negative CP_RNA
```
