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
