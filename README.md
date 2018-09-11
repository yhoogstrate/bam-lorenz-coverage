bam-lorenz-coverage
===================

Generate Lorenz plots and Coverage plots directly from BAM files

Implemented in:
 * Python3 + Matplotlib + Pysam


Installation:

```
$ git clone https://github.com/yhoogstrate/bam-lorenz-coverage.git
$ cd bam-lorenz-coverage
$ virtualenv -p python3 .venv
$ source .venv/bin/activate
$ python setup.py install

$ bam-lorenz-coverage --help
```

Usage:

```
Usage: bam-lorenz-coverage [OPTIONS] INPUT_ALIGNMENT_FILE

Options:
  --version                  Show the version and exit.
  -l, --lorenz-table TEXT    Output table Lorenz-curve (for stdout use: -)
  -x, --roc                  Output Lorenz-curve ROC to $lorenz_table.roc.txt
                             [ requires --lorenz-table to be set to file ]
  -c, --coverage-table TEXT  Output table Coverage-graph (for stdout use: -)
  -L, --lorenz-svg TEXT      Output figure Lorenz-curve (SVG).
  -C, --coverage-svg TEXT    Output figure Coverage-graph (SVG).
  --help                     Show this message and exit.

```

The lowercase arguments (-l, -c) allow extraction of the raw data tables for custom plotting. The uppercase arguments (-L, -C) directly generate a plot. The implemented plot only contains one sample per plot. For multi-sample plots, use the column tables and your imagination.