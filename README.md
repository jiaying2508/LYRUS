# LYRUS
Python program for pathogenicity prediction of human missense variants.

This script is tested using Python3.7.6

## Running Instructions
Clone this repository and run the following command within the downloaded directory
```console
$ python inputFile outputDir fathmmFile
```

The inputFile should contain 2 column:
  1. UniProt ID
  2. Single amino acid variant: [aa_ref][aa_pos][aa_var]

Example input is:  
```
Q9NQZ7 V363G
P11245 E203D
Q6XZF7 R1101Q
B1AL17 A139V
Q9NTN9-2 R423H
Q92887 T486I
............
```

The outputDir should be a full path to the desired directory to store the outputs

The fathmmFile(full path) should contain the output from FATHMM. To get the FATHMM output,
go to http://fathmm.biocompute.org.uk/inherited.html and run using the same InputFile.

It might be necessary to manually install the DSSP program, for instance
by typing on Linux:
```console
$ sudo apt install dssp
```

## Install from source
Rhapsody is written in pure Python so no local compilation is needed.

To install all needed dependencies, we strongly suggest to use Conda and create
a new environment with:
```console
$ conda create -n rhapsody python=3 numpy scikit-learn requests pyparsing matplotlib biopython tqdm
$ conda activate rhapsody
$ pip install prody
$ conda install -c salilab dssp
```

After cloning/forking the Rhapsody repository, you can permanently add the
repository path to the conda environment with:
```console
$ conda develop path/to/local/repository
```

If not using Conda, you can manually install all dependencies and then add
the repository location to the `PYTHONPATH` environmental variable. For
example, on Linux simply add the following line to your `~/.bashrc`:
```console
export PYTHONPATH="path/to/local/repository/:$PYTHONPATH"
```

If you are running on Windows, please follow this
[tutorial](https://stackoverflow.com/a/4855685).

## Running initial setup

After installation, please run:
```console
import rhapsody as rd
rd.initialSetup()
```
