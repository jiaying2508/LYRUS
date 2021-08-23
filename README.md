# LYRUS
Python program for pathogenicity prediction of human missense variants.

LYRUS is built on top of several existing Python libraries as well as other Software, and is tested using Python3.7.6

## Required python packages
Python packages (most of which can be installed using pip) needed to run LYRUS include:
  *requests
  *skbio
  pandas
  numpy
  scipy
  xgboost
  sklearn
  Bio: https://biopython.org/wiki/Download
  BeautifulSoup: https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-beautiful-soup
  evcouplings: http://prody.csb.pitt.edu/downloads/
  prody: http://prody.csb.pitt.edu/downloads/
  rhapsody: http://rhapsody.csb.pitt.edu/download.php
  pyrosetta: https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download

The instruction to install **pyrosetta** can be found at https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download

## Running Instructions
Clone this repository and run the following command within the downloaded directory
```console
$ python inputFile outputDir fathmmFile
```

The **inputFile** should contain 2 column:
  1. UniProt ID
  2. Single amino acid variant: [aa_ref][aa_pos][aa_var]

Example **inputFile** contains:  
```
Q9NQZ7 V363G
P11245 E203D
Q6XZF7 R1101Q
B1AL17 A139V
Q9NTN9-2 R423H
Q92887 T486I
............
```

The **outputDir** should be a **full path** to the desired directory to store the outputs

The **fathmmFile**(full path) should contain the output from FATHMM. To get the FATHMM output,
go to http://fathmm.biocompute.org.uk/inherited.html and run using the **InputFile**.
