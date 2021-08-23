# LYRUS
Python program for pathogenicity prediction of human missense variants.

LYRUS is built on top of several existing Python libraries as well as other Software, and is tested using Python3.7.6
## Files included in LYRUS
- lyrus.py
  This is the script to run LYRUS
- train.csv
  Training file for the XGBoost Classifier
- gene.txt (optional)
- refseqs.txt (optional)

## Other data files
The **data** folder that includes pre-computed variation number and EVMutation score (using the same orthologs as the variation number; differs from the ones provided by the Marks Lab https://marks.hms.harvard.edu/evmutation/downloads.html) can be downloaded at https://drive.google.com/drive/folders/1bFMi78D4LqjGMDZiP_X6OzBBcsttSoSy?usp=sharing. If you decided to use the pre-computed scores, please put the **data** folder in the **LYRUS** directory.

## Required python packages
Python packages (most of which can be installed using pip) needed to run LYRUS include:
- skbio: http://scikit-bio.org
- pandas: https://pandas.pydata.org/docs/getting_started/install.html
- numpy: https://numpy.org/install/
- scipy: https://www.scipy.org/install.html
- xgboost: https://xgboost.readthedocs.io/en/latest/install.html
- sklearn: https://scikit-learn.org/stable/install.html
- Bio: https://biopython.org/wiki/Download
- BeautifulSoup: https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-beautiful-soup
- evcouplings: http://prody.csb.pitt.edu/downloads/
- prody: http://prody.csb.pitt.edu/downloads/
- rhapsody: http://rhapsody.csb.pitt.edu/download.php
- pyrosetta: https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download

## Required external packages
LYRUS also depends on the following external packages:
Install **command line** version for:
1. Clustal Omega: http://www.clustal.org/omega/
2. PAUP: http://phylosolutions.com/paup-test/

Install the following files and put it in the **LYRUS** directory:
1. plmc-master: https://github.com/debbiemarkslab/plmc
2. FoldX: http://foldxsuite.crg.eu
3. FreeSASA: https://freesasa.github.io
4. MAESTRO: https://pbwww.services.came.sbg.ac.at/?page_id=477
5. P2Rank: https://github.com/rdk/p2rank

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
