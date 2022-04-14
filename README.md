# LYRUS: A Machine Learning Model for Predicting the Pathogenicity of Missense Variants
LYRUS incorporates five sequence-based, six structure-based, and four dynamics-based features. Uniquely, LYRUS includes a newly-proposed sequence co-evolution feature called variation number. LYRUS was trained using a dataset that contains 4,363 protein structures corresponding to 22,639 SAVs from the ClinVar database.

LYRUS is built on top of several existing Python libraries as well as other Software, and is tested using Python3.7.4

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
LYRUS also depends on the following external packages:<br/><br/>
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
Clone this repository and run the following command within the downloaded directory, with python version 3.7.4 or higher.

```console
import os
from LYRUS.lyrusClass import lyrusClass, lyrusPredict

gene = 'A1BG'
uniprot = 'P04217'
currDir = os.getcwd()
outputDir = '{}/test'.format(currDir)
try:
    os.mkdir(outputDir)
except:
    print('Output directory already exist')

#load model
lyrusModel = lyrusClass(gene, uniprot, outputDir, savFile=None)

#download orthologs from NCBI
lyrusModel.getFasta()

#download PDB from SWISS-MODEL
lyrusModel.getPDB()

#calculate all the parameters except for fathmm
lyrusModel.getParameters(maestroDir='MAESTRO_OSX_x64',p2rankDir='p2rank_2.2')
```

The **fathmmFile** should contain the output from FATHMM. To get the FATHMM output,
go to http://fathmm.biocompute.org.uk/inherited.html and run using the **fathmmInput.txt** available in the output directory.

```console
fathmmFile = 'test/fathmm.txt'

#calculate lyrus probability
lyrusPredict(gene, fathmmFile, outputDir, uniprot)
```