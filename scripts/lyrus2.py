import os
import sys
from lyrusClass import lyrusClass
import time

gene = sys.argv[1]
uniprot = sys.argv[2]
acc = sys.argv[3]
outputDir = sys.argv[4]
lyrusModel = lyrusClass(gene, uniprot, outputDir)
lyrusModel.setAccession(acc)
lyrusModel.getPDB()
#get fathmm input file
lyrusModel.getParameters(maestroDir='MAESTRO_linux_x64',p2rankDir='p2rank_2.2', EVparam=50)

