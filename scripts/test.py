import os
from LYRUS.lyrusClass import lyrusClass, lyrusPredict

gene = 'A1BG'
uniprot = 'P04217'
currDir = os.getcwd()
outputDir = '{}/test1'.format(currDir)
try:
    os.mkdir(outputDir)
except:
    print('Output directory already exist')

#load model
lyrusModel = lyrusClass(gene, uniprot, outputDir, savFile='{}/{}/input.txt'.format(outputDir,gene))

#download orthologs
lyrusModel.getFasta()

#download PDB
lyrusModel.getPDB()

#calculate all the parameters except for fathmm
lyrusModel.getParameters(maestroDir='MAESTRO_OSX_x64',p2rankDir='p2rank_2.2', EVparam=5)

#run fathmm at http://fathmm.biocompute.org.uk/inherited.html and save results in fathmm file
fathmmFile = 'test1/fathmm.txt'

lyrusPredict(gene, fathmmFile, outputDir, uniprot)

