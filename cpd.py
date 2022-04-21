from collections import defaultdict
from LYRUS.lyrusClass import lyrusClass, lyrusPredict
import os

savDict = defaultdict(list)
uniprotDict = {}
with open('potentialCS.csv') as f:
    for line in f:
        line = line.rstrip().split(',')
        if int(line[4]) > 1:
            uniprotDict[line[0]] = line[1]
            savDict[line[0]].append(line[2])
            savDict[line[0]].append(line[3][-1]+line[3][:-1]+line[5])
# print(savDict)

outputDir = '/Users/jiayinglai/LYRUS/cpd'
for gene in savDict:
    savDict[gene] = list(set(savDict[gene]))
    uniprot = uniprotDict[gene]
    dir = '{}/{}'.format(outputDir, gene)
    try:
        os.mkdir(dir)
    except:
        pass
    with open('{}/input.txt'.format(dir), 'w') as f:
        for sav in savDict[gene]:
            f.write('{}\n'.format(sav))
    lyrusModel = lyrusClass(gene, uniprot, outputDir, savFile='{}/{}/input.txt'.format(outputDir,gene))
    lyrusModel.getFasta()
    lyrusModel.getPDB()
    lyrusModel.getParameters(maestroDir='MAESTRO_OSX_x64',p2rankDir='p2rank_2.2', EVparam=5, dataDir='/Users/jiayinglai/LYRUS_data/LYRUS_local/data')
    fathmmFile = 'cpd_fathmm.txt'
    lyrusPredict(gene, fathmmFile, outputDir, uniprot)