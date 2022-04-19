import os
import sys
from LYRUS.lyrusClass import lyrusClass, lyrusPredict
import time
from shutil import copyfile

def runLyrus1(gene, uniprot, outputDir):
    lyrusModel = lyrusClass(gene, uniprot, outputDir)
    if lyrusModel.acc != None:
        return lyrusModel
    try:
        acc = lyrusModel.getFasta()
    except:
        time.sleep(10)
        try:
            acc = lyrusModel.getFasta()
        except:
            time.sleep(10)
            try:
                acc = lyrusModel.getFasta()
            except:
                pass

    if acc == None or len(acc) == 0:
        raise TypeError('the accession number is None. Either provide the accession number or set skip to be False.')
    
    return lyrusModel

def runLyrus2(gene, uniprot, acc, outputDir):
    # get all parameters
    lyrusModel = lyrusClass(gene, uniprot, outputDir, savFile='{}/{}/input1.txt'.format(outputDir,gene))
    lyrusModel.setAccession(acc)
    lyrusModel.getPDB()
    #get fathmm input file
    lyrusModel.getParameters(maestroDir='MAESTRO_OSX_x64',p2rankDir='p2rank_2.2', EVparam=5)

def runLyrus3(outputDir, geneList):
    #get fathmm input
    file = open('{}/final.txt'.format(outputDir), 'w')
    outputFile = open('{}/fathmmInput.txt'.format(outputDir), 'w')
    for gene in geneList:
        print(gene)
        files = os.listdir('{}/{}'.format(outputDir, gene))
        if 'finalParam.csv' in files:
            with open('{}/{}/fathmmInput.txt'.format(outputDir, gene)) as f:
                outputFile.write(f.readline())
                outputFile.write('\n')
            file.write('{}\n'.format(gene))
    outputFile.close()
    file.close()

def runLyrus4(gene, fathmmFile, outputDir, uniprot):
    lyrusPredict(gene, fathmmFile, outputDir, uniprot)

def main():
    currDir = os.getcwd()
    # dataDir = '{}/data'.format(currDir)
    outputDir = '{}/{}'.format(currDir, sys.argv[2])
    try:
        os.mkdir(outputDir)
    except:
        pass
    existingGene = []
    accDict = {}
    uniprotDict = {}
    try:
        with open('{}/accession.txt'.format(outputDir)) as f:
            for line in f:
                line = line.rstrip().split(',')
                existingGene.append(line[0])
                accDict[line[0]] = line[2]
                uniprotDict[line[0]] = line[1]
    except:
        pass
    
    if sys.argv[3] == '1':
        geneFile = sys.argv[1]
        with open(geneFile) as f:
            for line in f:
                line = line.rstrip().split()
                gene = line[0]
                uniprot = line[1]
                
                if gene not in existingGene:
                    print(gene)
                    file = open('{}/accession.txt'.format(outputDir), 'a')
                    errorFile = open('{}/error.txt'.format(outputDir), 'a')
                    try:
                        lyrusModel = runLyrus1(gene, uniprot, outputDir)
                        if lyrusModel.acc != None:
                            file.write('{},{},{}\n'.format(gene,uniprot,lyrusModel.acc))
                        else:
                            errorFile.write('{}\n'.format(gene))
                    except:
                        errorFile.write('{}\n'.format(gene))
                    errorFile.close()
                    file.close()
    elif sys.argv[3] == '2':
        for gene in accDict:
            runLyrus2(gene, uniprotDict[gene], accDict[gene], outputDir)

            # try:
            #     runLyrus2(gene, uniprotDict[gene], accDict[gene], outputDir)
            # except:
            #     pass
    elif sys.argv[3] == '3':
        runLyrus3(outputDir, existingGene)
    else:
        #sys.argv[4] is fathmm file name
        genes = []
        with open('{}/final.txt'.format(outputDir)) as f:
            for line in f:
                line = line.rstrip()
                genes.append(line)
        for gene in genes:
            runLyrus4(gene, sys.argv[3], outputDir, uniprotDict[gene])

if __name__ == '__main__':
    main()

