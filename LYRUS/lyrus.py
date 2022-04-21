#run XGBoost Classifier to predict SAV disease pathogenicity
'''
python lyrus.py -i <inputFile> -o <outputDir>; need to be full path> -f <fathmmFile>

python version:
    python3.7.6

args:
    -i inputFile:
        uniprot followed by single amino acid variant, separate by space
        example:
            Q9NQZ7 V363G

    -o outputDir:
        full path of the desired output directory

    -f fathmmFile:
        use the same inputFile to get fathmm output and provide the fathmm output file name here; use full path 
        FATHMM can be run at http://fathmm.biocompute.org.uk/inherited.html

required python packages:
    requests
    Bio
    skbio
    pandas
    BeautifulSoup
    numpy
    scipy
    prody
    pyrosetta
    evcouplings
    rhapsody
    xgboost
    sklearn

required external packages:
    plmc-master
    foldx
    freesasa
    maestro
    p2rank
    clustal omega
    paup

skipped process if file exist:
    ev_ and vn_ in data
    orthologs files

required files:
    refseqs.txt
    gene.txt
    train.csv
    data from https://drive.google.com/drive/folders/1bFMi78D4LqjGMDZiP_X6OzBBcsttSoSy?usp=sharing
    
'''

__author__ = "Jiaying Lai"
__date__ = "August, 2021"
__maintainer__ = "Jiaying Lai"
__email__ = "jiaying_lai@brown.edu"

import os
import re
import sys
import requests
import time
import shutil
from shutil import copyfile
from Bio import Entrez
from skbio import TreeNode
from io import StringIO
from Bio import Phylo
from Bio import SeqIO
import numpy as np
from scipy import stats
import random
import math

import json
from bs4 import BeautifulSoup
import argparse
import urllib.parse
import urllib.request
import multiprocessing as mp
import prody
from pyrosetta import *
from pyrosetta import PyMOLMover
from pyrosetta.toolbox import cleanATOM
from pyrosetta.teaching import *
from pyrosetta.toolbox import mutate_residue
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.core.fragment import *
from pyrosetta.rosetta.protocols.moves import *
from pyrosetta.rosetta.protocols.rigid import *
from pyrosetta.rosetta.protocols.docking import *

import pandas as pd
from evcouplings.couplings import CouplingsModel
from evcouplings.mutate import predict_mutation_table, single_mutant_matrix
import rhapsody as rd

from xgboost import XGBClassifier
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from scipy.stats import spearmanr, pearsonr

'''
run imputation for missing feature values
'''
def runImpute(file, outputDir, inputFile):
	f = open(file,'r')
	f.readline()
	matrix = []

	for line in f:
		line=line.strip().split(',')
		matrix.append(line[4:])
	f.close()
	matrix = [[float(i) if i!='nan' else np.nan for i in row]for row in matrix]
    
	imp = IterativeImputer(max_iter=100, random_state=0, n_nearest_features=15)
	matrix=imp.fit_transform(matrix)
	np.savetxt('{}/{}_LYRUS_imputed.csv'.format(outputDir, inputFile.split('.')[0]),matrix,delimiter=',')

'''
download gene names using uniprot id from uniprot
'''
def getGeneName(uniprotList, outputDir):

    print('\n# getting gene names\n')

    url = 'https://www.uniprot.org/uploadlists/'

    query = []
    for uniprot in uniprotList:
        if '-' in uniprot:
            query.append(uniprot.split('-')[0])
        else:
            query.append(uniprot)

    query = ' '.join(query)
    params = {
    'from': 'ACC+ID',
    'to': 'GENENAME',
    'format': 'tab',
    'query': query
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()

    output = open('{}/gene1.txt'.format(outputDir), 'w')
    output.write(response.decode('utf-8'))
    output.close()

'''
download pdb file
'''
def getPDB(uniprot, pdbDir):
    print('# downloading {} from SwissModel'.format(uniprot))
    upfilename = uniprot+'.pdb'
    url = 'https://swissmodel.expasy.org/repository/uniprot/'+upfilename+'?provider=swissmodel'
    r = requests.get(url, allow_redirects=True)
    open('{}/{}'.format(pdbDir, upfilename), 'wb').write(r.content)
    uppdb_unclean = open('{}/{}'.format(pdbDir, upfilename))
    uppdb_clean = open('{}/{}_clean.pdb'.format(pdbDir, uniprot),'w')
    for atom in uppdb_unclean:
        atom = atom.rstrip()
        if atom[:4] == 'ATOM':
            uppdb_clean.write(atom)
            uppdb_clean.write('\n')
    uppdb_unclean.close()
    uppdb_clean.close()

'''
get pdb matched single chain
'''
def getSingleChain(uniprot, pdbDir, savDict, logFile):
    AA = {'C':'CYS','D':'ASP','S':'SER','Q':'GLN','K':'LYS','I':'ILE','P':'PRO','T':'THR','F':'PHE','N':'ASN','G':'GLY','H':'HIS','L':'LEU','R':'ARG','W':'TRP','A':'ALA','V':'VAL','E':'GLU','Y':'TYR','M':'MET'}
    wt = AA[savDict[uniprot][0][0].upper()]
    res = savDict[uniprot][0][1:-1]
    oldfile=open('{}/{}_clean.pdb'.format(pdbDir, uniprot),'r')
    
    found = False
    chain = ''
    for line in oldfile:
        line=line.strip()
        l = [char for char in line]
        residue = ''.join(l[17:20]).strip()
        chain_x = l[21].strip()
        pos = ''.join(l[22:26]).strip()

        if residue == wt and pos == res:
            found = True
            chain = chain_x
            break
        
    oldfile.close()
    
    if not found:
        lFile = open(logFile, 'a')
        lFile.write('uniprot {} error, SAV {} not match\n'.format(uniprot, savDict[uniprot][0]))
        print('\n##############uniprot {} error, SAV {} not match#############\n'.format(uniprot, savDict[uniprot][0]))
        lFile.close()

    oldfile = open('{}/{}_clean.pdb'.format(pdbDir, uniprot),'r')
    newfile = open('{}/{}_clean_single.pdb'.format(pdbDir, uniprot),'w')
    for line in oldfile:
        line = line.strip()
        l = [char for char in line]
        residue = ''.join(l[17:20]).strip()
        chain_x = l[21].strip()
        pos = ''.join(l[22:26]).strip()

        temp = line.split()
        if chain_x == chain:
            newfile.write(line)
            newfile.write('\n')

    oldfile.close()
    newfile.close()        

'''
replace the amino acid with X if it is not among the 20 amino acid codes
'''
def replaceAA(sequence):
    seq = ''
    for char in sequence:
        if char in 'ARNDCQEGHILKMFPSTWYV':
            seq = seq + char
        else:
            seq = seq + 'X'
    return seq

'''
write a sequence to fasta file
'''
def write_fasta(file, accession, sequence):
    file.write('>{}\n'.format(accession))
    if len(sequence) <= 70:
        file.write('>{}\n'.format(sequence))
    else:
        line_num_exact = len(sequence) / 70
        line_num = math.floor(line_num_exact)
        for i in range(line_num):
            start = i * 70
            stop = i * 70 + 70
            file.write('{}\n'.format(sequence[start:stop]))

        start = line_num * 70
        file.write('{}\n'.format(sequence[start:]))

def writeNexus(fileName, seqDict, sequencetype):
    nexus_file = open(fileName, 'w')
    nexus_file.write('#NEXUS\n\n')
    nexus_file.write('BEGIN TAXA;\n')
    nexus_file.write('\tTITLE Taxa;\n')
    nexus_file.write('\tDIMENSIONS NTAX={};\n'.format(len(seqDict)))

    taxa = ' '.join(sorted(seqDict.keys()))
    dimension = 0
    for key in seqDict.keys():
        dimension = len(seqDict[key])
        break

    nexus_file.write('\tTAXLABELS\n')
    nexus_file.write('\t\t{}\n\t;\n\n'.format(taxa))
    nexus_file.write('END;\n\n')

    nexus_file.write('BEGIN CHARACTERS;\n')
    nexus_file.write('\tTITLE Character_Matrix;\n')
    nexus_file.write('\tDIMENSIONS NCHAR={};\n'.format(dimension))
    nexus_file.write('\tFORMAT DATATYPE = {} GAP = - MISSING = ?;\n'.format(sequencetype))
    nexus_file.write('\tMATRIX\n')

    for key in sorted(seqDict):
        nexus_file.write('\t{} {}\n'.format(key, seqDict[key]))

    nexus_file.write('\n;\nEND;')
    nexus_file.close()

def getTreeCMD(nexusFileName, outputFile, dir):
    #TREE SEARCH METHOD
    #tree search method for simultaneous analysis tree & support tests
    treeSearchMethod = 'hsearch nreps=1000 swap=tbr multrees=no'

    #BOOTSTRAP TREE SEARCH METHOD
    #   specify method for bootstrap tree searhc
    #   bootstrapMethod = "search=bandb"
    bootstrapMethod = 'search=heuristic nreps=50'

    createTreeCmdFile = open(outputFile, 'w')
    createTreeCmdFile.write('#NEXUS\n\n')
    createTreeCmdFile.write('set warnReset = no;\n')
    createTreeCmdFile.write('set increase = auto;\n')
    createTreeCmdFile.write('set datastorage=full;\n')
    createTreeCmdFile.write('set criterion=parsimony;\n')
    createTreeCmdFile.write('execute {};\n'.format(nexusFileName))
    createTreeCmdFile.write('{};\n'.format(treeSearchMethod))
    createTreeCmdFile.write('filter best;\n')

    treeFile = outputFile.replace('_getTree.cmd', '_trees.nex')

    createTreeCmdFile.write('savetrees file={} format=nexus replace=yes root=yes;\n'.format(treeFile))
    createTreeCmdFile.write('quit warnTsave=no;\n')
    createTreeCmdFile.close()

'''
functions to calculate variation number
'''
def findSet(l, i, seqDict):
    ll = []
    for j in l:
        ll.append(str(seqDict[j][i]))
    ll = list(set(ll))
    return ll

def updateVN(node, child, variation_number, seqDict, length):

    allClades = re.findall(r'[a-zA-Z0-9]+', str(node))
    for c in child:
        cClades = re.findall(r'[a-zA-Z0-9]+', str(c))
        subClades = list(set(allClades) - set(cClades))
        for i in range(length):
            cSet = findSet(cClades, i, seqDict)
            subSet = findSet(subClades, i, seqDict)
            for item in cSet:
                if item not in subSet:
                    variation_number[i] = variation_number[i] + 1
    return variation_number

def generateVN(tree, seqDict, seqLength):
    variation_number = np.zeros((seqLength,), dtype=int)
    queue = []
    queue.append(tree)

    while len(queue) > 0:
        node = queue.pop()
        if node.is_tip():
            continue
        else:
            child = node.children
            variation_number = updateVN(node, child, variation_number, seqDict, seqLength)
            #print('variation')
            #print(variation_number)
            for c in child:
                queue.append(c)
    return variation_number

def processVN(index, inputDir, outputDir, file):
    if os.path.isfile('{}/vn_{}.txt'.format(outputDir, file)) and os.path.isfile('{}/evmutation_{}_result.txt'.format(outputDir, file)):
        return (index, -1)

    accession_full = file
    
    # os.system('julia evolutionAnalysis.jl {} {} {} protein {}'.format(inputDir, file, accession_full, outputDir))
    '''
    the comment out portion below is equivalent to the python version of evolutionAnalysis.jl
    '''
    #'''
    ################################################################################
    #process homologous protein sequence file
    ################################################################################
    try:
        in_file = open('{}/{}'.format(inputDir, file, 'r'))

        print('# Generating Variation Number for {}'.format(file))

        accession = ''
        sequence = ''
        fastaDict = {}
        for line in in_file:
            line = line.rstrip()
            if len(line) == 0:
                continue

            if '>' in line:
                #if the first accession, initiate
                if accession == '':
                    accession = re.search('>([A-Za-z0-9_]+.[0-9]+)', line)
                    accession = accession.group(0)[1:]
                    continue

                #store accession-sequence pair in fastaDict
                sequence = replaceAA(sequence)
                fastaDict[accession] = sequence
                sequence = ''
                accession = re.search('>([A-Za-z0-9_]+.[0-9]+)', line)
                accession = accession.group(0)[1:]

            else:
                sequence = sequence + line

        sequence = replaceAA(sequence)
        fastaDict[accession] = sequence
        in_file.close()

        ################################################################################
        #         convert the homologous sequence file into fasta format
        ################################################################################
        file1 = open('{}/{}.fasta'.format(outputDir, accession_full), 'w')
        for key in fastaDict.keys():
            write_fasta(file1, key, fastaDict[key])
        file1.close()
        fastaDict.clear()

        ################################################################################
        #                      run cluster omega on the fasta file
        ################################################################################
        os.system('clustalo -i {}/{}.fasta -o {}/{}_aligned.fasta \
            --auto -v --force >/dev/null'.format(outputDir, accession_full, outputDir, accession_full))
            
        #'''
        ################################################################################
        #  read the aligned fasta file, and clean the identifier
        ################################################################################
        file1 = open('{}/{}_aligned.fasta'.format(outputDir, accession_full), 'r')
        seqDict = {}
        accession = ''
        seq = ''
        for line in file1:

            line = line.rstrip()
            if len(line) == 0:
                continue

            if '>' in line:
                if len(accession) > 0:
                    seqDict[accession] = seq
                    seq = ''

                accession = line.replace('_', '')
                accession = accession.replace('.', '')
                accession = accession.replace('>', '')

            else:
                seq = seq + line

        seqDict[accession] = seq
        file1.close()

        ################################################################################
        #                      convert fasta to nexus file
        ################################################################################
        nexusFile = '{}/{}.nex'.format(outputDir, accession_full)
        writeNexus(nexusFile, seqDict, 'protein')

        ################################################################################
        #                                   run paup
        ################################################################################
        getTreeCMD(nexusFile, '{}/{}_getTree.cmd'.format(outputDir, accession_full), outputDir)
        os.system('paup {}/{}_getTree.cmd >/dev/null'.format(outputDir, accession_full))

        ################################################################################
        #                              Generate variation number
        ################################################################################
        with open('{}/{}_aligned.fasta'.format(outputDir,accession_full)) as f:
            alist = [line.rstrip() for line in f]
        seqDict = {}

        accession = ''
        seq = ''
        for line in alist:
            if '>' in line:
                if accession != '':
                    seqDict[accession] = seq
                accession = line.replace('>', '')
                accession = accession.replace('_', '')
                accession = accession.replace('.', '')
                seq = ''
            else:
                seq = seq + line
        seqDict[accession] = seq

        homoAccession = accession_full
        homoAccession = homoAccession.replace('_', '')
        homoAccession = homoAccession.replace('.', '')
        homoSeq = seqDict[homoAccession]
        seqLength = len(homoSeq)

        ################################################################################
        #                       convert phylogenetic tree
        ################################################################################
        Phylo.convert('{}/{}_trees.nex'.format(outputDir, accession_full), 'nexus', '{}/{}_tree.tree'.format(outputDir, accession_full), 'newick')

        f = open('{}/{}_tree.tree'.format(outputDir, accession_full), 'r')
        tree = f.readline()

        tree = re.sub(r':\d.\d+', '', tree)
        tree = TreeNode.read(StringIO(tree))

        variation_number = generateVN(tree, seqDict, seqLength)

        homoIndexList = []
        f_vn = []
        for i in range(len(homoSeq)):
            if str(homoSeq[i]) != '-':
                homoIndexList.append(i)
                f_vn.append(variation_number[i])

        outputFile = open("{}/vn_{}.txt".format(outputDir,accession_full), 'w')

        vn_max = max(f_vn)
        vn_min = min(f_vn)

        for i in range(len(homoIndexList)):
            j = i + 1
            vn = variation_number[homoIndexList[i]]
            vn = (vn - vn_min) / (vn_max - vn_min)
            outputFile.write('{}\t{}\t{}\n'.format(str(j), homoSeq[homoIndexList[i]], str(vn)))

        outputFile.close()

        ################################################################################
        #             generate evmutation input file
        ################################################################################
        outputFile = open("{}/evmutation_{}.txt".format(outputDir,accession_full), "w")
        seq = seqDict[homoAccession]
        seq1 = ''
        for i in homoIndexList:
            seq1 = seq1 + seq[i]
        write_fasta(outputFile, homoAccession, seq1)

        for key in seqDict.keys():
            if key == homoAccession:
                continue
            seq = seqDict[key]

            seq1 = ''
            for i in homoIndexList:
                seq1 = seq1 + seq[i]
            write_fasta(outputFile, key, seq1)
        outputFile.close()

        return (index, homoIndexList)
    except:
        return (index, -1)

def processEVMutation(workingDir, file, homoIndexList):
    proteinAccession = file
    if homoIndexList == -1:
        return
    ################################################################################
    #                       calculate evmutation
    ################################################################################
    try:
        existingFiles = os.listdir(workingDir)
        evmutationfile = '{}.params'.format(proteinAccession)

        ll = (len(homoIndexList) - 1) * 0.2
        print('# performaing EVMutation Analysis for {}'.format(proteinAccession))

        if evmutationfile in existingFiles:
            os.system('plmc-master/bin/plmc -o {}/{}.params -le {} -lh 0.01 -m 100 {}/evmutation_{}.txt &>/dev/null'.format(workingDir, proteinAccession, ll, workingDir, proteinAccession))
        else:
            os.system('plmc-master/bin/plmc -o {}/{}.params -le {} -lh 0.01 -m 100 {}/evmutation_{}.txt &>/dev/null'.format(workingDir, proteinAccession, ll, workingDir, proteinAccession))
            # pass
        c = CouplingsModel('{}/{}.params'.format(workingDir, proteinAccession))
        singles = single_mutant_matrix(c, output_column='effect_prediction_epistatic')

        pd.DataFrame(singles).to_csv('{}/evmutation_{}_result.txt'.format(workingDir, proteinAccession))
    except:
        pass

def getFasta(refseqID, var, refseqL, refseqDir):
    print('\n# process {} refseq'.format(refseqID))
    page = requests.get('https://www.ncbi.nlm.nih.gov/protein/?term={}'.format(refseqID))
    geneID = ''
    homo_protein = ''
    if page.status_code == 200:
        soup = BeautifulSoup(page.text, 'html.parser')
        soup.prettify()

        for line in soup:
            geneID1 = re.search('<li>\s?Gene\s?ID\s?:\s?(\d+)\s?<\/li>', str(line))
            if geneID1:
                geneID = geneID1.group(1)
                break
        time.sleep(0.34)
        page = requests.get('https://www.ncbi.nlm.nih.gov/gene/{}/ortholog/?term={}'.format(geneID,refseqID))
        if page.status_code == 200:
            refseqList = []
            soup = BeautifulSoup(page.text, 'html.parser')
            soup.prettify()
            c = 0
            for line in soup:
                if 'appData' in str(line):
                    match = re.search('appData\.genes\s?=\s?([^;]+);', str(line))
                    if match:
                        genes = match.group(1)
                        try:
                            res = json.loads(genes)
                        except:
                            match = re.search('appData\.genes\s?=\s?(.+}]);', str(line))
                            genes = match.group(1)
                            res = json.loads(genes)
                        count = 0
                        for gene in res:
                            count += 1
                            if gene['tax_id'] == 9606:
                                for item in gene['refseq_accessions']:
                                    try:
                                        homo_p = item['protein_acc']
                                        homo_t= item['transcript_acc']
                                    except:
                                        continue
                                    handle = Entrez.efetch(db="protein", id=homo_p, rettype = 'fasta', retmode='text')
                                    txt = handle.read()
                                    handle.close()
                                    seq = ''
                                    lCount = 0
                                    for l in txt.split('\n'):
                                        if lCount == 0:
                                            lCount = 1
                                            continue
                                        seq = seq + l
                                    pp = int(var[1:-1]) - 1
                                    try:
                                        if seq[pp] == var[0]:
                                            print('# match found', homo_p)
                                            homo_protein = homo_p
                                            refseqList.append(homo_protein)
                                            if homo_protein in refseqL:
                                                print('# {} already exist'.format(homo_p))
                                                return homo_protein
                                            break
                                    except:
                                        pass
                            else:
                                try:
                                    refseqList.append(gene['refseq_accessions'][0]['protein_acc'])
                                except:
                                    pass

            if len(homo_protein) == 0:
                print('not match found for gene {} SAV {}'.format(refseqID, var))
                return ''

            outputFile = open('{}/{}'.format(refseqDir, homo_protein), 'w')
            print('# downloading orthologs files')
            for gene in refseqList:
                time.sleep(0.34)
                handle = Entrez.efetch(db="protein", id=gene, rettype = 'fasta', retmode='text')
                txt = handle.read()
                handle.close()
                seq = ''
                lCount = 0
                for l in txt.split('\n'):
                    if lCount == 0:
                        lCount = 1
                        continue
                    seq = seq + l
                outputFile.write('>{}\n{}\n'.format(gene, seq))
            outputFile.close()
            return homo_protein
    return ''

def collect_result(result):
    global results
    results.append(result)

def runFreeSASA(freeSASADir, pdbDir, uniprot):
    freeSASAfile = '{}/freesasa_{}.txt'.format(freeSASADir,uniprot)
    pdbfile = '{}/{}_clean_single.pdb'.format(pdbDir, uniprot)
    os.system('freesasa --format=seq --shrake-rupley -n 200 --probe-radius 1.2 --n-threads 4 {} > {}'.format(pdbfile, freeSASAfile))

def runFoldX(uniprot, workingDir, pdbDir, currDir):
    pdbfile = '{}/{}_clean_single.pdb'.format(pdbDir, uniprot)
    existingFiles = os.listdir(workingDir)
    repairPDBfile = '{}_clean_single_Repair.pdb'.format(uniprot)
    
    if repairPDBfile in existingFiles:
        '''
        uncomment below
        '''
        print('# running foldX for {}'.format(uniprot))
        copyfile(pdbfile, '{}/{}_clean_single.pdb'.format(currDir, uniprot))
        print('# Using FoldX to generate {}'.format(repairPDBfile))
        os.system('./foldx --command=RepairPDB --pdb={}_clean_single.pdb --output-dir={} &>/dev/null'.format(uniprot,workingDir))
        os.remove('{}/{}_clean_single.pdb'.format(currDir, uniprot))
        pass
    else:
        print('# running foldX for {}'.format(uniprot))
        copyfile(pdbfile, '{}/{}_clean_single.pdb'.format(currDir, uniprot))
        print('# Using FoldX to generate {}'.format(repairPDBfile))
        os.system('./foldx --command=RepairPDB --pdb={}_clean_single.pdb --output-dir={} &>/dev/null'.format(uniprot,workingDir))
        os.remove('{}/{}_clean_single.pdb'.format(currDir, uniprot))
    print('# Using FoldX to calculate free energy {}'.format(uniprot))
    foldx_out = '{}/Dif_{}_clean_Repair.fxout'.format(workingDir,uniprot)

    try:
        os.remove(foldx_out)
    except:
        pass

    copyfile('{}/{}'.format(workingDir, repairPDBfile), '{}/{}'.format(currDir, repairPDBfile))
    os.system('./foldx --command=BuildModel --pdb={} --mutant-file={}/individual_list_{}.txt --output-dir={} --numberOfRuns=5 --out-pdb=false &>/dev/null'.format(repairPDBfile,workingDir,uniprot, workingDir))
    os.remove('{}/{}'.format(currDir, repairPDBfile))

def runMaestro(sav, maestroDir, chainDict, pdbDir):
    global maestro
    k = sav.split(',')
    uniprot = k[0]
    pos = k[1]
    wt = k[2]
    var = k[3]
    pdbfile = '{}/{}_clean_single.pdb'.format(pdbDir, uniprot)
    try:
        print('# running MAESTRO for {}'.format(sav))
        mut = '{}{}.{}{{{}}}'.format(wt,pos,chainDict[uniprot],var)
        maestro_file = '{}/{}_{}{}{}_maestro.txt'.format(maestroDir, uniprot, wt, pos, var)
        os.system('{}/maestro {}/config.xml {} --evalmut=\'{}\' --bu > {}'.format(maestro, maestro, pdbfile, mut, maestro_file))
    except:
        pass

def runPrody(uniprot, pdbDir, prodyDir):
    print('# running Prody for {}'.format(uniprot))
    '''
    ANM
    '''
    filename = '{}/{}_clean_single.pdb'.format(pdbDir, uniprot)
    prot = prody.parsePDB(filename)
    anm,atoms = prody.calcANM(prot, selstr='calpha')
    flucts = prody.calcSqFlucts(anm)

    pdbIndexDict = {}
    file = open(filename)
    index = -1
    currPos = 0
    count = 0
    
    for line in file:
        line = line.strip()
        l = [char for char in line]
        pos = ''.join(l[22:26]).strip()
        line = line.split()
        if count == 0:
            index = 0
            currPos = int(pos)
            pdbIndexDict[currPos] = index
            count = 1
            continue
        if currPos == int(pos):
            continue
        index += 1
        currPos = int(pos)
        pdbIndexDict[currPos] = index
    file.close()

    file1 = open('{}/{}_anm.txt'.format(prodyDir, uniprot), 'w')
    for pos in pdbIndexDict.keys():
        file1.write('{},{}\n'.format(pos, flucts[pdbIndexDict[pos]]))
    file1.close()

    #mechanical_stiffness
    gfp = prody.parsePDB(filename)
    calphas = gfp.ca
    anm = prody.ANM('ANM analysis')
    anm.buildHessian(calphas, cutoff=13.0)
    anm.calcModes(n_modes='all')
    stiffness = prody.calcMechStiff(anm, calphas)
    ms = []
    for row in stiffness:
        ms.append(str(sum(row)/len(row)))

    file1 = open('{}/{}_ms.txt'.format(prodyDir, uniprot), 'w')
    for pos in pdbIndexDict.keys():
        file1.write('{},{}\n'.format(pos, ms[pdbIndexDict[pos]]))
    file1.close()

    #Effectieness & sensitivity
    ampar_ca = prody.parsePDB(filename, subset='ca')

    anm_ampar = prody.ANM('AMPAR')
    anm_ampar.buildHessian(ampar_ca)
    anm_ampar.calcModes()
    prs_mat, eff, sen = prody.calcPerturbResponse(anm_ampar)
    eff = eff.tolist()
    sen = sen.tolist()

    file1 = open('{}/{}_effectiveness.txt'.format(prodyDir, uniprot), 'w')
    for pos in pdbIndexDict.keys():
        file1.write('{},{}\n'.format(pos, eff[pdbIndexDict[pos]]))
    file1.close()

    file1 = open('{}/{}_sensitivity.txt'.format(prodyDir, uniprot), 'w')
    for pos in pdbIndexDict.keys():
        file1.write('{},{}\n'.format(pos, sen[pdbIndexDict[pos]]))
    file1.close()

def processProdyFile(file, pos):
    try:
        f = open(file)
        for line in f:
            line = line.rstrip()
            l = line.split(',')
            if l[0] == pos:
                return float(l[1])
        f.close()
    except:
        return 'nan'
    return 'nan'

def getPrody(key, prodyDir):
    l = []
    uniprot = key.split(',')[0]
    pos = key.split(',')[1]
    l.append(processProdyFile('{}/{}_anm.txt'.format(prodyDir, uniprot), pos))
    l.append(processProdyFile('{}/{}_ms.txt'.format(prodyDir, uniprot), pos))
    l.append(processProdyFile('{}/{}_effectiveness.txt'.format(prodyDir, uniprot), pos))
    l.append(processProdyFile('{}/{}_sensitivity.txt'.format(prodyDir, uniprot), pos))
    return l

def pyrosettaScore(index, key, chainDict, pdbDir, pyrosettaDir):
    try:
        init()
        uniprot = key.split(',')[0]
        pos = key.split(',')[1]
        var = key.split(',')[3]

        filename = '{}/{}_clean_single.pdb'.format(pdbDir, uniprot)
        pose = pose_from_pdb(filename)
        scorefxn = get_fa_scorefxn()
        scorefxn.set_weight(fa_atr,0)
        scorefxn.set_weight(fa_rep,0)
        scorefxn.set_weight(fa_intra_rep,0)
        scorefxn.set_weight(fa_sol,0)
        scorefxn.set_weight(lk_ball_wtd,0)
        scorefxn.set_weight(fa_intra_sol,0)
        scorefxn.set_weight(fa_elec,0)
        scorefxn.set_weight(hbond_lr_bb,0)
        scorefxn.set_weight(hbond_sr_bb,0)
        scorefxn.set_weight(hbond_bb_sc,0)
        scorefxn.set_weight(hbond_sc,0)
        scorefxn.set_weight(dslf_fa13,0)
        scorefxn.set_weight(rama_prepro,0)
        scorefxn.set_weight(p_aa_pp,0)
        scorefxn.set_weight(fa_dun,0)
        scorefxn.set_weight(omega,0)
        scorefxn.set_weight(pro_close,0)
        scorefxn.set_weight(yhh_planarity,0)
        scorefxn.set_weight(fa_intra_sol_xover4,0)
        wt_score=scorefxn(pose)
        #print (scorefxn.show(pose))
        mutate_residue(pose,pose.pdb_info().pdb2pose(chainDict[uniprot],int(pos)),var)
        mutant_score=scorefxn(pose)
        diff = mutant_score - wt_score
        file = open('{}/{}_{}_{}.txt'.format(pyrosettaDir, uniprot, pos, var), 'w')
        file.write('{} {}\n'.format(mutant_score, diff))
        file.close()
        # print(uniprot, ' pyrosetta ', index, mutant_score, diff)
        return (index, mutant_score, diff)
    except:
        return (index, 'nan', 'nan')

def runP2rank(uniprot, p2rankDir, pdbDir):
    global p2rank
    print('# running P2Rank for {}'.format(uniprot))
    p2rank_file = '{}/{}_clean_single.pdb_residues.csv'.format(p2rankDir, uniprot)
    try:
        os.remove(p2rank_file)
    except:
        pass
    os.system('{}/prank predict -f {}/{}_clean_single.pdb &>/dev/null'.format(p2rank, pdbDir, uniprot))
    shutil.move('{}/test_output/predict_{}_clean_single/{}_clean_single.pdb_residues.csv'.format(p2rank, uniprot, uniprot), p2rankDir)

def trainModel():
    data, target = [], []
    f = open('train.csv', 'r')
    f.readline()
    for line in f:
        line = line.strip().split(',')
        line = [float(i) for i in line]
        data.append(line[1:])
        target.append(0) if line[0] < 3 else target.append(1)
    f.close()
    data = np.array(data)
    target = np.array(target)

    exported_pipeline = XGBClassifier(learning_rate=0.1, max_depth=9, min_child_weight=1,
                                      n_estimators=100, nthread=1, subsample=0.9500000000000001)
    
    exported_pipeline.fit(data, target)

    return exported_pipeline

def runLYRUS(file):
    data = []
    f = open(file, 'r')
    for line in f:
        line = line.strip().split(',')
        line = [float(i) for i in line]
        data.append(line)
    f.close()
    
    data = np.array(data)

    exported_pipeline = trainModel()
    pred = exported_pipeline.predict(data)
    pred_proba = exported_pipeline.predict_proba(data)
    return pred, pred_proba

parser = argparse.ArgumentParser(description='LYRUS: A Machine Learning Model for Predicting the Pathogenicity of Missense Variants')
parser.add_argument('-i', '--inputFile', type=str, metavar='', required=True, help='input file')
parser.add_argument('-o', '--outputDir', type=str, metavar='', required=True, help='output directory; full path')
parser.add_argument('-f', '--fathmm', type=str, metavar='', required=True, help='fathmm file')

args = parser.parse_args()

if __name__ == '__main__':
    currDir = os.getcwd()
    global results
    results = []
    global maestro
    global p2rank
    dirs = os.listdir(currDir)
    for dir in dirs:
        if 'p2rank_' in dir.lower():
            p2rank = dir
        elif 'maestro_' in dir.lower():
            maestro = dir

    inputFile = args.inputFile
    outputDir = args.outputDir
    try:
        os.mkdir(outputDir)
    except:
        pass

    try:
        os.mkdir('{}/data'.format(currDir))
    except:
        pass
    
    cpu_count = max(1, mp.cpu_count() - 2)
    ################################################################################
    '''
    read inputFile and process input
    '''
    ################################################################################
    uniprotList = []
    savDict = {} #key: uniprot; val: [SAVs]
    savList = []
    savAllDict = {}
    file = open(inputFile)
    for line in file:
        line = line.rstrip()
        if len(line) == 0:
            continue
        l = line.split()
        savList.append([l[0], l[1][1:-1], l[1][0], l[1][-1]])
        savAllDict[','.join([l[0], l[1][1:-1], l[1][0], l[1][-1]])] = []
        if l[0] not in uniprotList:
            uniprotList.append(l[0])
            savDict[l[0]] = [l[1]]
        else:
            savDict[l[0]].append(l[1])
    file.close()

    ################################################################################
    '''
    get gene names from uniprot
    '''
    ################################################################################
    getGeneName(uniprotList, outputDir)
    geneDict1 = {} #key: uniprot; val: geneName
    file = open('{}/gene1.txt'.format(outputDir))
    file.readline()
    for line in file:
        line = line.rstrip()
        l = line.split()
        try:
            geneDict1[l[0]] = l[1]
        except:
            pass
    file.close()

    geneDict = {}
    for uniprot in uniprotList:
        if uniprot in geneDict1.keys():
            geneDict[uniprot] = geneDict1[uniprot]
        elif '-' in uniprot:
            if uniprot.split('-')[0] in geneDict1.keys():
                geneDict[uniprot] = geneDict1[uniprot.split('-')[0]]

    '''
    add additional gene information
    '''
    try:
        file = open('gene.txt')
        for line in file:
            line = line.rstrip()
            line = line.split(',')
            geneDict[line[1]] = line[0]
        file.close()
    except:
        pass

    ################################################################################
    '''
    download and clean all the pdb files
    '''
    ################################################################################
    pdbDir = '{}/PDB'.format(outputDir)
    chainDict = {}
    try:
        os.mkdir(pdbDir)
    except:
        pass
    
    logFile = '{}/errorLog.txt'.format(outputDir)
    for uniprot in uniprotList:
        if os.path.isfile('{}/{}_clean_single.pdb'.format(pdbDir, uniprot)):
            getPDB(uniprot, pdbDir)
            getSingleChain(uniprot, pdbDir, savDict, logFile)
        else:
            getPDB(uniprot, pdbDir)
            getSingleChain(uniprot, pdbDir, savDict, logFile)
        try:
            file = open('{}/{}_clean_single.pdb'.format(pdbDir, uniprot))
            ll = file.readline().strip()
            l = [char for char in ll]
            chain_x = l[21].strip()
            chainDict[uniprot] = chain_x
            file.close()
        except:
            pass

    ################################################################################
    '''
    write foldX input files
    '''
    ################################################################################
    foldXDir = '{}/foldXDir'.format(outputDir)
    try:
        os.mkdir(foldXDir)
    except:
        pass

    for uniprot in savDict.keys():
        try:
            savs = savDict[uniprot]
            file = open('{}/individual_list_{}.txt'.format(foldXDir, uniprot), 'w')
            for sav in savs:
                file.write('{}{}{};\n'.format(sav[0], chainDict[uniprot], sav[1:]))
            file.close()
        except:
            pass

    ################################################################################
    '''
    download all the orthologs files
    '''
    ################################################################################
    refseqDir = '{}/refseqs'.format(outputDir)
    try:
        os.mkdir(refseqDir)
    except:
        pass
    
    refseqDict = {} #key: gene; val: acc
    refseqList = []

    try:
        refseqFile = open('refseqs.txt')
        for line in refseqFile:
            line = line.rstrip()
            l = line.split()
            refseqDict[l[0]] = l[1]
        refseqFile.close()
    except:
        pass

    try:
        refseqFile = open('{}/refseqs1.txt'.format(outputDir))
        for line in refseqFile:
            line = line.rstrip()
            l = line.split()
            refseqDict[l[0]] = l[1]
        refseqFile.close()
    except:
        pass

    for file in os.listdir(refseqDir):
        if '.DS_Store' in file:
            continue
        refseqList.append(file)
    
    Entrez.email = ''
    refseqFile = open('{}/refseqs1.txt'.format(outputDir), 'a')
    for uniprot in uniprotList:
        try:
            if geneDict[uniprot] in refseqDict.keys():
                continue
            homoProtein = getFasta(geneDict[uniprot], savDict[uniprot][0], refseqList, refseqDir)
            if len(homoProtein) == 0:
                file = open(logFile, 'a')
                file.write('refseq for {} not found: {}\n'.format(geneDict[uniprot], savDict[uniprot][0]))
                file.close()
            else:
                refseqDict[geneDict[uniprot]] = homoProtein
                refseqFile.write('{} {}\n'.format(geneDict[uniprot], homoProtein))
        except:
            pass
    refseqFile.close()

    vnDir = '{}/vn'.format(outputDir)
    try:
        os.mkdir(vnDir)
    except:
        pass

    ################################################################################
    '''
    find variation number and EVMutation
    '''
    ################################################################################

    '''
    uncomment below
    '''
    refseqList = []
    for file in os.listdir(refseqDir):
        if '.DS_Store' in file:
            continue
    refseqList.append(file)

    pool = mp.Pool(cpu_count)
    for i, ref in enumerate(sorted(refseqList)):
        pool.apply_async(processVN, args=(i, refseqDir, vnDir, ref), callback=collect_result)
    pool.close()
    pool.join()
    
    results.sort(key=lambda x: x[0])
    results_final = [r for i, r in results]
    pool = mp.Pool(cpu_count)
    for ref, homoIndex in zip(sorted(refseqList), results_final):
        pool.apply_async(processEVMutation, args=(vnDir, ref, homoIndex))
    pool.close()
    pool.join()

    for key in savAllDict.keys():
        uniprot = key.split(',')[0]
        try:
            acc = refseqDict[geneDict[uniprot]]
            copyfile('data/vn_{}.txt'.format(acc), '{}/vn_{}.txt'.format(vnDir, acc))
        except:
            pass

        try:
            acc = refseqDict[geneDict[uniprot]]
            copyfile('data/evmutation_{}_result.txt'.format(acc), '{}/evmutation_{}_result.txt'.format(vnDir, acc))
        except:
            pass

        try:
            acc = refseqDict[geneDict[uniprot]]
            file = open('{}/vn_{}.txt'.format(vnDir, acc))
            for line in file:
                line = line.rstrip()
                l = line.split()
                if l[0] == key.split(',')[1]:
                    savAllDict[key].append(float(l[2]))
                    break
            file.close()
        except:
            savAllDict[key].append('nan')

        if len(savAllDict[key]) == 0:
            savAllDict[key].append('nan')

        try:
            acc = refseqDict[geneDict[uniprot]]
            file = open('{}/evmutation_{}_result.txt'.format(vnDir, acc))
            for line in file:
                line = line.rstrip()
                l = line.split(',')
                if l[3] == key.split(',')[1] and l[5] == key.split(',')[3]:
                    savAllDict[key].append(float(l[8]))
                    break
            file.close()
        except:
            savAllDict[key].append('nan')
        
        if len(savAllDict[key]) == 1:
            savAllDict[key].append('nan')

    ################################################################################
    '''
    use rhapsody to find polyphen2
    '''
    ################################################################################
    test_SAVs = []#['{} {} {} {}'.format(uniprot, pos, wt_1codon, var_1codon)]
    
    for sav in savList:
        test_SAVs.append(' '.join(sav))

    rh = rd.Rhapsody()
    rh.queryPolyPhen2(test_SAVs)

    polyphen2_file = 'pph2-full.txt'
    file = open(polyphen2_file)
    for line in file:
        line = line.rstrip()
        if len(line) == 0:
            continue
        line = line.split('\t')
        if '#' in line[0]:
            continue
        for i in range(len(line)):
            line[i] = line[i].replace(' ', '')
        
        k = ','.join(line[0:4])
        try:
            savAllDict[k].insert(1, float(line[22]))
        except:
            savAllDict[k].insert(1, 'nan')
        try:
            savAllDict[k].insert(2, float(line[23]))
        except:
            savAllDict[k].insert(2, 'nan')
    file.close()
    
    for key in savAllDict.keys():
        if len(savAllDict[key]) == 2:
            savAllDict[key].insert(1, 'nan')
            savAllDict[key].insert(2, 'nan')

    pp2Dir = '{}/pp2'.format(outputDir)
    try:
        os.mkdir(pp2Dir)
    except:
        pass
    
    pp2Files = ['rhapsody-SAVs.txt', 'pph2-completed.txt', 'pph2-full.txt', 'pph2-log.txt', 'pph2-short.txt', 'pph2-snps.txt', 'pph2-started.txt']
    for file in pp2Files:
        try:
            os.remove('{}/{}'.format(pp2Dir, file))
        except:
            pass
        shutil.move(file, pp2Dir)

    ################################################################################
    '''
    process fathmm score
    '''
    ################################################################################
    fathmm_file = args.fathmm
    file = open(fathmm_file)
    for line in file:
        line = line.rstrip()
        line = line.split('\t')
        key = line[2] + ',' + line[3][1:-1] + ',' + line[3][0] + ',' + line[3][-1]
        try:
            savAllDict[key].append(float(line[5]))
        except:
            pass
    file.close()
    for key in savAllDict.keys():
        if len(savAllDict[key]) == 4:
            savAllDict[key].append('nan')

    ################################################################################
    '''
    run FoldX
    '''
    ################################################################################
    pool = mp.Pool(cpu_count)
    for uniprot in uniprotList:
        pool.apply_async(runFoldX, args=(uniprot, foldXDir, pdbDir, currDir))
    pool.close()
    pool.join()
    
    for key in savAllDict.keys():
        uniprot = key.split(',')[0]
        try:
            index = 0
            file = open('{}/individual_list_{}.txt'.format(foldXDir, uniprot))
            for line in file:
                index += 1
                line = line.rstrip()
                l = key.split(',')[2] + chainDict[uniprot] + key.split(',')[1] + key.split(',')[3] + ';'
                if l == line:
                    break
            file.close()

            file = open('{}/Dif_{}_clean_single_Repair.fxout'.format(foldXDir, uniprot))
            total = 0
            count = 0
            for line in file:
                line = line.rstrip()
                if count == 5:
                    break
                if '{}_clean_single_Repair_{}_'.format(uniprot, index) in line:
                    l = line.split()
                    total += float(l[1])
                    count += 1
            total = total / count
            file.close()

            savAllDict[key].append(total)
        except:
            savAllDict[key].append('nan')
        
        if len(savAllDict[key]) == 5:
            savAllDict[key].append('nan')
    
    ################################################################################
    '''
    run freeSASA
    '''
    ################################################################################
    freeSASADir = '{}/freeSASA'.format(outputDir)
    try:
        os.mkdir(freeSASADir)
    except:
        pass

    pool = mp.Pool(cpu_count)
    for uniprot in uniprotList:
        pool.apply_async(runFreeSASA, args=(freeSASADir, pdbDir, uniprot))
    pool.close()
    pool.join()

    for key in savAllDict.keys():
        uniprot = key.split(',')[0]
        try:
            file = open('{}/freesasa_{}.txt'.format(freeSASADir, uniprot))
            for line in file:
                line = line.rstrip()
                l = line.split()
                if l[2] == key.split(',')[1]:
                    savAllDict[key].append(float(l[5]))
                    break
            file.close()
        except:
            savAllDict[key].append('nan')

        if len(savAllDict[key]) == 6:
            savAllDict[key].append('nan')
    ################################################################################
    '''
    run MAESTRO
    '''
    ################################################################################
    maestroDir = '{}/maestro'.format(outputDir)
    try:
        os.mkdir(maestroDir)
    except:
        pass
    
    pool = mp.Pool(cpu_count)
    for sav in savAllDict.keys():
        pool.apply_async(runMaestro, args=(sav, maestroDir, chainDict, pdbDir))
    pool.close()
    pool.join()

    for key in savAllDict.keys():
        sav = key.split(',')
        try:
            file = open('{}/{}_{}{}{}_maestro.txt'.format(maestroDir, sav[0], sav[2], sav[1], sav[3]))
            line = file.readline().rstrip().split()
            score_index = 0
            m_index = 0
            for i in range(len(line)):
                if line[i] == 'score':
                    score_index = i
                elif line[i] == 'mutation':
                    m_index = i

            for line in file:
                line = line.rstrip()
                if len(line) == 0:
                    continue
                l = line.split()
                if sav[1] in l[m_index]:
                    savAllDict[key].append(float(l[score_index]))
                    break
            file.close()
        except:
            savAllDict[key].append('nan')

        if len(savAllDict[key]) == 7:
            savAllDict[key].append('nan')

    ################################################################################
    '''
    run prody
    ANM, Mechanical stiffness, Effectiveness, Sensitivity
    '''
    ################################################################################
    prodyDir = '{}/prody'.format(outputDir)
    try:
        os.mkdir(prodyDir)
    except:
        pass

    pool = mp.Pool(cpu_count)
    for uniprot in savDict:
        pool.apply_async(runPrody, args=(uniprot, pdbDir, prodyDir))
    pool.close()
    pool.join()
    
    for key in savAllDict.keys():
        savAllDict[key].extend(getPrody(key, prodyDir))

    ################################################################################
    '''
    pyrosetta
    '''
    ################################################################################

    pyrosettaDir = '{}/pyrosetta'.format(outputDir)
    try:
        os.mkdir(pyrosettaDir)
    except:
        pass

    pool = mp.Pool(cpu_count)
    for i, key in enumerate(savAllDict.keys()):
        pool.apply_async(pyrosettaScore, args=(i, key, chainDict, pdbDir, pyrosettaDir))
    pool.close()
    pool.join()

    for key in savAllDict.keys():
        uniprot = key.split(',')[0]
        pos = key.split(',')[1]
        var = key.split(',')[3]
        try:
            file = open('{}/{}_{}_{}.txt'.format(pyrosettaDir, uniprot, pos, var))
            l = file.readline().rstrip()
            l = l.split()
            savAllDict[key].extend([float(l[0]), float(l[1])])
            file.close()
        except:
            savAllDict[key].extend(['nan', 'nan'])
    
    ################################################################################
    '''
    P2Rank
    '''
    ################################################################################
    p2rankDir = '{}/p2rank'.format(outputDir)
    try:
        os.mkdir(p2rankDir)
    except:
        pass
    
    pool = mp.Pool(cpu_count)
    for uniprot in savDict.keys():
        pool.apply_async(runP2rank, args=(uniprot, p2rankDir, pdbDir))
    pool.close()
    pool.join()

    for key in savAllDict.keys():
        sav = key.split(',')
        file = open('{}/{}_clean_single.pdb_residues.csv'.format(p2rankDir, sav[0]))
        for line in file:
            line = line.rstrip()
            if len(line) == 0:
                continue
            line = line.replace(' ', '')
            l = line.split(',')
            if sav[1] == l[1]:
                if float(l[5]) < 0.5:
                    savAllDict[key].append(0)
                elif float(l[5]) >= 0.5:
                    savAllDict[key].append(1)
                break
        file.close()

        if len(savAllDict[key]) == 14:
            savAllDict[key].append('nan')

    '''
    write file for all savs and features
    '''
    outputFile = open('{}/{}_LYRUS_input.csv'.format(outputDir, inputFile.split('.')[0]), 'w')
    outputFile.write('uniprot,gene,accession,sav,variation_number,dScore,Score1,evmutation,FATHMM,FoldX,SASA,maestro,ANM,MS,effectiveness,sensitivity,pyrosetta_mutant,pyrosetta_diff,active_site\n')
    for key in savAllDict.keys():
        sav = ''.join([key.split(',')[2], key.split(',')[1], key.split(',')[3]])
        uniprot = key.split(',')[0]
        acc = ''
        gene = ''
        try:
            gene = geneDict[uniprot]
        except:
            pass
        try:
            acc = refseqDict[geneDict[uniprot]]
        except:
            pass

        outputFile.write('{},{},{},{}'.format(uniprot, gene, acc , sav))
        for k in savAllDict[key]:
            outputFile.write(',{}'.format(k))
        outputFile.write('\n')
    outputFile.close()

    
    '''
    imput missing values
    '''
    file = '{}/{}_LYRUS_input.csv'.format(outputDir, inputFile.split('.')[0])
    runImpute(file, outputDir, inputFile)
    
    '''
    run LYRUS
    '''
    file = '{}/{}_LYRUS_imputed.csv'.format(outputDir, inputFile.split('.')[0])
    runLYRUS(file)
    pred, pred_proba = runLYRUS(file)

    '''
    write final output file
    '''
    file = open('{}/{}_LYRUS_imputed.csv'.format(outputDir, inputFile.split('.')[0]))
    outputFile = open('{}/{}_LYRUS_prediction.csv'.format(outputDir, inputFile.split('.')[0]), 'w')
    outputFile.write('uniprot,gene,accession,sav,prediction,prediction_proba,variation_number,dScore,Score1,evmutation,FATHMM,FoldX,SASA,maestro,ANM,MS,effectiveness,sensitivity,pyrosetta_mutant,pyrosetta_diff,active_site\n')
    
    for key, line, p, proba in zip(savAllDict.keys(), file, pred, pred_proba):
        sav = ''.join([key.split(',')[2], key.split(',')[1], key.split(',')[3]])
        uniprot = key.split(',')[0]
        acc = ''
        gene = ''
        try:
            gene = geneDict[uniprot]
        except:
            pass
        try:
            acc = refseqDict[geneDict[uniprot]]
        except:
            pass

        line = line.rstrip()

        outputFile.write('{},{},{},{},{},{},{}\n'.format(uniprot, gene, acc, sav, p, proba[1], line))
    file.close()
    outputFile.close()
