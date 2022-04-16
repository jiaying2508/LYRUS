from numpy import outer
import variation_number as vn
import os
import requests
import math
from evcouplings.couplings import CouplingsModel
from evcouplings.mutate import predict_mutation_table, single_mutant_matrix
import pandas as pd
import rhapsody as rd
from shutil import copyfile
import freesasa
import prody
import shutil
import numpy as np
from xgboost import XGBClassifier
import sys

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

class lyrusClass():
    def __init__(self, gene, uniprot, workingDir, email='', savFile=None):
        self.gene = gene
        self.uniprot = uniprot
        self.workingDir = workingDir
        self.seqType = 'protein'
        self.email = email
        outputDir = '{}/{}'.format(workingDir, gene)
        try:
            os.mkdir(outputDir)
        except:
            pass
        self.outputDir = outputDir
        self.acc = None
        try:
            with open('{}/acc.txt'.format(self.outputDir)) as f:
                self.acc = f.readline()
                # print(self.acc)
        except:
            pass
        self.savFile = savFile

    def getFasta(self, acc=None, skip=False):
        """
        Download orthologs sequences from NCBI Orthologs Database

        Parameters
        ----------
        acc(optional): NCBI accession for the gene
        skip(optional): skip the downloading process if skip == True
        """
        if skip:
            if acc is None:
                raise TypeError(
                    'the accession number is None. Either provide the accession number or set skip to be False.')
            self.acc = acc
        else:
            try:
                acc = vn.getFasta(self.gene, self.outputDir,
                                  self.seqType, refseqID=acc, email=self.email)
                self.acc = acc
            except:
                raise TypeError(
                    'the accession number is None. Either provide the accession number or set skip to be False.')

        if acc is None:
                raise TypeError('The accession number is None. Tyr again.')
        file = open('{}/acc.txt'.format(self.outputDir), 'w')
        file.write(acc)
        file.close()
        return acc

    def _readFasta(self, file):
        # read in fasta file and only keep positions if human has a letter
        # dict: sequence dict
        dict = {}
        seq = ''
        acc = ''
        f = open(file)
        for line in f:
            line = line.rstrip()
            if len(line) == 0:
                continue
            if '>' in line:
                if len(acc) > 0:
                    dict[acc] = seq
                seq = ''
                acc = line[1:]
            else:
                seq = seq + line
        dict[acc] = seq
        f.close()

        return dict

    def _generateAAV(self):
        """generate amino acid variation file
        """
        codon3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        seqDict = self._readFasta('{}/{}'.format(self.outputDir, self.gene))
        homoSeq = seqDict[self.acc]
        self.homoSeq = homoSeq
        aaList = codon3to1.values()
        mutationList = []

        for pos in self.savDict:
            wt = self.savDict[pos]
            if homoSeq[pos-1] != wt:
                raise ValueError(
                    'refSeq sequence not match PDB sequence, please provide an alternative uniprot ID or accession number when using getFasta method')

        for i in self.savDict:
            # j = str(i+1)
            for aa in aaList:
                if aa != self.savDict[i]:
                    mutationList.append('{}{}{}'.format(
                        self.savDict[i], str(i), aa))
        file = open('{}/input.txt'.format(self.outputDir), 'w')
        for mutation in mutationList:
            file.write('{}\n'.format(mutation))
        file.close()
        self.savList = mutationList
        
    def _setSAV(self):
        """set custom set of SAVs  and fathmm/polyphen2 input file
        if not provided, use the input.txt file generated using generateAAV
        """
        if self.savFile is None:
            self._generateAAV()
        else:
            self.savList = []
            with open(self.savFile) as f:
                for line in f:
                    if len(line) == 0:
                        continue
                    line = line.rstrip()
                    self.savList.append(line)
        mutationList = ','.join(self.savList)
        file = open('{}/fathmmInput.txt'.format(self.outputDir), 'w')
        file.write('{} {}'.format(self.uniprot, mutationList))
        file.close()

    def setAccession(self, acc):
        """set accession number"""
        self.acc = acc

    def _getPDB(self):
        """download PDB file from SwissModel"""
        upfilename = self.uniprot + '.pdb'
        url = 'https://swissmodel.expasy.org/repository/uniprot/' + \
            upfilename + '?provider=swissmodel'
        r = requests.get(url, allow_redirects=True)
        open('{}/{}'.format(self.outputDir, upfilename), 'wb').write(r.content)
        uppdb_unclean = open('{}/{}'.format(self.outputDir, upfilename))
        uppdb_clean = open(
            '{}/{}_clean.pdb'.format(self.outputDir, self.uniprot), 'w')
        for atom in uppdb_unclean:
            atom = atom.rstrip()
            if atom[:4] == 'ATOM':
                uppdb_clean.write(atom)
                uppdb_clean.write('\n')
        uppdb_unclean.close()
        uppdb_clean.close()

    def getPDB(self):
        """get single chain PDB file"""
        codon3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        self._getPDB()
        chain = 'A'
        oldfile = open(
            '{}/{}_clean.pdb'.format(self.outputDir, self.uniprot), 'r')
        newfile = open(
            '{}/{}_clean_single.pdb'.format(self.outputDir, self.uniprot), 'w')
        lineCount = 0
        savDict = {}
        for line in oldfile:
            if len(line) == 0:
                continue
            line = line.strip()
            l = [char for char in line]
            chain_x = l[21].strip()
            if lineCount == 0:
                chain = chain_x
                lineCount = 1

            if chain_x == chain:
                newfile.write(line)
                newfile.write('\n')

            pos = int(line[22:26].strip())
            wt = codon3to1[line[17:20]]
            savDict[pos] = wt
        self.chain = chain
        self.savDict = savDict
        oldfile.close()
        newfile.close()

    def _foldXInput(self):
       file = open('{}/individual_list.txt'.format(self.outputDir), 'w')
       for sav in self.savList:
           file.write('{}{}{};\n'.format(sav[0], self.chain, sav[1:]))
       file.close()    

    def _write_fasta(self, file, accession, sequence):
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
    
    def _EVMutationInput(self):
        file = open('{}/EVMutationInput.txt'.format(self.outputDir), 'w')
        seqDict = self._readFasta('{}/{}_aligned.fasta'.format(self.outputDir, self.acc))
        homoIndexList = []
        seq = ''
        for i in range(len(seqDict[self.acc])):
            if seqDict[self.acc][i] != '-':
                seq = seq + seqDict[self.acc][i]
                homoIndexList.append(i)
        self._write_fasta(file, self.acc, seq)

        for key in seqDict.keys():
            if key != self.acc:
                seq = ''
                for i in homoIndexList:
                    seq = seq + seqDict[key][i]
                self._write_fasta(file, key, seq)
        file.close()
        return homoIndexList
        
    def _processEVMutation(self, homoIndexList, param=100):
        ################################################################################
        #                       calculate evmutation
        ################################################################################
        ll = (len(homoIndexList) - 1) * 0.2
        print('# performing EVMutation Analysis for {}'.format(self.gene))

        os.system('plmc-master/bin/plmc -o {}/{}.params -le {} -lh 0.01 -m {} {}/EVMutationInput.txt &>/dev/null'.format(self.outputDir, self.acc, ll, param, self.outputDir))
        c = CouplingsModel('{}/{}.params'.format(self.outputDir, self.acc))
        singles = single_mutant_matrix(c, output_column='effect_prediction_epistatic')

        pd.DataFrame(singles).to_csv('{}/evmutation_result.txt'.format(self.outputDir))
        os.remove('{}/{}.params'.format(self.outputDir, self.acc))
  
    def _runFoldX(self):
        currDir = os.getcwd()
        pdbfile = '{}/{}_clean_single.pdb'.format(self.outputDir, self.uniprot)
        print('# running foldX for {}'.format(self.uniprot))
        copyfile(pdbfile, '{}/{}_clean_single.pdb'.format(currDir, self.uniprot))
        repairPDBfile = '{}_clean_single_Repair.pdb'.format(self.uniprot)
        # print('./foldx --command=RepairPDB --pdb={}_clean_single.pdb --output-dir={} &>/dev/null'.format(self.uniprot,self.outputDir))
        os.system('./foldx --command=RepairPDB --pdb={}_clean_single.pdb --output-dir={} &>/dev/null'.format(self.uniprot,self.outputDir))
        os.remove('{}/{}_clean_single.pdb'.format(currDir, self.uniprot))

        foldx_out = '{}/Dif_{}_clean_single_Repair.fxout'.format(self.outputDir,self.uniprot)
        try:
            os.remove(foldx_out)
        except:
            pass
        
        copyfile('{}/{}'.format(self.outputDir, repairPDBfile), '{}/{}'.format(currDir, repairPDBfile))
        # print('./foldx --command=BuildModel --pdb={} --mutant-file={}/individual_list.txt --output-dir={} --numberOfRuns=5 --out-pdb=false &>/dev/null'.format(repairPDBfile,self.outputDir,self.outputDir))
        os.system('./foldx --command=BuildModel --pdb={} --mutant-file={}/individual_list.txt --output-dir={} --numberOfRuns=3 --out-pdb=false &>/dev/null'.format(repairPDBfile,self.outputDir,self.outputDir))
        os.remove('{}/{}'.format(currDir, repairPDBfile))
    
    def _runFreesasa(self):
        print('running freesasa')
        structure = freesasa.Structure('{}/{}_clean_single.pdb'.format(self.outputDir, self.uniprot))
        result = freesasa.calc(structure,
        freesasa.Parameters({'algorithm' : freesasa.ShrakeRupley,
        'n-slices' : 200, 'probe-radius' : 1.2}))
        file = open('{}/freesasa.txt'.format(self.outputDir), 'w')
        for key in result.residueAreas()[self.chain]:
            # print(key, result.residueAreas()[self.chain][key].total)
            file.write('{} {:0.2f}\n'.format(key, result.residueAreas()[self.chain][key].total))
        file.close()
    
    def _runMaestro(self, maestro):
        print('# running MAESTRO for {}'.format(self.gene))
        pdbfile = '{}/{}_clean_single.pdb'.format(self.outputDir, self.uniprot)
        file = open('{}/maestro.txt'.format(self.outputDir), 'w')
        for sav in self.savList:
            wt = sav[0]
            var = sav[-1]
            pos = sav[1:-1]
            mut = '{}{}.{}{{{}}}'.format(wt,pos,self.chain,var)
            maestro_file = '{}/{}_maestro.txt'.format(self.outputDir, sav)
            os.system('{}/maestro {}/config.xml {} --evalmut=\'{}\' --bu > {}'.format(maestro, maestro, pdbfile, mut, maestro_file))
            with open('{}/{}_maestro.txt'.format(self.outputDir, sav)) as f:
                line = f.readline().rstrip().split()
                score_index = 0
                m_index = 0
                for i in range(len(line)):
                    if line[i] == 'score':
                        score_index = i
                    elif line[i] == 'mutation':
                        m_index = i
                
                for line in f:
                    if len(line) == 0:
                        continue
                    line = line.rstrip()
                    if len(line) == 0:
                        continue
                    l = line.split()
                    if pos in l[m_index]:
                        file.write('{} {}\n'.format(sav, l[score_index]))
                        break
            os.remove('{}/{}_maestro.txt'.format(self.outputDir, sav))
        file.close()

    def _runPrody(self):
        filename = '{}/{}_clean_single.pdb'.format(self.outputDir, self.uniprot)
        print('# running Prody for {}'.format(self.uniprot))
        prot = prody.parsePDB(filename)
        #ANM
        anm, atoms = prody.calcANM(prot, selstr='calpha')
        flucts = prody.calcSqFlucts(anm)

        pdbIndexDict = {}
        file = open(filename)
        index = -1
        currPos = 0
        count = 0
        
        for line in file:
            if len(line) == 0:
                continue
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

        file1 = open('{}/{}_anm.txt'.format(self.outputDir, self.uniprot), 'w')
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

        file1 = open('{}/{}_ms.txt'.format(self.outputDir, self.uniprot), 'w')
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

        file1 = open('{}/{}_effectiveness.txt'.format(self.outputDir, self.uniprot), 'w')
        for pos in pdbIndexDict.keys():
            file1.write('{},{}\n'.format(pos, eff[pdbIndexDict[pos]]))
        file1.close()

        file1 = open('{}/{}_sensitivity.txt'.format(self.outputDir, self.uniprot), 'w')
        for pos in pdbIndexDict.keys():
            file1.write('{},{}\n'.format(pos, sen[pdbIndexDict[pos]]))
        file1.close()

    def _runPyrosetta(self):
        print("running PyRosetta for {}".format(self.gene))
        file = open('{}/pyrosetta.txt'.format(self.outputDir), 'w')
        filename = '{}/{}_clean_single.pdb'.format(self.outputDir, self.uniprot)
        uniprot = self.uniprot
        
        for sav in self.savList:
            init()
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
            wt_score = scorefxn(pose)
            # print (scorefxn.show(pose))
        
            pos = sav[1:-1]
            var = sav[-1]
            mutate_residue(pose,pose.pdb_info().pdb2pose(self.chain,int(pos)),var)
            mutant_score = scorefxn(pose)
            diff = mutant_score - wt_score
            file.write('{} {} {}\n'.format(sav, mutant_score, diff))
        file.close()

    def _runP2rank(self, p2rank):
        print('# running P2Rank for {}'.format(self.uniprot))
        p2rank_file = '{}/{}_clean_single.pdb_residues.csv'.format(self.outputDir, self.uniprot)
        try:
            os.remove(p2rank_file)
        except:
            pass
        os.system('{}/prank predict -f {}/{}_clean_single.pdb &>/dev/null'.format(p2rank, self.outputDir, self.uniprot))
        shutil.move('{}/test_output/predict_{}_clean_single/{}_clean_single.pdb_residues.csv'.format(p2rank, self.uniprot, self.uniprot), self.outputDir)
        shutil.rmtree('{}/test_output/predict_{}_clean_single'.format(p2rank, self.uniprot))

    def _runPolyphen2(self):
        curr = os.getcwd()
        os.chdir(self.outputDir)
        # savList = ' '.join(self.savList)
        test_SAVs = []
        for sav in self.savList:
            s = self.uniprot + ' ' + sav[1:-1] + ' ' + sav[0] + ' ' + sav[-1]
            test_SAVs.append(s)
        rh = rd.Rhapsody()
        rh.queryPolyPhen2(test_SAVs)
        os.chdir(curr)

    def _cleanDir(self):
        # os.remove('{}/taxid.txt'.format(self.outputDir))
        os.remove('{}/{}.pdb'.format(self.outputDir, self.uniprot))
        # os.remove('{}/{}_aligned.fasta'.format(self.outputDir, self.acc))
        os.remove('{}/{}_getTree.cmd'.format(self.outputDir, self.acc))
        os.remove('{}/{}_tree.tree'.format(self.outputDir, self.acc))
        os.remove('{}/{}_trees.nex'.format(self.outputDir, self.acc))
        os.remove('{}/{}.fasta'.format(self.outputDir, self.acc))
        os.remove('{}/{}.nex'.format(self.outputDir, self.acc))
        os.remove('{}/EVMutationInput.txt'.format(self.outputDir))
        os.remove('{}/Average_{}_clean_single_Repair.fxout'.format(self.outputDir, self.uniprot))
        os.remove('{}/PdbList_{}_clean_single_Repair.fxout'.format(self.outputDir, self.uniprot))
        os.remove('{}/Raw_{}_clean_single_Repair.fxout'.format(self.outputDir, self.uniprot))
        os.remove('{}/{}_clean_single_Repair.fxout'.format(self.outputDir, self.uniprot))
        os.remove('{}/pph2-short.txt'.format(self.outputDir))
        os.remove('{}/pph2-log.txt'.format(self.outputDir))
        os.remove('{}/pph2-snps.txt'.format(self.outputDir))
        os.remove('{}/pph2-started.txt'.format(self.outputDir))
        os.remove('{}/rhapsody-SAVs.txt'.format(self.outputDir))

    def _finalParam(self):
        vnDict = {} #key:M1
        with open('{}/vn_{}.txt'.format(self.outputDir, self.acc)) as f:
            for line in f:
                if len(line) == 0:
                    continue
                line = line.rstrip().split()
                key = line[1] + line[0]
                # print(key)
                vnDict[key] = line[2]
                # exit()
        p2Dict = {} #key: M1A
        p21Dict = {}
        with open('{}/pph2-full.txt'.format(self.outputDir)) as f:
            for line in f:
                if line.startswith(self.uniprot):
                    line = line.rstrip().split('\t')
                    l = [x.strip() for x in line]
                    k = l[2] + l[1] + l[3]
                    p2Dict[k] = l[22]
                    p21Dict[k] = l[23]
        evDict = {} #key:1A
        with open('{}/evmutation_result.txt'.format(self.outputDir)) as f:
            f.readline()
            for line in f:
                if len(line) == 0:
                    continue
                line = line.rstrip().split(',')
                key = line[3] + line[5]
                evDict[key] = line[8]
        foldxDict = {} #key:M1A
        foldxIndexDict = {}
        with open('{}/individual_list.txt'.format(self.outputDir)) as f:
            index = 0
            for line in f:
                if len(line) == 0:
                    continue
                line = line.rstrip()
                foldxIndexDict[index] = line[0] + line[2:-1]
                index += 1
        with open('{}/Dif_{}_clean_single_Repair.fxout'.format(self.outputDir, self.uniprot)) as f:
            index = 0
            count = 0
            total = 0
            countT = 3
            for line in f:
                if line.startswith(self.uniprot):
                    if count == countT:
                        count = 0
                        total = total / countT
                        foldxDict[foldxIndexDict[index]] = str(total)
                        total = 0
                        index += 1
                    line = line.rstrip().split()
                    total += float(line[1])
                    count += 1
            total = total / countT
            foldxDict[foldxIndexDict[index]] = str(total)
        freesasaDict = {} #key:1
        with open('{}/freesasa.txt'.format(self.outputDir)) as f:
            for line in f:
                if len(line) == 0:
                    continue
                line = line.rstrip().split()
                freesasaDict[line[0]] = line[1]
        maestroDict = {} #key: M1A
        with open('{}/maestro.txt'.format(self.outputDir)) as f:
            for line in f:
                if len(line) == 0:
                    continue
                line = line.rstrip().split()
                maestroDict[line[0]] = line[1]
        anmDict = {} #key: 1
        with open('{}/{}_anm.txt'.format(self.outputDir, self.uniprot)) as f:
            for line in f:
                if len(line) == 0:
                    continue
                line = line.rstrip().split(',')
                anmDict[line[0]] = line[1]
        msDict = {}
        with open('{}/{}_ms.txt'.format(self.outputDir, self.uniprot)) as f:
            for line in f:
                line = line.rstrip().split(',')
                msDict[line[0]] = line[1]
        effDict = {}
        with open('{}/{}_effectiveness.txt'.format(self.outputDir, self.uniprot)) as f:
            for line in f:
                if len(line) == 0:
                    continue
                line = line.rstrip().split(',')
                effDict[line[0]] = line[1]
        senDict = {}
        with open('{}/{}_sensitivity.txt'.format(self.outputDir, self.uniprot)) as f:
            for line in f:
                if len(line) == 0:
                    continue
                line = line.rstrip().split(',')
                senDict[line[0]] = line[1]
        pymutDict = {} #key: M1A
        pydiffDict = {}
        with open('{}/pyrosetta.txt'.format(self.outputDir)) as f:
            for line in f:
                if len(line) == 0:
                    continue
                line = line.rstrip().split()
                pymutDict[line[0]] = line[1]
                pydiffDict[line[0]] = line[2]
        p2rankDict = {} #key: 1
        with open('{}/{}_clean_single.pdb_residues.csv'.format(self.outputDir, self.uniprot)) as f:
            f.readline()
            for line in f:
                if len(line) == 0:
                    continue
                line = line.rstrip().split(',')
                val = '0'
                if float(line[5]) >= 0.5:
                    val = '1'
                p2rankDict[line[1].strip()] = val

        outputFile = open('{}/finalParam.csv'.format(self.outputDir), 'w')
        outputFile.write('sav,vn,dScore,Score1,evmutation,foldx,freesasa,maestro,anm,ms,effectiveness,sensitivity,prMutant,prDiff,p2rank\n')

        for sav in self.savList:
            try:
                out = [sav]
                out.append(vnDict[sav[:-1]])
                out.append(p2Dict[sav])
                out.append(p21Dict[sav])
                out.append(evDict[sav[1:]])
                out.append(foldxDict[sav])
                out.append(freesasaDict[sav[1:-1]])
                out.append(maestroDict[sav])
                out.append(anmDict[sav[1:-1]])
                out.append(msDict[sav[1:-1]])
                out.append(effDict[sav[1:-1]])
                out.append(senDict[sav[1:-1]])
                out.append(pymutDict[sav])
                out.append(pydiffDict[sav])
                out.append(p2rankDict[sav[1:-1]])
                out = ','.join(out)
                outputFile.write('{}\n'.format(out))
            except:
                pass
        outputFile.close()

    def getParameters(self, maestroDir='MAESTRO_OSX_x64',p2rankDir='p2rank_2.2', EVparam=100):
        """generate input parameters for LYRUS"""
        self._setSAV()

        #variation_number
        vn.processVN(self.gene, self.outputDir, self.acc, 'protein')

        #evmutation
        homoIndexList = self._EVMutationInput()
        self._processEVMutation(homoIndexList, param=EVparam)

        #ployphen2
        self._runPolyphen2()
        
        #p2rank
        self._runP2rank(p2rankDir)

        #pyrosetta
        self._runPyrosetta()

        #prody
        self._runPrody()

        #maestro
        self._runMaestro(maestroDir)

        #freesasa
        self._runFreesasa()

        #foldx
        self._foldXInput()
        self._runFoldX()

        self._finalParam()
        self._cleanDir()

def _trainModel():
    # print(os.path.realpath(__file__))
    # print(os.path.abspath(os.path.dirname(__file__)))
    data, target = [], []
    f = open('{}/train.csv'.format(os.path.abspath(os.path.dirname(__file__))), 'r')
    f.readline()
    for line in f:
        if len(line) == 0:
            continue
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

def lyrusPredict(gene, fathmmFile, outputDir, uniprot):
    
    file = '{}/{}/finalParam.csv'.format(outputDir, gene)
    fathmmDict = {}
    with open(fathmmFile) as f:
        f.readline()
        for line in f:
            if len(line) == 0:
                continue
            # print(line)
            line = line.rstrip()
            line = line.split('\t')
            # print(line)
            if line[2] == uniprot:
                try:
                    fathmmDict[line[3]] = line[5]
                except:
                    pass
    
    outputFile = open('{}/{}/lyrusPredict.csv'.format(outputDir, gene), 'w')
    outputFile.write('sav,lyrusScore,lyrusProbability,vn,dScore,Score1,evmutation,fathmm,foldx,freesasa,maestro,anm,ms,effectiveness,sensitivity,prMutant,prDiff,p2rank\n')
    data = []
    savList = []
    with open(file) as f:
        f.readline()
        for line in f:
            if len(line) == 0:
                continue
            line = line.rstrip().split(',')
            sav = line[0]
            if sav in fathmmDict.keys():
                line.insert(5, fathmmDict[sav])
                l = [float(i) for i in line[1:]]
                # print(l)
                data.append(l)
                savList.append(line)
    data = np.array(data)
    # print(data)
    exported_pipeline = _trainModel()
    pred = exported_pipeline.predict(data)
    pred_proba = exported_pipeline.predict_proba(data)
    for p1, p2, sav in zip(pred, pred_proba, savList):
        sav.insert(1, str(p1))
        sav.insert(2, str(p2[1]))
        s = ','.join(sav)
        outputFile.write('{}\n'.format(s))
    outputFile.close()
    # exit()
