import sys
import os
import time
outputDir = sys.argv[1]
bashDir = '{}/bash'.format(outputDir)
try:
    os.mkdir(bashDir)
except:
    pass
node = 1
count = 0
with open(sys.argv[2]) as f:
    for line in f:
        line = line.rstrip().split(',')
        gene = line[0]
        uniprot = line[1]
        acc = line[2]
        dir = '{}/{}'.format(outputDir, gene)
        files = os.listdir(dir)
        # print(files)
        if 'finalParam.csv' in files:
            continue
        elif 'lyrus2.txt' in files:
            continue
        else:
            file = open('{}/lyrus2.txt'.format(dir), 'w')
            file.close()
            print(gene)
            if count == 60:
                count = 0
                time.sleep(96400)
            count += 1
            bashFile = open('{}/{}.sh'.format(bashDir, gene), 'w')
            if node < 3:
                node += 1
                bashFile.write('#!/bin/bash\n#SBATCH -n 1\n#SBATCH --mem=120g\n#SBATCH -t 96:00:00\n#SBATCH --account=ccmb-condo\n')
            else:
                node = 1
                bashFile.write('#!/bin/bash\n#SBATCH -n 1\n#SBATCH --mem=120g\n#SBATCH -t 96:00:00\n')
            bashFile.write('#SBATCH -o {}/{}.out'.format(bashDir, gene))
            bashFile.write('\nmodule load python/3.7.4\nmodule load paup/4.0a168\nmodule load clustal_omega/1.2.4\nsource ~/my_cool_science/bin/activate\n')

            bashFile.write('\npython lyrus2.py {} {} {} {}\n'.format(gene, uniprot, acc, outputDir))
            bashFile.write('\ndeactivate')
            bashFile.close()
            os.system('sbatch {}/{}.sh'.format(bashDir, gene))
