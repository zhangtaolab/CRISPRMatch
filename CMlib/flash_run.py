import os
from subprocess import Popen
from CMlib import subprocesspath
import pandas as pd


def flash_merge(flashbin, infofile, output, threadnumber):


    # flashbin = subprocesspath.subprocesspath(flashbin)

    ##/Users/Forrest/SVN/bwa/bwa mem -O 0 -B 0 -E 0 -k 5 ../DM_404.fa oligo_tmp2.fa
    info = pd.read_table(infofile, index_col="Index")
    documentdir = os.path.dirname(os.path.abspath(infofile))
    flashbin = subprocesspath.subprocesspath(flashbin)
    outputdir = os.path.join(os.path.dirname(os.path.abspath(output)),output+'/')
    #print(outputdir)

    for idx in info.index:
        #    print(info.loc[idx].Note, info.loc[idx].gene_name, info.loc[idx].start, info.loc[idx].end)
        fastq1 = os.path.join(documentdir, info.loc[idx]['Sample'] + '/' + info.loc[idx]['Sample'] + '_1.fastq')
        fastq2 = os.path.join(documentdir, info.loc[idx]['Sample'] + '/' + info.loc[idx]['Sample'] + '_2.fastq')
        name = info.loc[idx]['Sample']
        print("Merging",info.loc[idx].Note)


        outfile = os.path.join(documentdir, info.loc[idx]['Sample'] + '/')
        flashcmd = ' '.join([flashbin, '-o', name, '-t', str(threadnumber), '-d', outputdir, fastq1, fastq2, '2>&1 | tee', outputdir + name + '_flash.log'])
        print(flashcmd)
        runflash = Popen(flashcmd, shell=True)
        runflash.communicate()

        singlefqraw = os.path.join(outputdir,info.loc[idx]['Sample'] + '.extendedFrags.fastq')
        singlefqnow = os.path.join(documentdir, info.loc[idx]['Sample'] + '/' + info.loc[idx]['Sample'] + '.extendedFrags.fastq')

        cpcmd = ' '.join(['cp', singlefqraw,  singlefqnow])
        #print(cpcmd)
        os.system(cpcmd)
    return True
