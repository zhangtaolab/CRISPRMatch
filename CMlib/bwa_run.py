import os
from glob import glob
import pandas as pd
from pyfasta import Fasta
from subprocess import Popen

def prepare(infofile, refname, output, bwabin, samtoolsbin, picardbin):
    datainfo=pd.read_table(infofile,index_col="Index")
    outputname = os.path.join(output, 'bwa_run.sh')
    documentdir = os.path.dirname(infofile)
    outio = open(outputname,"w")
    for idx in datainfo.index:
        fqname = documentdir+'/'+datainfo.ix[idx]['Sample']+'/'+datainfo.ix[idx]['Sample']+'.extendedFrags.fastq'
        bamfile = datainfo.ix[idx]['Note'] + '.bam'
        print(bwabin,' mem ', os.path.basename(refname), ' ', fqname, ' | ',picardbin,' SortSam I=/dev/stdin O=', bamfile,
              ' SO=coordinate', sep='', file=outio)
        print(samtoolsbin,' index ',bamfile, file=outio)
        #print('bwa mem ', os.path.basename(refname), ' ', fqname, ' | picard SortSam I=/dev/stdin O=', bamfile, ' SO=coordinate', sep='',
        #      file=outio)
        #print('samtools index ', bamfile, file=outio)
    outio.close()
    print("bwa command load!")

    ###run bwa mem
    bwacmd="bash bwa_run.sh"
    print(bwacmd)
    runbwaalign = Popen(bwacmd, shell=True, cwd=output)
    runbwaalign.communicate()

    print("bwa mem finished")

    return True