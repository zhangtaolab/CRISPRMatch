import os
import pysam
from pyfasta import Fasta
import matplotlib
from scipy import stats
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def caldel(samfilename, start, end, genename):
    n = 0

    mutateinfor = dict()

    deletelent = dict()
    samfile = pysam.AlignmentFile(samfilename, 'r')
    for read in samfile.fetch(genename):

        # print(read.cigartuples, read.cigarstring, read.reference_start, read.cigartuples[0][1], read.cigartuples[0][1]+read.reference_start)


        nowsite = read.reference_start
        #     print(read.cigarstring)
        for cigarnow in read.cigartuples:
            # print(cigarnow)
            cigartype = cigarnow[0]
            #         print(cigartype)
            cigarlenght = cigarnow[1]

            cigarend = nowsite + cigarlenght

            if start < nowsite < end:

                if cigartype == 2:

                    if cigarlenght < (end - start):

                        if cigarlenght in deletelent:

                            deletelent[cigarlenght] += 1
                        else:
                            deletelent[cigarlenght] = 1

            for i in range(nowsite, cigarend):

                if i in mutateinfor:

                    if cigartype in mutateinfor[i]:

                        mutateinfor[i][cigartype] += 1

                    else:

                        mutateinfor[i][cigartype] = 1

                else:

                    mutateinfor[i] = dict()

                    mutateinfor[i][cigartype] = 1

            nowsite += cigarlenght

        n += 1

    mutateinforpd = pd.DataFrame.from_dict(mutateinfor, orient='index').fillna(value=0)
    mutateinforpd['sum'] = mutateinforpd.sum(axis=1)
    #    print(mutateinforpd[2],mutateinforpd['sum'])
    mutateinforpd['delrate'] = mutateinforpd[2]/mutateinforpd['sum']
    deletelentpd = pd.DataFrame.from_dict(deletelent, orient='index')

    return (mutateinforpd, deletelentpd)

def barchart(infofile,refname, output, bamdir):
    datainfo = pd.read_table(infofile, index_col="Index")
    fa = Fasta(refname)
    for idx in datainfo.index:
        note = datainfo.loc[idx].Note
        if (re.search("gRNA", datainfo.loc[idx].Note)):
            start = datainfo.loc[idx].start - 10
            end = datainfo.loc[idx].end + 10
        elif (re.search("crRNA", datainfo.loc[idx].Note)):
            start = datainfo.loc[idx].start
            end = datainfo.loc[idx].end + 30
        #print(start, end)
        bamfile = os.path.join(bamdir, note + '.bam')
        pdffile = os.path.join(output, note + '.pdf')
        #print(bamfile)
        genename = datainfo.loc[idx].gene_name
        seq = fa[genename][start - 1:end].upper()
        seqlist = list()

        for nt in seq:
            seqlist.append(nt)

        (mutateinforpd, deletelentpd) = caldel(samfilename=bamfile, start=start, end=end, genename=genename)
        reg = mutateinforpd.loc[start:end]
        #print(reg)
        fig, ax = plt.subplots()
        ax.bar(reg.index, reg.delrate, color='blue')
        ax.set_title(note)
        ax.set_xticks(reg.index)
        ax.set_xticklabels(seqlist)
        #    plt.show()
        plt.savefig(pdffile, dpi=300, format="pdf")
        plt.close(fig)
        print(pdffile, "done!")