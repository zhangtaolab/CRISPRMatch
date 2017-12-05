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
    """

    :param samfilename: path and name of each bam file
    :param start: setting the start site of deletion calculation
    :param end: setting the end site of deletion calculation
    :param genename: the name of genome-editing target region
    :return:
    """
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

def barchart(infofile,groupinfo,refname, output, bamdir):
    """

    :param infofile: a description file of details of each sample, example: sample_infor.txt
    :param groupinfo: a description file of details of each group, example: group_infor.txt
    :param refname:  a fasta format of the sequence in the target region, exaple:Samples_gene.fa
    :param output: folder of final result
    :param bamdir: folder of temporary files
    :return:
    """
    datainfo = pd.read_table(infofile, index_col="Index")
    groupinfor = pd.read_table(groupinfo)
    stranddict = dict()
    for idy in groupinfor.index:
        stranddict[groupinfor.loc[idy].rep1] = groupinfor.loc[idy].strand
        stranddict[groupinfor.loc[idy].rep2] = groupinfor.loc[idy].strand
        stranddict[groupinfor.loc[idy].control] = groupinfor.loc[idy].strand
    fa = Fasta(refname)
    for idx in datainfo.index:
        note = datainfo.loc[idx].Note
        strand = stranddict[note]
        type = ''
        if (re.search("gRNA", datainfo.loc[idx].Note)):
            if strand == '+':
                start = datainfo.loc[idx]['start'] - 10
                end = datainfo.loc[idx]['end'] + 10
                type = "gf"
            else:
                start = datainfo.loc[idx]['start'] - 10
                end = datainfo.loc[idx]['end'] + 10
                type = "gr"
        elif (re.search("crRNA", datainfo.loc[idx].Note)):
            if strand == '+':
                start = datainfo.loc[idx]['start']
                end = datainfo.loc[idx]['end'] + 30
                type = "cf"
            else:
                start = datainfo.loc[idx]['start'] - 30
                end = datainfo.loc[idx]['end']
                type = "cr"
        # if (re.search("gRNA", datainfo.loc[idx].Note)):
        #     start = datainfo.loc[idx].start - 10
        #     end = datainfo.loc[idx].end + 10
        # elif (re.search("crRNA", datainfo.loc[idx].Note)):
        #     start = datainfo.loc[idx].start
        #     end = datainfo.loc[idx].end + 30
        #print(start, end)
        bamfile = os.path.join(bamdir, note + '.bam')
        pdffile = os.path.join(output, note + '.pdf')
        #print(bamfile)
        genename = datainfo.loc[idx].gene_name
        seq = fa[genename][start - 1:end].upper()
        seqlist = list()
        seqlistPAM = list()
        seqlistother = list()

        for nt in seq:
            seqlist.append(nt)

        (mutateinforpd, deletelentpd) = caldel(samfilename=bamfile, start=start, end=end, genename=genename)
        reg = mutateinforpd.loc[start:end]
        regPAM = list()
        regother = list()
        if type == 'gf':
            seqlistPAM = seqlist[-13:-10]
            seqlistother = seqlist
            seqlistother[-13:-10] = ['', '', '']
            regPAM = reg[-13:-10]
        if type == 'gr':
            seqlistPAM = seqlist[10:13]
            seqlistother = seqlist
            seqlistother[10:13] = ['', '', '']
            regPAM = reg[10:13]
        if type == 'cf':
            seqlistPAM = seqlist[0:4]
            seqlistother = seqlist
            seqlistother[0:4] = ['', '', '', '']
            regPAM = reg[0:4]
        if type == 'cr':
            seqlistPAM = seqlist[-4:]
            seqlistother = seqlist
            seqlistother[-4:] = ['', '', '', '']
            regPAM = reg[-4:]
        #print(reg)
        fig, ax = plt.subplots()
        ax.bar(reg.index, reg.delrate, color='blue')
        ax.set_title(note)
        ax.set_xticks(reg.index, minor=True)
        ax.set_xticklabels(seqlistother, color="black", minor=True)  # minor=True表示次坐标轴
        ax.set_xticks(regPAM.index)
        ax.set_xticklabels(seqlistPAM, color="red")
        # ax.set_xticks(reg.index)
        # ax.set_xticklabels(seqlist)
        #    plt.show()
        plt.savefig(pdffile, dpi=300, format="pdf")
        plt.close(fig)
        print(pdffile, "done!")