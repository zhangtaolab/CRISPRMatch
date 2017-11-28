import os
import pysam
from pyfasta import Fasta
import matplotlib
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
from glob import glob

def alnfile(infofile, refname, output, bamdir):
    #fa = Fasta(refname)
    info = pd.read_table(infofile, index_col="Index")
    for idx in info.index:
        bamname = os.path.join(bamdir, info.loc[idx].Note + '.bam')
        outfile = os.path.join(output, info.loc[idx].Note + '_aln.fa')
        alnfile = os.path.join(output, info.loc[idx].Note + '_aln.txt')
        print("output", info.loc[idx].Note)
        outfa = open(outfile, 'w')
        outlan = open(alnfile, 'w')
        if (re.search("gRNA", info.loc[idx].Note)):
            start = info.loc[idx].start - 10
            end = info.loc[idx].end + 10
        elif (re.search("crRNA", info.loc[idx].Note)):
            start = info.loc[idx].start
            end = info.loc[idx].end + 30
            #start = info.loc[idx].start - 10
            #end = info.loc[idx].end - 10
        gene = info.loc[idx].gene_name
        samfile = pysam.AlignmentFile(bamname, "rb")
        mtreads = set()
        totalcov = 0
        covage = 0

        replace = set()

        insert = set()

        deletion = set()

        reads = dict()

        for pileupcolumn in samfile.pileup(gene, max_depth=50000):

            # print (pileupcolumn.pos, pileupcolumn.n)



            totalcov += pileupcolumn.n
            #        print(pileupcolumn.pos, pileupcolumn.n)

            if end >= pileupcolumn.pos >= start:

                for pileupread in pileupcolumn.pileups:
                    #                print(pileupcolumn.pos, pileupcolumn.n)

                    if pileupread.alignment.query_name not in reads:
                        #                    print(pileupread.alignment.query_name)
                        reads[pileupread.alignment.query_name] = ''

                    if not pileupread.is_del and not pileupread.is_refskip:
                        reads[pileupread.alignment.query_name] += pileupread.alignment.query_sequence[
                            pileupread.query_position]
                        #                    print(reads[pileupread.alignment.query_name])

                        #                    print(pileupread.query_position)
                        #                     querybase = pileupread.alignment.query_sequence[pileupread.query_position]

                        #         #             refbase = pileupread.alignment.get_reference_sequence()[pileupread.query_position]
                        #                     refbase = fa[gene][pileupcolumn.pos].upper()
                        #                     if querybase !=refbase :
                        # #                         replace += 1
                        #                         mtreads.add(pileupread.alignment.query_name)
                        #                         replace.add(pileupread.alignment.query_name)

                        #                 if pileupread.indel > 0:

                        # #                     insert += 1
                        #                     mtreads.add(pileupread.alignment.query_name)
                        #                     insert.add(pileupread.alignment.query_name)
                        #                     print()

                    if pileupread.indel < 0:
                        reads[pileupread.alignment.query_name] += '-' * abs(pileupread.indel)
                        #                    print(reads[pileupread.alignment.query_name])
                        #                    print(reads)
                        # #                     deletion += 1
                        #                     mtreads.add(pileupread.alignment.query_name)
                        #                     deletion.add(pileupread.alignment.query_name)
        lt = end - start + 1
        #    print(lt)
        typdict = dict()
        for i in reads:
            if len(reads[i]) == lt:
                #            print(reads[i])
                if reads[i] in typdict:
                    typdict[reads[i]] += 1
                else:
                    typdict[reads[i]] = 1
        for mutype in typdict:
            print('>', typdict[mutype], sep='', file=outfa)
            print(mutype, file=outfa)
            print(typdict[mutype], '\t'.join(mutype), sep='\t', file=outlan)
        outfa.close()
        outlan.close()

def alnpdf(infofile, output):
    info = pd.read_table(infofile, index_col="Index")
    for idx in info.index:

        alnfile = os.path.join(output, info.loc[idx].Note + '_aln.txt')
        outfile = os.path.join(output, info.loc[idx].Note + '_aln.pdf')
        #print("start output", alnfile, "figure")
        if os.path.getsize(alnfile):  ## check aln file
            print("start output", alnfile, "figure")
        else:
            print("error", alnfile, "figure")
            continue
        data = pd.read_table(alnfile, header=None)  # nrows=400,只读前400行，usecols=(0,1,2,5,6)只提取0，1，2，5，6列
        #    print(len(data.columns)) ##统计列数
        #    print(len(data.index)) ##统计行数
        withset = len(data.columns) * 2 + 10
        heightset = len(data.index) * 2 + 10
        #    print(withset)
        #    print(heightset)
        fig, ax = plt.subplots()
        fig.set_size_inches(12, 18)
        ax.set_title(info.loc[idx].Note)
        ax.set_ylim(0, heightset)
        ax.set_xlim(0, withset)
        ax.set_yticks([])  ##去掉刻度线
        ax.set_xticks([])
        ax.spines['left'].set_visible(False)  ##设置边框可见性 ax.spines['left'].set_linewidth(0)可设置边框粗细
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ypos = 5
        for x in data.index:  # 逐行读取txt 序列
            #        print(x,len(data.loc[x]))
            n = 1
            xpos = 5
            ax.text(4, ypos + 1, data.loc[x][0], size=1, horizontalalignment='right', verticalalignment='center', )

            while n < len(data.loc[x]):
                if (data.loc[x][n] == "A"):
                    color = "red"
                elif (data.loc[x][n] == "T"):
                    color = "blue"
                elif (data.loc[x][n] == "G"):
                    color = "green"
                elif (data.loc[x][n] == "C"):
                    color = "orange"
                else:
                    color = "white"
                # print("n=",n, "data=",data.loc[x][n], "color=", color)

                ax.broken_barh([(xpos, 2)], (ypos, 2), facecolors=color, alpha=0.2)
                ax.text(xpos + 1, ypos + 1, data.loc[x][n], size=1, horizontalalignment='center',
                        verticalalignment='center')
                n += 1
                xpos += 2
            ypos += 2
            # plt.show()
        plt.savefig(outfile, dpi=300, format="pdf")
        plt.close(fig)
        print(outfile, "have finished")