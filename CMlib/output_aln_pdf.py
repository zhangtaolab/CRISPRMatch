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

def alnpdftest(infofile, output, refname, groupinfo):
    """

    :param infofile: a description file of details of each sample, example: sample_infor.txt
    :param output: a description file of details of each group, example: group_infor.txt
    :return:
    """
    info = pd.read_table(infofile, index_col="Index")
    fa = Fasta(refname)
    groupinfor = pd.read_table(groupinfo)
    groupinfor.ix[:, pd.isnull(groupinfor).all()] = "UNKNOWN"
    groupinfor = groupinfor.fillna("UNKNOWN")  ##填充表格中NaN处
    stranddict = dict()
    for idy in groupinfor.index:
        stranddict[groupinfor.loc[idy].rep1] = groupinfor.loc[idy].strand
        stranddict[groupinfor.loc[idy].rep2] = groupinfor.loc[idy].strand
        stranddict[groupinfor.loc[idy].control] = groupinfor.loc[idy].strand


    for idx in info.index:

        #bamname =  os.path.join(bamdir, info.loc[idx].Note+'.bam')
        #print("Calculating",bamname)

        note = info.loc[idx].Note
        genename = info.loc[idx]['gene_name']
        strand = stranddict[note]


        if (re.search("gRNA", info.loc[idx].Note)):
            if strand == '+':
                start = info.loc[idx]['start'] - 10
                end = info.loc[idx]['end'] + 10

            else:
                start = info.loc[idx]['start'] - 10
                end = info.loc[idx]['end'] + 10

        elif (re.search("crRNA", info.loc[idx].Note)):
            if strand == '+':
                start = info.loc[idx]['start']
                end = info.loc[idx]['end'] + 30

            else:
                start = info.loc[idx]['start'] - 30
                end = info.loc[idx]['end']



        alnfile = os.path.join(output, info.loc[idx].Note + '_aln.txt')
        outfile = os.path.join(output, info.loc[idx].Note + '_aln.test.pdf')
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
        fig.set_size_inches(0.01 * withset, 0.01 * heightset)
        #fig.set_size_inches(12, 18)
        ax.set_title(info.loc[idx].Note, size=2,fontdict={'family': 'sans-serif'})
        ax.set_ylim(0, heightset)
        ax.set_xlim(0, withset)
        ax.set_yticks([])  ##去掉刻度线
        ax.set_xticks([])
        ax.spines['left'].set_visible(False)  ##设置边框可见性 ax.spines['left'].set_linewidth(0)可设置边框粗细
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ypos = 5

        seq = fa[genename][start - 1:end].upper()  ##reference sequence
        seqlist = list()

        for nt in seq:
            seqlist.append(nt)
        for x in data.index:  # 逐行读取txt 序列
            #        print(x,len(data.loc[x]))
            n = 1
            xpos = 5
            ax.text(4, ypos + 1, data.loc[x][0], size=1, horizontalalignment='right', verticalalignment='center', )

            while n < len(data.loc[x]):
                if (data.loc[x][n] == seq[n-1]):
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
                else:
                    color = "white"

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