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

def rate_cal(infofile, groupinfo, refname, output, bamdir):
    """

    :param infofile: a description file of details of each sample, example: sample_infor.txt
    :param groupinfo: a description file of details of each group, example: group_infor.txt
    :param refname: a fasta format of the sequence in the target region, exaple:Samples_gene.fa
    :param output: folder of final result
    :param bamdir: folder of temporary files
    :return:
    """


    info=pd.read_table(infofile,index_col="Index")
    groupinfor = pd.read_table(groupinfo)
    #print(groupinfor)
    #groupinfor = groupinfor.dropna(axis=0, how='any',thresh=7)  ##过滤表哥中没填满的行，thresh=7表示至少7个数不是NA
    groupinfor.ix[:,pd.isnull(groupinfor).all()] = "UNKNOWN"
    groupinfor=groupinfor.fillna("UNKNOWN") ##填充表格中NaN处
    stranddict = dict()
    for idy in groupinfor.index:
        stranddict[groupinfor.loc[idy].rep1] = groupinfor.loc[idy].strand
        stranddict[groupinfor.loc[idy].rep2] = groupinfor.loc[idy].strand
        stranddict[groupinfor.loc[idy].control] = groupinfor.loc[idy].strand

    outputname = os.path.join(output, 'mut_rate.txt')
    outio = open(outputname, "w")
    fa = Fasta(refname)
    print("start calculation!")
    print("Sample\tmuation\treplace\tinsertion_only\tdeletion_only\tinsert&deletion", file=outio)


    for idx in info.index:
        #    print(info.loc[idx].Note, info.loc[idx].gene_name, info.loc[idx].start, info.loc[idx].end)
        bamname =  os.path.join(bamdir, info.loc[idx].Note+'.bam')
        print("Calculating",bamname)

        note = info.loc[idx].Note
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


        # if (re.search("gRNA", info.loc[idx].Note)):
        #     start = info.loc[idx].start - 10
        #     end = info.loc[idx].end + 10
        # # print(info.loc[idx].Note, "orignal-start",info.loc[idx].start, "after:",start, "orignal-end",info.loc[idx].end, "after:",end)
        # elif (re.search("crRNA", info.loc[idx].Note)):
        #     start = info.loc[idx].start
        #     end = info.loc[idx].end + 30
        # print(info.loc[idx].Note, "orignal-start",info.loc[idx].start, "after:",start, "orignal-end",info.loc[idx].end, "after:",end)
        #    start = info.loc[idx].start
        #    end = info.loc[idx].end
        gene = info.loc[idx].gene_name
        samfile = pysam.AlignmentFile(bamname, "rb", check_sq=False)
        #print(samfile.count())
        mtreads = set()
        totalcov = 0
        #covage = samfile.count()/100
        covage = set()

        replace = set()

        insert = set()

        deletion = set()

        mthreads_all = set()

        insert_deletion = set()

        insert_only = set()

        deletion_only = set()

        for pileupcolumn in samfile.pileup(gene, max_depth=50000):

            # print (pileupcolumn.pos, pileupcolumn.n)  pos代表该位点的坐标，n代表它的coverage,pileups代表对应的reads

            totalcov += pileupcolumn.n

            if end >= pileupcolumn.pos >= start:


                for pileupread in pileupcolumn.pileups:

                    covage.add(pileupread.alignment.query_name)

                    if not pileupread.is_del and not pileupread.is_refskip:
                        #             print(pileupread.query_position)
                        querybase = pileupread.alignment.query_sequence[pileupread.query_position]

                        #             refbase = pileupread.alignment.get_reference_sequence()[pileupread.query_position]

                        refbase = fa[gene][pileupcolumn.pos].upper()

                        if querybase != refbase:
                            #                         replace += 1
                            mtreads.add(pileupread.alignment.query_name)

                            replace.add(pileupread.alignment.query_name)

                            # pileupread.indel: 在当前的pileup位点之后的位置的indel长度。如果下一个位点是insertion,indel>0；如果下一个位点是deletion,indel<0
                    if pileupread.indel > 0:
                        #                     insert += 1
                        mtreads.add(pileupread.alignment.query_name)
                        insert.add(pileupread.alignment.query_name)
                        mthreads_all.add(pileupread.alignment.query_name)

                    if pileupread.indel < 0:
                        #                     deletion += 1
                        mtreads.add(pileupread.alignment.query_name)
                        deletion.add(pileupread.alignment.query_name)
                        mthreads_all.add(pileupread.alignment.query_name)

                        #         print(pileupcolumn.pos, pileupcolumn.n, replace, insert, deletion)
        insert_deletion = insert & deletion
        insert_only = insert - deletion
        deletion_only = deletion - insert
        #mthreads_all = mtreads - replace
        # print(info.loc[idx].Note, end='\t', file=outio)
        # print(len(mtreads)/500, end='\t', file=outio)
        # print(len(replace)/500, end='\t', file=outio)
        # print(len(insert_only)/500, end='\t', file=outio)
        # print(len(deletion_only)/500, end='\t', file=outio)
        # print(len(insert_deletion)/500, end='\n', file=outio)

        #print(len(covage))
        print(info.loc[idx].Note, end='\t', file=outio)
        #print(len(mtreads)/len(covage)*100, end='\t', file=outio)
        print(len(mthreads_all)/len(covage) * 100, end='\t', file=outio)
        print(len(replace)/len(covage)*100, end='\t', file=outio)
        print(len(insert_only)/len(covage)*100, end='\t', file=outio)
        print(len(deletion_only)/len(covage)*100, end='\t', file=outio)
        print(len(insert_deletion)/len(covage)*100, end='\n', file=outio)

        # print(info.loc[idx].Note.'\t'.len(mtreads)/500,len(replace)/500,'\t',len(insert_only)/500,'\t',len(deletion_only)/500,'\t',len(insert_deletion)/500, file=outio)
        # print(info.loc[idx].Note, len(mtreads)/500, len(replace)/500, len(insert)/500, len(deletion)/500 )
        samfile.close()
    outio.close()

def display(groupinfo, output):
    """
    :param groupinfo: a description file of details of each group, example: group_infor.txt
    :param output: folder of final result
    :return:
    """

    mutfile = os.path.join(output, 'mut_rate.txt')
    mut_rate = pd.read_table(mutfile, sep='\t')
    groupinfor = pd.read_table(groupinfo)
    groupinfor = groupinfor.dropna(axis=0, how='any',thresh=6)  ##过滤表哥中没填满的行，thresh=7表示至少7个数不是NA,控制treatment和CK至少有一个
    #groupinfor.ix[:,pd.isnull(groupinfor).all()] = "UNKNOWN"
    groupinfor=groupinfor.fillna("UNKNOWN") ##填充表格中NaN处
    #print(groupinfor)
    mut_result = dict()
    for idx in mut_rate.index:
        mut_result[mut_rate.loc[idx].Sample] = mut_rate.values[idx]  ##读入mutation信息
    #    mut_result['OsPDS-RZ-gRNA1_Rep1'][2]

    ## prepare for display
    #replace = list()
    #replace_yerr = list()
    mutation=list()
    mutation_yerr = list()
    insertO = list()
    insertO_yerr = list()
    deletionO = list()
    deletionO_yerr = list()
    insert_deletion = list()
    insert_deletion_yerr = list()
    glist = list()
    ck_glist = list()

    ck_mutation = list()
    ck_insertO = list()
    ck_deletionO = list()
    ck_insert_deletion = list()
    for idy in groupinfor.index:
        rep1 = groupinfor.loc[idy].rep1
        rep2 = groupinfor.loc[idy].rep2
        ck = groupinfor.loc[idy].control
        if (mut_result.__contains__(rep1) and mut_result.__contains__(rep2)):

            # replace_mean = np.mean([mut_result[rep1][2], mut_result[rep2][2]])  ##np.mean([1,2,3,4,5])
            # #    print(group_mean)
            # replace.append(replace_mean)
            # replace_std = np.std([mut_result[rep1][2], mut_result[rep2][2]])  ## 标准差
            # #    print("std", group_var)
            # replace_yerr.append(replace_std)
            mutation_mean = np.mean([mut_result[rep1][1], mut_result[rep2][1]])  ##np.mean([1,2,3,4,5])
            #    print(group_mean)
            mutation.append(mutation_mean)
            mutation_std = np.std([mut_result[rep1][1], mut_result[rep2][1]])  ## 标准差
            #    print("std", group_var)
            mutation_yerr.append(mutation_std)

            insertO_mean = np.mean([mut_result[rep1][3], mut_result[rep2][3]])
            insertO.append(insertO_mean)
            insertO_std = np.std([mut_result[rep1][3], mut_result[rep2][3]])
            insertO_yerr.append(insertO_std)

            deletionO_mean = np.mean([mut_result[rep1][4], mut_result[rep2][4]])
            deletionO.append(deletionO_mean)
            deletionO_std = np.std([mut_result[rep1][4], mut_result[rep2][4]])
            deletionO_yerr.append(deletionO_std)

            insert_deletion_mean = np.mean([mut_result[rep1][5], mut_result[rep2][5]])
            insert_deletion.append(insert_deletion_mean)
            insert_deletion_std = np.std([mut_result[rep1][5], mut_result[rep2][5]])
            insert_deletion_yerr.append(insert_deletion_std)
        elif mut_result.__contains__(rep1):
            print("The group:",groupinfor.loc[idy].group, ": Rep2 is missing.")
            mutation.append(mut_result[rep1][1])
            mutation_yerr.append(0)
            insertO.append(mut_result[rep1][3])
            insertO_yerr.append(0)
            deletionO.append(mut_result[rep1][4])
            deletionO_yerr.append(0)
            insert_deletion.append(mut_result[rep1][5])
            insert_deletion_yerr.append(0)
        elif mut_result.__contains__(rep2):
            print("The group:", groupinfor.loc[idy].group, ": Rep1 is missing.")
            mutation.append(mut_result[rep2][1])
            mutation_yerr.append(0)
            insertO.append(mut_result[rep2][3])
            insertO_yerr.append(0)
            deletionO.append(mut_result[rep2][4])
            deletionO_yerr.append(0)
            insert_deletion.append(mut_result[rep2][5])
            insert_deletion_yerr.append(0)
        else:
            print("All repetitions in group:", groupinfor.loc[idy].group, " is missing.")
            mutation.append(0)
            mutation_yerr.append(0)
            insertO.append(0)
            insertO_yerr.append(0)
            deletionO.append(0)
            deletionO_yerr.append(0)
            insert_deletion.append(0)
            insert_deletion_yerr.append(0)

        if ck=='UNKNOWN':
            print("The group:",groupinfor.loc[idy].group, ": CK is missing.")
            ck_mutation.append(0)
            ck_insertO.append(0)
            ck_deletionO.append(0)
            ck_insert_deletion.append(0)
        else:
            ck_mutation.append(mut_result[ck][1])
            ck_insertO.append(mut_result[ck][3])
            ck_deletionO.append(mut_result[ck][4])
            ck_insert_deletion.append(mut_result[ck][5])


        glist.append(groupinfor.loc[idy].group)
        ck_glist.append(groupinfor.loc[idy].control)
    ## prepare for display

    ## print out pdf
    mutfile = os.path.join(output, 'mut_result.pdf')
    fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True)
    fig.set_size_inches(20, 9)
    width = 0.15
    bar1 = ax0.bar(groupinfor.index, mutation, width, color="#CC79A7")
    #bar1 = ax0.bar(groupinfor.index, replace, width, color='pink', yerr=replace_yerr, elinewidth=0.1, capsize=1.5)
    ax0.errorbar(groupinfor.index, mutation, yerr=mutation_yerr, fmt='', elinewidth=0.5, capsize=2, capthick=0.5, ls='None', ecolor='black')
    bar2 = ax0.bar(groupinfor.index + width, deletionO, width, color="#D55E00")
    #bar2 = ax0.bar(groupinfor.index + width, insertO, width, color='green', yerr=insertO_yerr, linewidth=0.5,capsize=1.5)
    ax0.errorbar(groupinfor.index+ width, deletionO, yerr=deletionO_yerr, fmt='', elinewidth=0.5, capsize=2, capthick=0.5, ls='None', ecolor='black')
    #bar3 = ax0.bar(groupinfor.index + width * 2, deletionO, width, color='blue', yerr=deletionO_yerr, linewidth=0.5,capsize=1.5)
    bar3 = ax0.bar(groupinfor.index + width * 2, insertO, width, color="#0072B2")
    ax0.errorbar(groupinfor.index + width * 2, insertO, yerr=insertO_yerr, fmt='', elinewidth=0.5, capsize=2, capthick=0.5,ls='None', ecolor='black')
    bar4 = ax0.bar(groupinfor.index + width * 3, insert_deletion, width, color="#009E73")
    #bar4 = ax0.bar(groupinfor.index + width * 3, insert_deletion, width, color='orange', yerr=insert_deletion_yerr,linewidth=0.5, capsize=1.5)
    ax0.errorbar(groupinfor.index + width * 3, insert_deletion, yerr=insert_deletion_yerr, fmt='', elinewidth=0.5, capsize=2,capthick=0.5,ls='None', ecolor='black')

    # ax.bar(reg.index, reg.delrate, color='blue')
    ax0.set_title('Treatment', size=15,fontdict = {'family': 'Times New Roman'})
    ax0.set_ylabel('All Mutation (%)', size=15,fontdict = {'family': 'Times New Roman'})
    ax0.set_xticks(groupinfor.index + 1.5 * width)
    #ax0.set_xticklabels(glist, rotation=35, size=6)
    ax0.set_xticklabels(glist, rotation=35, fontdict = {'family': 'Arial'}, size = 5)
    ax0.legend((bar1[0], bar2[0], bar3[0], bar4[0]), ('mutation_all', 'deletion_only', 'insert_only','insert&&deletion'))

    bar5 = ax1.bar(groupinfor.index, ck_mutation, width, color="#CC79A7")
    bar6 = ax1.bar(groupinfor.index + width, ck_deletionO, width, color="#D55E00")
    bar7 = ax1.bar(groupinfor.index + width * 2, ck_insertO, width, color="#0072B2")
    bar8 = ax1.bar(groupinfor.index + width * 3, ck_insert_deletion, width, color="#009E73")
    # ax.bar(reg.index, reg.delrate, color='blue')
    ax1.set_title('Control', size=15,fontdict = {'family': 'Times New Roman'})
    ax1.set_ylabel('All Mutation (%)', size=15,fontdict = {'family': 'Times New Roman'})
    ax1.set_xticks(groupinfor.index + 1.5 * width)
    #ax1.set_xticklabels(ck_glist, rotation=35, size=6)
    ax1.set_xticklabels(ck_glist, rotation=35, fontdict = {'family': 'Arial'}, size = 5)
    ax1.legend((bar5[0], bar6[0], bar7[0], bar8[0]), ('mutation_all', 'deletion_only', 'insert_only', 'insert&&deletion'))
    # plt.show()
    plt.savefig(mutfile, dpi=300, format="pdf")
    plt.close(fig)
    ## print out pdf
