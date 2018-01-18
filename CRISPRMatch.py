import argparse
import sys
from CMlib import bwa
import os
import os.path
from multiprocessing import Pool
from math import log
from pyfasta import Fasta
import pandas as pd
from CMlib import bwa_run
from CMlib import mut_rate
from CMlib import output_aln_fa
from CMlib import plot_each_bam
from CMlib import plot_pdf
from CMlib import output_aln_pdf
from subprocess import Popen
from subprocess import PIPE
import re

def main():
    args = check_options(get_options())

    # ?build bwa index
    bwaindexfile = os.path.basename(args.genome)

    bwatestindex = os.path.join(args.saved, bwaindexfile+'.sa')

    bwaindex = os.path.join(args.saved, bwaindexfile)




    bwabuild = True

    # if os.path.isfile(bwatestindex):
    #
    #     if not args.docker:
    #
    #         print('find:', bwatestindex)
    #
    #         bwamess = "Found bwa index file " + bwatestindex + ". Do you want rebuild it? Press Y or N to continue:"
    #
    #         print(bwamess)
    #
    #         while True:
    #
    #             char = getch()
    #
    #             if char.lower() in ("y", "n"):
    #
    #                 print(char)
    #
    #                 if char == 'y':
    #
    #                     bwabuild = True
    #
    #                 elif char == 'n':
    #
    #                     bwabuild = False
    #
    #                 break
    print("bwabuild:", bwabuild, "threads:", args.threads)
    #print("genomesize:", genomesize, "kmer:", kmer, "jfkmerfile:",
    #      jfkmerfile, "kmerbuild:", kmerbuild, "bwabuild:", bwabuild, "threads:", args.threads)

    # ?Build Gene_fasta index

    if bwabuild:
    
        bwa.bwaindex(args.bwa, args.genome, args.saved)
        print("## Step 1:")
    
        print("bwa index build finished ...")
    
    else:
    
        print("Use", bwatestindex)
    
    
    print("bwa index finshed!!")
    
    # ?run bwa alignment
    print("## Step 2:")
    print("loading fastq files...!")
    bwa_run.prepare(args.input, args.genome, args.saved, args.bwa, args.samtools, args.picard)
    print("bwa mem finished!")
    # end run bwa alignment


    # ?mutation ration calculation
    print("## Step 3:")
    mut_rate.rate_cal(args.input, args.groupinfo, args.genome, args.result, args.saved)
    print("Mutation calculation finished!")
    # end mutation


    # ?mutation result display
    mut_rate.display(args.groupinfo, args.result)
    print("Mutation calculation result have been displayed!")
    # end display


    # ?output aln and fa file
    print("## Step 4:")
    output_aln_fa.alnfile(args.input, args.groupinfo, args.genome, args.result, args.saved)
    print("Alignment files were output!")
    # end output aln and fa file
    
    
    # ?output aln figure
    print("Starting to plot each alignment...")
    output_aln_fa.alnpdf(args.input, args.result)
    print("Alignment figures were done!")
    # end output aln figure
    
    
    # ?plot each bam
    print("## Step 5:")
    print("Starting to plot each bam...")
    plot_each_bam.barchart(args.input, args.groupinfo,args.genome, args.result, args.saved)
    print("plot each bam finished!")
    # end plot each bam
    
    # ?plot pdf
    print("## Step 6:")
    print("Starting to plot pdf...")
    plot_pdf.plotpdf(args.groupinfo, args.genome, args.result, args.saved)
    print("plot pdf finished!")
    # end plot pdf

    # # ?output aln and fa file
    # print("## Step test:")
    # output_aln_pdf.alnpdftest(args.input, args.result, args.genome,args.groupinfo)
    # print("Alignment files were output!")
    # # end output aln and fa file


def check_options(parser):

    args = parser.parse_args()

    # Start check samtools
    if args.samtools:

        if not os.path.exists(args.samtools):

            print("Can not locate samtools, please input full path of samtools\n")

            parser.print_help()

            sys.exit(1)

    else:

        samtoolspath = which('samtools')

        if samtoolspath:

            samtoolsversion=samtools('samtools')
            if samtoolsversion == 'None':

                print("Can not locate samtools, please input full path of samtools\n")

                parser.print_help()

                sys.exit(1)

            else:

                args.samtools = samtoolspath[0]

        else:

            print("Can not locate samtools, please input full path of samtools\n")

            parser.print_help()

            sys.exit(1)

    # End check samtools

    # Start check picard
    if args.picard:

        if not os.path.exists(args.picard):

            print("Can not locate picard, please input full path of picard\n")

            parser.print_help()

            sys.exit(1)

    else:

        picardpath = which('picard')

        if picardpath:

            picardversion=picard('picard')
            if picardversion == 'None':

                print("Can not locate picard, please input full path of picard\n")

                parser.print_help()

                sys.exit(1)

            else:

                args.picard = picardpath[0]

        else:

            print("Can not locate picard, please input full path of picard\n")

            parser.print_help()

            sys.exit(1)

    # End check picard

    # Start check bwa
    if args.bwa:

        if not os.path.exists(args.bwa):

            print("Can not locate bwa, please input full path of bwa\n")

            parser.print_help()

            sys.exit(1)

        bwaversion = bwa.bwaversion(args.bwa)

        if bwaversion == 'None':

            print("Can not locate bwa, please input full path of bwa\n")

            parser.print_help()

            sys.exit(1)

    else:

        bwapath = which('bwa')

        if bwapath:

            bwaversion = bwa.bwaversion(bwapath[0])

            if bwaversion == 'None':

                print("Can not locate bwa, please input full path of bwa\n")

                parser.print_help()

                sys.exit(1)

            else:

                args.bwa = bwapath[0]

        else:

            print("Can not locate bwa, please input full path of bwa\n")

            parser.print_help()

            sys.exit(1)

    # End check bwa

    if not os.path.exists(args.genome):

        print("Can not locate genome file, please input genome file.\n")

        parser.print_help()

        sys.exit(1)

    # Start check saved folder
    if not os.path.exists(args.saved):

        os.mkdir(args.saved)

    #End check saved folder

    # Start check result folder
    if not os.path.exists(args.result):
        os.mkdir(args.result)

    # End check result folder

    # # Start check saved folder
    # if os.path.exists(args.saved):
    #
    #     if not args.docker:
    #
    #         print(args.saved, "exists. Everything in this folder will be remove. Press Y or N to continue: ")
    #
    #         while True:
    #
    #             char = getch()
    #
    #             if char.lower() in ("y", "n"):
    #
    #                 print(char)
    #
    #                 if char == 'n':
    #
    #                     sys.exit(1)
    #
    #                 break
    #
    # else:
    #
    #     os.mkdir(args.saved)
    # End check saved folder

    # Print Checked information
    print("#"*40)

    print("bwa version:", args.bwa, bwaversion)

    print("samtools version:", args.samtools, samtoolsversion)

    print("picard version:", args.picard, picardversion)

    #print("jellyfish version:", args.jellyfish, jellyfishversion)

    print("genome file:", args.genome)

    #print("input file:", args.input)

    #print("5\' labeled R primer:", args.primer)

    print("tmp output folder:",  os.path.realpath(args.saved))
    print("result output folder:", os.path.realpath(args.result))

    print("threads number:", args.threads)

    #print("homology:", args.homology)

    #print("dtm:", args.dtm)

    print("#"*40)

    return args

def getch():
    """
    For yes/no choice
    """
    import sys, tty, termios
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch


def which(filename):
    """docstring for which"""
    locations = os.environ.get("PATH").split(os.pathsep)
    candidates = []
    for location in locations:
        candidate = os.path.join(location, filename)
        if os.path.isfile(candidate):
            candidates.append(candidate)
    return candidates

def samtools(filename):
    """
    :param filename:
    :return: samtools version
    """
    samtoolspath=which(filename)
    samtoolscmd = ' '.join([samtoolspath[0], '--version'])
    #location= samtoolspath[0]
    samtoolsrun = Popen(samtoolscmd, stdout=PIPE, stderr=PIPE, shell=True)
    i=samtoolsrun.stdout.readlines()[0]
    version = i.decode('utf-8').rstrip('\n')
    samtoolsrun.communicate()
    return version

def picard(filename):
    """
    :param filename:
    :return:
    """
    picardpath=which(filename)
    picardcmd = ' '.join([picardpath[0], 'ViewSam', '-h'])
    version = 'None'
    picardrun = Popen(picardcmd, stdout=PIPE, stderr=PIPE, shell=True)
    #print(picardcmd)
    for i in picardrun.stderr.readlines():

        i = i.decode('utf-8').rstrip('\n')

        if re.search('Version', i):
            (_, version) = i.split(' ')
            print(version)

    picardrun.communicate()
    return version

def get_options():

    parser = argparse.ArgumentParser(description="CRISPRMatch is for location finding", prog='CRISPRMatch')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    #parser.add_argument('-j', '--jellyfish', dest='jellyfish', help='jellyfish path')

    parser.add_argument('-b', '--bwa', dest='bwa', help='bwa path')

    parser.add_argument('-sm', '--samtools', dest='samtools', help='samtools path')

    parser.add_argument('-pi', '--picard', dest='picard', help='picard path')

    parser.add_argument('-g', '--genome', dest='genome', help='fasta format genome file', required=True)

    parser.add_argument('-i', '--input', dest='input', help='sample information input file', required=True)

    parser.add_argument('-gi', '--groupinfo', dest='groupinfo', help='group information input file', required=True)

    parser.add_argument('-s', '--save', dest='saved', help='tmp saved folder', default='tmpfiles')

    parser.add_argument('-r', '--result', dest='result', help='result saved folder', default='result')

    #parser.add_argument('-p', '--primer', dest='primer', help='5\' labeled R primer', default='')

    parser.add_argument('-t', '--threads', dest='threads', help='threads number or how may cpu you wanna use',
                        default=1, type=int)

    #parser.add_argument('-l', '--length', dest='length', help='probe length', default=45, type=int)

    #parser.add_argument('--homology', dest='homology', help='homology, from 50 to 100', default=75, type=float)

    #parser.add_argument('-d', '--dtm', dest='dtm', help='dTm, from 0 to 37', default=10, type=float)

    #parser.add_argument('--step', dest='step', help='step length, min=1', default=5, type=int)

    parser.add_argument('--docker', default=False)
    # parser.parse_args(['--version'])
    # args = parser.parse_args()

    return parser

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)
