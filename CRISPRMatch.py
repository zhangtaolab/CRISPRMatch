import argparse
import sys
from CMlib import bwa
import os
import os.path
from multiprocessing import Pool
from math import log
from pyfasta import Fasta

def main():
    args = check_options(get_options())

    # ?build bwa index
    bwaindexfile = os.path.basename(args.genome)

    bwatestindex = os.path.join(args.saved, bwaindexfile+'.sa')

    bwaindex = os.path.join(args.saved, bwaindexfile)

    bwabuild = True

    if os.path.isfile(bwatestindex):

        if not args.docker:

            print('find:', bwatestindex)

            bwamess = "Found bwa index file " + bwatestindex + ". Do you want rebuild it? Press Y or N to continue:"

            print(bwamess)

            while True:

                char = getch()

                if char.lower() in ("y", "n"):

                    print(char)

                    if char == 'y':

                        bwabuild = True

                    elif char == 'n':

                        bwabuild = False

                    break
    print("bwabuild:", bwabuild, "threads:", args.threads)
    #print("genomesize:", genomesize, "kmer:", kmer, "jfkmerfile:",
    #      jfkmerfile, "kmerbuild:", kmerbuild, "bwabuild:", bwabuild, "threads:", args.threads)

def check_options(parser):

    args = parser.parse_args()

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
    if os.path.exists(args.saved):

        if not args.docker:

            print(args.saved, "exists. Everything in this folder will be remove. Press Y or N to continue: ")

            while True:

                char = getch()

                if char.lower() in ("y", "n"):

                    print(char)

                    if char == 'n':

                        sys.exit(1)

                    break

    else:

        os.mkdir(args.saved)
    # End check saved folder

    # Print Checked information
    print("#"*40)

    print("bwa version:", args.bwa, bwaversion)

    #print("jellyfish version:", args.jellyfish, jellyfishversion)

    print("genome file:", args.genome)

    #print("input file:", args.input)

    #print("5\' labeled R primer:", args.primer)

    print("result output folder:",  os.path.realpath(args.saved))

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


def get_options():

    parser = argparse.ArgumentParser(description="CRISPRMatch is for location finding", prog='CRISPRMatch')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    #parser.add_argument('-j', '--jellyfish', dest='jellyfish', help='jellyfish path')

    parser.add_argument('-b', '--bwa', dest='bwa', help='bwa path')

    parser.add_argument('-g', '--genome', dest='genome', help='fasta format genome file', required=True)

    #parser.add_argument('-i', '--input', dest='input', help='fasta format input file', required=True)

    parser.add_argument('-s', '--save', dest='saved', help='result saved folder', default='index')

    #parser.add_argument('-p', '--primer', dest='primer', help='5\' labeled R primer', default='')

    parser.add_argument('-t', '--threads', dest='threads', help='threads number or how may cpu you wanna use',
                        default=1, type=int)

    #parser.add_argument('-l', '--length', dest='length', help='probe length', default=45, type=int)

    #parser.add_argument('--homology', dest='homology', help='homology, from 50 to 100', default=75, type=float)

    #parser.add_argument('-d', '--dtm', dest='dtm', help='dTm, from 0 to 37', default=10, type=float)

    #parser.add_argument('--step', dest='step', help='step length, min=1', default=5, type=int)

    #parser.add_argument('--docker', default=False)
    # parser.parse_args(['--version'])
    # args = parser.parse_args()

    return parser

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)
