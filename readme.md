# CRISPRMatch
## Brief introduction
An automatic calculation and visualization tool for high-throughput CRISPR genome-editing data analysis
## I. <u>Requirements</u>
Anaconda</br>
python3</br>
bwa</br>
samtools</br>
picard</br>
FLASH</br>

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) **<font color=red>Note:</font>** Using `Anaconda` to Install all packages (`bwa,samtools,picard,FLASH`)

## II. <u>Manually Install</u>
CentOS Linux release 7.3.1611 (terminal)
1. Install Anaconda</br>
```
$ yum install wget git
$ mkdir /home/software
$ cd /home/software
$ wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
$ bash Anaconda3-5.0.1-Linux-x86_64.sh
```
2. Install required packages  
```
$ conda install bwa \  
                samtools \  
                picard \  
                flash \ 
                matplotlib \  
                pysam \  
                pandas \  
                argparse \  
                numpy \
```
- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) **<font color=red>Note:</font>** To ensure the tool working, please using `Anaconda` to install all packages (`bwa,samtools,picard,FLASH ...`)

3. Download CRISPRMatch and test
```
$ cd /home/software
$ git clone https://github.com/zhangtaolab/CRISPRMatch.git
$ python3 /home/software/CRISPRMatch/CRISPRMatch.py -h
  
  usage: CRISPRMatch [-h] [--version] [-b BWA] [-sm SAMTOOLS] [-pi PICARD] -g
                    GENOME -i INPUT -gi GROUPINFO [-s SAVED] [-r RESULT]
                    [-t THREADS] [--docker DOCKER]

  CRISPRMatch is for location finding

  optional arguments:
    -h, --help            show this help message and exit
    --version             show program's version number and exit
    -b BWA, --bwa BWA     bwa path
    -sm SAMTOOLS, --samtools SAMTOOLS
                          samtools path
    -pi PICARD, --picard PICARD
                          picard path
    -g GENOME, --genome GENOME
                          fasta format genome file
    -i INPUT, --input INPUT
                          sample information input file
    -gi GROUPINFO, --groupinfo GROUPINFO
                          group information input file
    -s SAVED, --save SAVED
                          tmp saved folder
    -r RESULT, --result RESULT
                          result saved folder
    -t THREADS, --threads THREADS
                          threads number or how may cpu you wanna use
```

## III. <u>Start running</u>
1. Files for mutation calculation  
- **File1**: Genome-editing target sequences  
[Fasta format example](https://github.com/zhangtaolab/CRISPRMatch/blob/master/sampledata/Samples_gene.fa)
- **File2**: NGS samples information  
![#f03c15](https://placehold.it/15/f03c15/000000?text=+) <font color=red>*note*:</font>   
For CRISPR-Cas9 system, the `'Note'` must contain `'gRNA'` label.  
For CRISPR-Cpf1 system, the `'Note'` must contain `'crRNA'` label.  
*example*:  
[sample information](https://github.com/zhangtaolab/CRISPRMatch/blob/master/sampledata/sample_infor.txt)  
- **File3**: NGS group information  
![#f03c15](https://placehold.it/15/f03c15/000000?text=+) <font color=red>*note*:</font> At present, two repeats are supported<br>
*example*:</br>
[group information](https://github.com/zhangtaolab/CRISPRMatch/blob/master/sampledata/group_info.txt)  
- **Note**: the information files `File1`, `File2` and `File3` are required!  
</br>
2. command line example:</br>

> (1) For single long reads

```
$ cd /home/software/CRISPRMatch/
$ python3 CRISPRMatch.py -g sampledata/Samples_gene.fa -i sampledata/sample_infor.txt -gi sampledata/group_info.txt -t 2
```
```diff
- Note: absolute path is preferred when using customer data
```
> (2) For paired-end reads
```
$ cd /home/software/CRISPRMatch/
$ python3 CRISPRMatch_paired.py -g sampledata2/Samples_gene.fa -i sampledata2/sample_infor.txt -gi sampledata2/group_info.txt -t 2
```
```diff
- Note: absolute path is preferred when using customer data
```