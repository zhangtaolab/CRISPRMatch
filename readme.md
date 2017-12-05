# CRISPRMatch
## Brief introduction
An automatic calculation and visualization tool for high-throughput CRISPR genome-editing data analysis
## Requirements
Anaconda</br>
python3</br>
bwa</br>
samtools</br>
picard</br>

## Manually Install
CentOS Linux release 7.3.1611 (terminal)
1. Install Anaconda</br>
$ yum install wget git</br>
$ mkdir /home/software</br>
$ cd /home/software</br>
$ wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh</br>
$ bash Anaconda3-5.0.1-Linux-x86_64.sh</br>
[Anaconda link](https://www.anaconda.com/download/)</br>
[other mirror](https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/)

2. Install dependent package  
$ conda install bwa \  
samtools \  
picard \  
matplotlib \  
pysam \  
pandas \  
argparse \  
numpy \  

3. Download CRISPRMatch  
$ cd /home/software</br>
$ git clone https://github.com/zhangtaolab/CRISPRMatch.git</br>
$ python3 /home/software/CRISPRMatch/CRISPRMatch.py -h</br>

```
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

## Start running
1. Files for mutation calculation  
- **File1**: Genome-editing target sequences  
[Fasta format example](https://github.com/zhangtaolab/CRISPRMatch/tree/master/document/Samples_gene.fa)
- **File2**: NGS samples information  
*note*:   
For CRISPR-Cas9 system, the 'Note' must contain 'gRNA' label.  
For CRISPR-Cpf1 system, the 'Note' must contain 'crRNA' label.  
*example*:  
[sample information](https://github.com/zhangtaolab/CRISPRMatch/tree/master/document/sample_infor.txt)  
- **File3**: NGS group information  
*note*: At present, two repeats are supported<br>
*example*:</br>
[group information](https://github.com/zhangtaolab/CRISPRMatch/tree/master/document/group_info.txt)  
**note**: the information files **File1**, **File2** and **File3** are required!  
</br>
2. command line example:</br>

```
python CRISPRMatch.py -g document/Samples_gene.fa -i document/sample_infor.txt -gi document/group_info.txt
```
```aaa
```