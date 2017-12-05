# CRISPRMatch
## Brief introduction
An automatic calculation and visualization tool for high-throughput CRISPR genome-editing data analysis
## Requirements
Anaconda</br>
python3</br>
bwa</br>
samtools</br>
picard</br>

## Getting Start
1. Install Anaconda, version python3+  
[Anaconda link](https://www.anaconda.com/download/) 
[other mirror](https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/)

2. Install required softwares  
*command line*</br></br>
$ conda install bwa</br>
$ conda install samtools</br>
$ conda install picard</br>

3. Clone CRISPRMatch   
*git clone https://github.com/zhangtaolab/CRISPRMatch.git*

4. Files for mutation calculation  
- **File1**: Genome-editing target sequences  
[Fasta format example](https://github.com/zhangtaolab/CRISPRMatch/tree/master/document/Samples_gene.fa)
- **File2**: NGS samples information  
*note*:   
For CRISPR-Cas9 system, the 'Note' must contain 'gRNA' label.  
For CRISPR-Cpf1 system, the 'Note' must contain 'crRNA' label.  
*example*:  
[sample information](https://github.com/zhangtaolab/CRISPRMatch/tree/master/document/sample_infor.txt)  
- **File3**: NGS group information  
*note*: At present, two repeats are supported<br></br>
*example*:  
[group information](https://github.com/zhangtaolab/CRISPRMatch/tree/master/document/group_info.txt)

## Start running
**note**: the information files **File1**, **File2** and **File3** are required!  
</br>
**command line example**:   
*python CRISPRMatch.py -g document/Samples_gene.fa -i document/sample_infor.txt -gi document/group_info.txt*