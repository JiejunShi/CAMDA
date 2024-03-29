# Concurrence of Active Methylation and De-methylAtion (CAMDA)
DNA methylation is introduced and maintained by DNMT family and removed by TET family. TET1 protein prevents de novo methyltransferases from binding to regulatory elements, DNMT3A also blocks TET1 binding, especially in promoter regions. Interestingly, these two ‘competing’ enzyme families are observed to be jointly associated with tumor malignancy. For example, the DNMT3A and TET2 double-knockout mice show worse survival than single-knockout counterparts; also, mutations in DNMT3A and TET2 significantly co-occur in human T-cell lymphoma. These findings suggest that the concurrence of methylation and demethylation processes is related to tumorigenesis. However, to what extent this competition contribute to cancer gene regulation remains largely unknown.

For years, DNA methylation levels are quantified in an ‘average’ manner. The increased average methylation level of CpG island (CGI), i.e., CGI hypermethylation, is a well-established mechanism for gene silencing. Besides the average methylation, DNA methylation has been quantified by its variation as ‘methylation heterogeneity’ or ‘epigenetic polymorphism'. Previous studies reveal that methylation variation is associated with global transcription variation. Moreover, higher methylation variation is linked to worse survival in acute myeloid leukemia, chronic lymphocytic leukemia, and diffuse large B-cell lymphomas, but not in glioblastoma. These studies highlight the importance of methylation variation in tumor evolution. However, neither average methylation nor methylation variation can delineate the degree of concurrence of active methylation and demethylation.

We quantified "methylation concurrence" based on the concurrence events captured by bisulfite sequencing. The methylation concurrence events are represented by the unmethylated CpGs in partially methylated reads (red circles in **Fig 1**).  
<div align=center><img src="https://github.com/JiejunShi/methylation_interruption/blob/master/images/CAMDA_schematic.png" /></div>  

**Fig 1. Schematic of methylation concurrence captured by bisulfite sequencing.**  
Bisulfite sequencing reads are dissected into three categories of fragments, i.e. methylated(***M***) fragments (consecutive solid circles in **Fig 1**), unmethylated(***U***) fragments (consecutive blank circles), and methylation-concurrence(***C***) fragments (consecutive red circles).  
<div align=center><img src="https://github.com/JiejunShi/methylation_interruption/blob/master/images/CAMDA_Equation.png" /></div>  

**Fig 2. Definition of methylation concurrence ratio.**  
The methylation concurrence ratio of a genomic region is defined as the sum of ***C*** fragments’ weights divided by the sum of all fragments’ weights in that region. Each fragment’s weight can be set as either its number of CpGs(in section **1.** below) or 1(unweighted, in section **2.** below).
‘***M***’, ‘***U***’, and ‘***C***’ represent the numbers of methylated fragments, unmethylated fragments and methylation-concurrence fragments, respectively. ‘ω<sub>m</sub>’, ‘ω<sub>u</sub>’ and ‘ω<sub>c</sub>’ are the weights for each fragment.
## Authors
- Jiejun Shi (jiejuns@uci.edu)
- Wei Li (wei.li@uci.edu)
## Dependencies
- Python3 with following packages
  - numpy
  - pandas
  - copy
  - collections
- R with following packages
  - getopt
  - dplyr
- samtools v0.1.19  
*Note: the '-X' option of 'samtools view' is required. So we suggest to use samtools 0.1.19.*
## Installation
No installation needed.
## Usage
There are two executable scripts in CAMDA toolkit, i.e. `./scripts/CAMDA.py` and `./scripts/ReadCT2CAMDA.r`. `./scripts/functions.py` contains the functions required by `./scripts/CAMDA.py`. Example files of all the input and output can be found in `./demo/`.

	$ python ./scripts/CAMDA.py
 	CAMDA Toolkit
 	For help information of each function, try:
		python CAMDA.py <Function> -h
	Availible Functions:
		CAMDA	Calculate Average Methylation Ratio(MethRatio) and Methylation Concurrence Ratio(CAMDA) of each CpG from BSMAP alignments.
		BedRatio	Calculate MethRatio or CAMDA of given regions from CpG's ratios generated by 'CAMDA' command.
		ReadCT	Generate ReadCT file from BSMAP alignments.   

	$ Rscript ./scripts/ReadCT2CAMDA.r -h
	Usage: Rscript ReadCT2CAMDA.r [-[-help|h]] [-[-ReadCT|i] <character>] [-[-Regions|r] <character>] [-[-UseStrand|s]] [-[-Weight|w] [<character>]] [-[-Output|o] [<character>]]
		-h|--help	useage
		-i|--ReadCT	ReadCT file. REQUIRED.
		-r|--Regions	Bed file of regions whose Methylation Concurrence Ratio(CAMDA) will be reported. REQUIRED.
		-s|--UseStrand	If -s is specified, strand infomation(6th column) in Regions file will be used.
		-w|--Weight	Weight applied to each sub-read fragment, either "cg" or "1". "cg" means weights equal to the CpG number of each fragment. "1" means no weight applied. [Default="cg"]
		-o|--Output	Output file report CAMDA and MethRatio of each region. [Default="Region_CAMDA.tsv"].

  - CAMDA tools take the BSMAP([by Yuanxin](https://sites.google.com/a/brown.edu/bioinformatics-in-biomed/bsmap-for-methylation)) alignments as inputs. And we suggust to sort bam file by coordinates before calculating CAMDA score.

### The calculations of weighted and unweighted versions of CAMDA are different. They are introduced below in Section 1 and 2, respectively.

### 1. CAMDA (***M***, ***U***, and ***C*** fragments are weighted by the number of CpGs)
#### 1.1 Generating Methylation Concurrence Ratio(CAMDA) and Average Methylation Ratio(MethRatio) of each CpG from BSMAP alignments

	$ python ./scripts/CAMDA.py CAMDA ./demo/example.bam hg19.fa -o example -w example -s "/path/to/samtools/v0.1.19" -x CG
	# time cost: ~5 min

  - This command will generate 4 outputs. `example_CpG_MethRatio.tsv` and `example_CpG_CAMDA.tsv` are the MethRatio and CAMDA of each CpG. MethRatio or CAMDA scores are in 5th column of the files. `example_CpG_MethRatio.wig` and `example_CpG_CAMDA.wig` are wiggle files for visulization. 

#### 1.2 Calculating CAMDA of given regions from CAMDA of CpG

	$ python ./scripts/CAMDA.py BedRatio ./demo/example.bed example_CpG_CAMDA.tsv -o example_CAMDA.tsv
	# time cost: ~1 min

#### 1.3 Calculating MethRatio of given regions from MethRatio of CpG

	$ python ./scripts/CAMDA.py BedRatio ./demo/example.bed example_CpG_MethRatio.tsv -o example_MethRatio.tsv
	# time cost: ~1 min

### 2. CAMDA (***M***, ***U***, and ***C*** fragments are unweighted)
#### 2.1 Generating ReadCT file from BSMAP alignments

	$ python ./scripts/CAMDA.py ReadCT ./demo/example.bam hg19.fa -o example_ReadCT.tsv -s "/path/to/samtools/v0.1.19" -x CG
	# time cost: ~5 min

  - ReadCT file saves all the CpG in each reads. "C" indicates methylated cytosine, and "T" indicates unmethylated cytosine. Each line is a BS-seq read. Columns of this file represent: "chr","strand","FirstCT","LastCT","CT_count","CT_pos","CT_seq".

#### 2.2 Calculating CAMDA and MethRatio from ReadCT file

	$ Rscript ./scripts/ReadCT2CAMDA.r -i example_ReadCT.tsv -r ./demo/example.bed -w 1 -o example_CAMDA_unweighted.tsv
	# time cost: ~5 min
	
  - If `-w 1`, it will generate unweighted version of CAMDA. This command can also generate weighted version of CAMDA if you replace the option `-w 1` with `-w cg`. Then in the output file `example_CAMDA_weighted.tsv`, the CAMDA value will be the same as in the output of step **1.2**. 
  - Because step **2.1** is very time-consuming for large datasets, we recommend to follow step **1.1** and **1.2** if you only want the weighted version of CAMDA, which performs better in terms of gene expression correlation.

## The CAMDA paper and its supplementary data
The methylation concurrence metric (CAMDA) was first introduced and published in paper below. 

Shi, J. *et al*. The Concurrence of DNA Methylation and Demethylation is Associated with Transcription Regulation. [***Nature Communications*** 12:5285 (2021).](https://www.nature.com/articles/s41467-021-25521-7)

The files under `./paper-data/` are the supplementary data for the above paper, which can be used to reproduce the main findings in this paper.
 
