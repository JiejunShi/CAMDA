# *M*ethylation-*I*nterruption-*E*valuated *L*ocal *D*iscordance (MIELD)
Besides global hypomethylation and focal hypermethylation, local DNA methylation discordance emerges as a new feature of tumor methylome. We quantified methylation discordance based on the methylation interruption events captured by bisulfite sequencing. The methylation interruption events are represented by the unmethylated CpGs in partially methylated reads (slashed circles in **Fig 1**).
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
## Equation of MIELD
![image](https://github.com/JiejunShi/methylation_interruption/blob/master/images/MIELD_schematic.png)  
**Fig 1**. Schematic of methylation interruption captured by bisulfite sequencing.  
Bisulfite sequencing reads are dissected into three categories of fragments, i.e. methylated fragments (consecutive solid circles in **Fig 1**), unmethylated fragments (consecutive blank circles), and methylation-interrupted fragments (consecutive slashed circles). Thus, MIELD score of a particular genomic region is measured by the following equation.  
![image](https://github.com/JiejunShi/methylation_interruption/blob/master/images/MIELD_Equation.jpg)  
‘M’, ‘U’ and ‘I’ represent the numbers of methylated fragments, unmethylated fragments and methylation-interrupted fragments, respectively. ‘ω_m’, ‘ω_u’ and ‘ω_i’ are the weights for each fragment. Optional weights can be **the CpG counts of each fragment** or **1**(unweighted).
## Usage
There are two scripts in MIELD toolkit. (`./src/functions.py` contains the functions required by `./src/MIELD.py`) Example files of all the input and output can be found in `./demo/`.

	$ python ./src/MIELD.py
 	MIELD Toolkit
 	For help information of each function, try:
		python MIELD.py <Function> -h
	Availible Functions:
		MIELD	Calculate Mean Methylation Ratio and Methylation-Interruption-Evaluated Local Discordance(MIELD) of each CpG from BSMAP alignments.
		BedRatio	Calculate Mean Methylation Ratio or MIELD of given regions from CpG's ratios generated by 'MIELD' command.
		ReadCT	Generate ReadCT file from BSMAP alignments. `  

	$ Rscript ./src/ReadCT2MIELD.r -h
	Usage: Rscript Scripts/MIELD/ReadCT2MIELD.r [-[-help|h]] [-[-ReadCT|i] <character>] [-[-Regions|r] <character>] [-[-UseStrand|s]] [-[-Weight|w] [<character>]] [-[-Output|o] [<character>]]
		-h|--help	useage
		-i|--ReadCT       ReadCT file. REQUIRED.
		-r|--Regions      Bed file of regions whose Methylation-Interruption-Evaluated Local Discordance(MIELD) will be reported. REQUIRED.
		-s|--UseStrand    If -s is specified, strand infomation(6th column) in Regions file will be used.
		-w|--Weight       Weight applied to each sub-read fragment, either "cg" or "1". "cg" means weight equal to the length of fragment. "1" means no weight applied. [Default="cg"]
		-o|--Output       Output file report MIELD and MethRatio of each region. [Default="Region_MIELD.tsv"].

  - MIELD tools take the BSMAP([by Yuanxin](https://sites.google.com/a/brown.edu/bioinformatics-in-biomed/bsmap-for-methylation)) alignments as inputs. And we suggust to sort bam file by coordinates before calculating MIELD score.

### 1. MIELD (weights equal to the CpG counts of each fragment)
#### 1.1 Generating MIELD and Mean Methylation Ratio(MethRatio) of each CpG from BSMAP alignments

	$ python ./src/MIELD.py MIELD ./demo/example.bam hg19.fa -o example -w example -s "/path/to/samtools/v0.1.19" -x CG

  - This command will generate 4 outputs. **example_CpG_MethRatio.tsv** and **example_CpG_MIELD.tsv** are the MethRatio and MIELD of each CpG. MethRatio or MIELD scores are in 5th column of the files. **example_CpG_MethRatio.wig** and **example_CpG_MIELD.wig** are wiggle files for visulization. 

#### 1.2 Calculating MIELD ratio of given regions from MIELD scores of CpG

	$ python ./src/MIELD.py BedRatio ./demo/example.bed example_CpG_MIELD.tsv -o example_MIELD.tsv

#### 1.3 Calculating MethRatio of given regions from MethRatio of CpG

	$ python ./src/MIELD.py BedRatio ./demo/example.bed example_CpG_MethRatio.tsv -o example_MethRatio.tsv

### 2. MIELD (weights equal to 1 for each fragment)
#### 2.1 Generating ReadCT file from BSMAP alignments

	$ python ./src/MIELD.py ReadCT ./demo/example.bam hg19.fa -o example_ReadCT.tsv -s "/path/to/samtools/v0.1.19" -x CG

  - ReadCT file saves all the CpG in each reads. "C" indicates methylated cytosine, and "T" indicates unmethylated cytosine. Each line is a BS-seq read. Columns of this file represent: "chr","strand","FirstCT","LastCT","CT_count","CT_pos","CT_seq".

#### 2.2 Calculating MIELD and MethRatio from ReadCT file

	$ Rscript ./src/ReadCT2MIELD.r -i example_ReadCT.tsv -r ./demo/example.bed -w 1 -o example_MIELD_weight_1.tsv
	
  - If `-w 1`, MIELD scores with weights equal to 1 are generated. If `-w cg`, MIELD scores with weights equal to CpG counts are generated. In the 2nd case, the MIELD score will be the same as the in step **1.2**. 
  - Because step **2.1** is time-consuming, we recommend to follow step **1.1** and **1.2** if you only want the MIELD score with weights equal to CpG counts, which performs better in terms of expression correlation.

