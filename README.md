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
‘M’, ‘U’ and ‘I’ represent the numbers of methylated fragments, unmethylated fragments and methylation-interrupted fragments, respectively. ‘ω_m’, ‘ω_u’ and ‘ω_i’ are the weights for each fragment. Optional weights can be the CpG counts of each fragment or 1(unweighted).
## Usage
**MIELD** toolkit has three functions:  
	$ python ./src/MIELD.py -h
 	MIELD Toolkit
 	For help information of each function, try:
		python MIELD.py <Function> -h
	Availible Functions:
		MIELD	Calculate Mean Methylation Ratio and Methylation-Interruption-Evaluated Local Discordance(MIELD) of each CpG from BSMAP alignments.
		BedRatio	Calculate Mean Methylation Ratio or MIELD of given regions from CpG's ratios generated by 'MIELD' command.
		ReadCT	Generate ReadCT file from BSMAP alignments. `  

**MIELD** tools take the BSMAP generated bam files as inputs. And we suggust to sort bam file by coordinates before calcuting MIELD score.
1. Generating Mean Methylation Ratio and MIELD of each CpG from BSMAP alignments
	$ python ./src/MIELD.py MIELD -h
	
