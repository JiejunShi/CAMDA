# *M*ethylation-*I*nterruption-*E*valuated *L*ocal *D*iscordance (MIELD)
Besides global hypomethylation and focal hypermethylation, local DNA methylation discordance emerges as a new feature of tumor methylome. We quantified methylation discordance based on the methylation interruption events captured by bisulfite sequencing. The methylation interruption events are represented by the unmethylated CpGs in partially methylated reads (slashed circles in Fig 1). 
## Equation of MIELD
![image](https://github.com/JiejunShi/methylation_interruption/blob/master/images/MIELD_schematic.png)
Bisulfite sequencing reads are dissected into three categories of fragments, i.e. methylated fragments (consecutive solid circles in Fig 1), unmethylated fragments (consecutive blank circles), and methylation interrupted fragments (consecutive slashed circles). Thus, MIELD score of a particular genomic region is measured by the following equation. 

![image](https://github.com/JiejunShi/methylation_interruption/blob/master/images/MIELD_Equation.jpg)
M, U and I represent the numbers of methylated fragments, unmethylated fragments and methylation interrupted fragments, respectively. ω_m, ω_u and ω_i are the weights for each fragment. Optional weights can be the CpG counts of each fragment or 1.
