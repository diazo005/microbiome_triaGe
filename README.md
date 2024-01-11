# microbiome_triaGe
R script for microbiome data visualization 
The goal of this project is to develop a convenient R Script for an initial visual exploration of animal microbiome data from longitudinal studies.  It will provide researchers with crucial microbiome plots to develop intuition of the data patterns present in their datasets and guide their further statistical analysis.

## Package and data types 
The script requires the installation of the following packages:
- Phyloseq
- Tidyverse
- gridExtra
- metagenomeSeq
- vegan
- table1
- microbiome

*The script accepts rds as input file and csv as parameters file.*

## How to use it
### Input
The script requires 3 items:
Parameters file: A csv file that requires the (1) name of the project in the 1st column; (2) the dataset type in the 2nd column (either ‘16S’ or ‘shotgun’); (3) the phyloseq object name as rds file in the 3rd column (e.g., my_ps_object.rds) ; (4) the taxonomic ranks in which the user wants to analyze the data (they need to be introduced with the first capital letter and comma without spaces. E.g., Phylum,Order,Genus) in the 4th column; (5) the working directory address where the input files are standing and the output files will be stored, in the 5th column; (6) the name of the column in the dataset that has the time variable, in the 6th column; (7) the name of the column in the dataset that has the treatment variable, in the 7th column (figure 1).

### Output
