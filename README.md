# microbiome_triaGe
R script for microbiome data visualization 
The goal of this project is to develop a convenient R Script for an initial visual exploration of animal microbiome data from **longitudinal studies**. It will provide researchers with crucial microbiome plots to develop intuition of the data patterns present in their datasets and guide their further statistical analysis. 
Questions: diazo005@umn.edu (Gerardo R Diaz)

## Package and data types 
The script requires the installation of the following packages:
- Phyloseq
- Tidyverse
- gridExtra
- metagenomeSeq
- vegan
- table1
- microbiome

*The script accepts **.rds** files as input and **.csv** files as parameters.*
I suggest using *DADA2* and *Kraken2* to produce the count matrices for 16S rRNA gene sequencing and shotgun metagenomic sequencing, respectively.

## Input
The script requires 3 items:
- **Parameters file:** A csv file that requires the (1) name of the project in the 1st column; (2) the dataset type in the 2nd column (either ‘16S’ or ‘shotgun’); (3) the phyloseq object name as rds file in the 3rd column (e.g., my_ps_object.rds) ; (4) the taxonomic ranks in which the user wants to analyze the data (they need to be introduced with the first capital letter and comma without spaces. E.g., Phylum,Order,Genus) in the 4th column; (5) the working directory address where the input files are standing and the output files will be stored, in the 5th column; (6) the name of the column in the dataset that has the time variable, in the 6th column; (7) the name of the column in the dataset that has the treatment variable, in the 7th column (figure 1).
- **Phyloseq object:** This object is produced with phyloseq package in R and contains the otu table, taxa table and sample data combined as a multidimensional matrix. The phyloseq object needs to be saved only as rds file in order to be appropriate as input for the script. 
- **Script modification:** The script Microbiome_triaGe_Diaz.R requires a simple modification in line 7 with the address where the parameters.csv file is located. It may be the folder that will be used as working directory (Figure 2).
## How to use it
1. Copy the **microbiome_triaGe.R** script and the **parameters.csv** files to your working directory.
2. Copy your phyloseq object saved as **.rds** file to your working directory (more on how to build phyloseq objects: https://vaulot.github.io/tutorials/Phyloseq_tutorial.html)
       You should have something like this: ![alt text](/figures/Microbiome1.png "You need 3 files in your working directory!")
4. Modify **parameters.csv** file as follows:
    - Column 1: Name of the project
    - Column 2: Specify dataset type. Write either *16S* or *shotgun*, for 16S rRNA sequencing or Shotgun Metagenomic sequencing. 
    - Column 3: Write the name of your phyloseq object saved as .rds file. For instance, *ps.long.rds* or *ps.rds*.
    - Column 4: Write the taxonomic ranks in which the figures will be generated. They should be written with the first capital letter and comma without spaces. For instance, *Phylum,Order,Genus* or *Order,Class,Genus*.
    - Column 5: Write the address for your working directory using the forward slash. Something like this: */Path/to/my/working/directory*
    - Column 6: Write the name of the column that contains the timepoints in your dataset. In other words your time variable name.
    - Column 7: Write the name of the column that contains the variable of interest in your dataset. In other words your treatment variable name.
    Here is an example of how to fill it: ![alt text](/figures/Microbiome2.png "You may have something like this!")
5. Modify **microbiome_triaGe.R** script: Write in line 7 the address where your **parameters.csv** is located.
   For example: ![alt text](/figures/Microbiome3.png "The address is in light green")

## Output
The script will output 4 types of files:
1.	Relative abundance plots stratified by time and treatment variable (png files at 300dpi resolution).
    Something like this: ![alt text](/figures/Microbiome4.png "Relative abundance plots example")
2.	Alpha diversity indices plotted as boxplots, stratified by time and treatment variable (png files at 300dpi resolution).
    Something like this: ![alt text](/figures/Microbiome5.png "Relative abundance plots example")
3.	Beta diversity plots stratified by time and treatment variable (png files at 300dpi resolution).
    Something like this: ![alt text](/figures/Microbiome6.png "Relative abundance plots example")
4.	Metadata of the microbiome dataset with the alpha indices calculated added to it (csv file).
    Something like this: ![alt text](/figures/Microbiome7.png "Relative abundance plots example")

Please note that a relative abundance table of top taxa will be produced in the viewer of RStudio. Also, all these output items will be located in the working directory provided as an address in the parameter file.
You may find your working directory populated with the plots. Something like this: ![alt text](/figures/Microbiome8.png "The address is in light green")

## But...How does the script work?
The script will conduct 5 tasks: 
1.	Load info: It will load the R packages required to work and will set the working directory, read the phyloseq object as rds file, and read the time and treatment variables names for plotting.
2.	Data formatting: It will subset only Bacteria and archaea from the microbiome dataset, normalize the data by the Cumulative Sum Scaling (CSS) method, and agglomerate the data to the taxonomic ranks that the user specified in the parameters.csv file for plotting. 
3.	Relative abundance plot: It will convert the normalized phyloseq objects to long format and add the adequate columns to create bar plots of relative abundance at the taxonomic ranks specified by the user. The plots will be stratified by time and treatment variable. It will also produce summary tables with the abundances of top taxa across samples in the viewer of RStudio. Note that the relative abundance plots will be outputted as png files at 300dpi resolution in the working directory. 
4.	Alpha diversity plot: it will agglomerate the non-normalized microbiome data to different taxonomic ranks specified by the user and calculate 3 alpha diversity indices at each taxonomic rank: Richness, Shannon’s index and Pielou’s evenness. It will produce combined boxplots of the 3 diversity indices at each taxonomic rank and output them as png files at 300dpi resolution in the working directory. The plots will be stratified by time and treatment variable. Note that it will also produce a csv file named df_with_alpha_div.csv that will contain the microbiome samples data (metadata) combined with the alpha diversity indices calculated previously.
5.	Beta diversity plot: It will calculate the ordination of the normalized data at each taxonomic rank and produce beta diversity plots stratified by time and treatment variable. The figures will be outputted as png files at 300dpi resolution in the working directory. 

