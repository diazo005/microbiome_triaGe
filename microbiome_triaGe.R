# Title: microbiome_triaGe.R
# Description: Script to preliminary explore a microbiome dataset from a longitudinal study
# Author: Gerardo R Diaz
# Date: 12/15/23
# sidenote (How to save ps_object as rds): saveRDS(myPS, 'myPS.rds')
# Load parameters: CHANGE THE LOCATION OF YOUR PARAMETER FILE
params <- read.csv('/Users/gerardodiaz/Desktop/Gerardo/UMN/PhD/1 Fall 2023/CSCI5465_Computing Biology/Final project/parameters.csv')

##### 1. Load info #####
# Load packages
library('phyloseq')
library("tidyverse")
library('gridExtra')
library('metagenomeSeq')
library('vegan')
library('table1')
library('microbiome')
# Set working directory 
setwd(params$wd)
# Load dataset
ps <- readRDS(params$ps_object)
# assigning variables from parameters
tx <- params$tx_variable
tp <- params$time_variable
sample_data(ps)[[tx]] <- factor(sample_data(ps)[[tx]])
sample_data(ps)[[tp]] <- factor(sample_data(ps)[[tp]])

##### 2. Data formatting ##### 
# Subsetting only samples (only ones part of the trial)
ps_controls <- subset_samples(ps, is.na(sample_data(ps)[[tx]]))
ps <- subset_samples(ps, !is.na(sample_data(ps)[[tx]]))
# Subsetting only archaea and bacteria
if (params$X16S_shotgun == 'shotgun') {
  ps <- subset_taxa(ps, Domain %in% c("Bacteria","Archaea"))
} else {
  ps <- subset_taxa(ps, !Kingdom %in% c("NA"))
}
# Normalization
metaseq <- phyloseq_to_metagenomeSeq(ps)
metaseq.norm<- cumNorm(metaseq, p=cumNormStat(metaseq))
CSS_metaseq <- MRcounts(metaseq.norm, norm = TRUE)
ps_CSS <- merge_phyloseq(otu_table(CSS_metaseq,taxa_are_rows=T),sample_data(ps),tax_table(ps))
# Agglomeration 
ranks <- strsplit(params$tax_ranks, ",")[[1]]
for (i in ranks) {
  object_name <- paste0("ps_CSS_",i)
  assign(object_name, tax_glom(ps_CSS, i))
}

##### 3. Relative abundance plots #####
# Formatting long for ggplot
ps_CSS_objects <- ls()
ps_objects <- grep("^ps_CSS_", ps_CSS_objects, value = TRUE)
for (i in ps_objects) {
  # obtaining the taxonomic rank to pass in other functions
  name_parts <- strsplit(i, "_")[[1]]
  tax_rank <- tail(name_parts, 1)
  # obtaining dynamic object name for long 
  object_name <- paste0("long_",tax_rank)
  # getting the object from environment
  ps_name <- get(i)
  # passing the object into long formatting for ggplot
  relative <- transform_sample_counts(ps_name, function(x) (x / sum(x))*100 )
  relative_long <- psmelt(relative)
  relative_long <- relative_long %>%
    group_by(!!as.name(tax_rank)) %>%
    mutate(mean_relative_abund = mean(Abundance))
  relative_long[[tax_rank]] <- as.character(relative_long[[tax_rank]])
  relative_long$mean_relative_abund <- as.numeric(relative_long$mean_relative_abund)
  relative_long[[tax_rank]][relative_long$Abundance < 5] <- "Others (< 5%)"
  assign(object_name, relative_long)
  # EXTRA: Obtaining a table for relative abundances (Comment out if not needed)
  print(table1(~ Abundance | relative_long[[tax_rank]], data=relative_long, transpose = T))
  # plotting to a png file 
  plot<-relative_long %>%
    ggplot(aes(x = Sample, y = Abundance, fill = !!as.name(tax_rank))) +
    geom_bar(stat = "identity",width = 1,color="black") +
    geom_bar(stat = "identity", alpha=1)+theme_bw()+
    facet_wrap(as.formula(paste("~",tp, "+",tx)), nrow=1, scales="free_x")+
    theme(axis.text.x = element_blank())+
    labs(title = paste0(tax_rank,'-level RA plot stratified by time and treatment'), x= 'Samples', y="Relative abundance (%)")+
    ggthemes::scale_fill_tableau("Tableau 20")
  ggsave(paste0('RA_', tax_rank, '_level', '.png'), plot,width = 10, height = 6, dpi=300)
}

##### 4. Alpha diversity #####
# Agglomerating to different taxonomic ranks
ranks <- strsplit(params$tax_ranks, ",")[[1]]
Alpha_list <- list()
for (i in ranks) {
  # agglomerating to different tax ranks
  glom <- tax_glom(ps, i)
  div <- estimate_richness(glom, measures=c("Observed", "Shannon"))
  Alpha <- cbind(div,(evenness(glom, 'pielou')))
  names(Alpha) <- paste0(names(Alpha),'_', i)
  # creating dataframe with alpha-div indices
  Alpha_df <- cbind(sample_data(ps),Alpha)
  # plotting richness
  plot1 <- ggplot(Alpha_df, aes(x = !!as.name(tp), y = !!as.name(paste0('Observed_',i)), color=!!as.name(tx)))+
    geom_boxplot(position=position_dodge(0.8),lwd=1) +
    geom_jitter(position=position_dodge(0.8), size=2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+ 
    labs(title = paste0(i,'-level Richness stratified by time and treatment'), x = "Time-point", y = "Richness")+  
    guides(color=guide_legend(title="Treatment"))
  # plotting shannon's index
  plot2 <- ggplot(Alpha_df, aes(x = !!as.name(tp), y = !!as.name(paste0('Shannon_',i)), color=!!as.name(tx)))+
    geom_boxplot(position=position_dodge(0.8),lwd=1) +
    geom_jitter(position=position_dodge(0.8), size=2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+ 
    labs(title = paste0(i,"-level Shannon's index stratified by time and treatment"), x = "Time-point", y = "Shannon's index")+  
    guides(color=guide_legend(title="Treatment"))
  # plotting pielou's index
  plot3 <- ggplot(Alpha_df, aes(x = !!as.name(tp), y = !!as.name(paste0('pielou_',i)), color=!!as.name(tx)))+
    geom_boxplot(position=position_dodge(0.8),lwd=1) +
    geom_jitter(position=position_dodge(0.8), size=2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+ 
    labs(title = paste0(i,"-level Pielou's evenness stratified by time and treatment"), x = "Time-point", y = "Pielou's evenness")+  
    guides(color=guide_legend(title="Treatment"))
  # combining the 3 plots
  plots <- grid.arrange(plot1, plot2, plot3, nrow = 3)
  # plotting the 3 plots to a png file
  ggsave(paste0(i,'-level_alpha_div_plots','.png'), plots,width = 6, height = 12, dpi=300)
  # appending all the alpha-div at different tax ranks
  Alpha_list[[i]] <- Alpha
}
# outputting table csv with the dataframe and the alpha-div indices
Alpha_df <- do.call(cbind, Alpha_list)
Alpha_df <- cbind(sample_data(ps), Alpha_df)
write.csv(Alpha_df, "df_with_alpha_div.csv")

##### 5. Beta diversity #####
for (i in ps_objects) {
  # obtaining the taxonomic rank to pass in other functions
  name_parts <- strsplit(i, "_")[[1]]
  tax_rank <- tail(name_parts, 1)
  # getting the object from environment
  ps_name <- get(i)
  set.seed(143)
  # calculating ordination
  ordinated <- ordinate(ps_name, "NMDS", "bray")
  # plotting to a png file
  Beta<-plot_ordination(ps_name, ordinated, type="sample", color=tx)+
    geom_point(size=4.5)+stat_ellipse()+theme_classic()+labs(title = paste0(tax_rank,"-level Beta div stratified by time"))+
    facet_wrap(as.formula(paste("~",tp)), nrow=1, scales="free_x")
  Beta$layers
  Beta$layers <- Beta$layers[-1]
  ggsave(paste0(tax_rank,'-level_beta_div_plot','.png'), Beta ,width = 12, height = 6, dpi=300)
}
