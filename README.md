# Re-annotation of the hyperthermophilic archaeon *Pyrococcus furiosus* using genome and RNA sequencing  
Felix Grünberger (Institute of Microbiology and Archaea Center, University of Regensburg, Regensburg, Germany)  

This Github repository contains code, figures and documentation for reproducing our work and may be used for future similar projects.  

## Structure  
The repository is structured in the folders `code` (R code and Rmd documents), `data` (genomic data and output from ANNOgesic/ISESfinder, etc.. ) and `figures` (raw R plots as PDFs).    

## Workflow  
Different sequencing data sets used in this study are summarized in the following figure:  
<!--![draft workflow](figures/draft_workflow.tiff)-->

## *Pyrococcus assemblies*   
Script and figures for the comparison of different *Pyrococcus furiosus* assemblies have a `dnaseq_`.. prefix.  

## Nanopore sequencing  
This repository contains all necessary code for analyzing the data, starting with the output from <a target="_blank" href = "https://canu.readthedocs.io/en/latest/">`canu`</a> and <a target ="_blank" href = "https://nanopolish.readthedocs.io/en/latest/">`nanopolish`</a>:  
- *.fasta* data are stored in <a target="_blank" href = "https://github.com/felixgrunberger/pyrococcus_annotation/tree/master/data/genome_data">`data/genome_data`</a>  
- Documentation of *scripts* can be found at <a target="_blank" href= "https://github.com/felixgrunberger/pyrococcus_annotation/tree/master/code/nanopore_code">`data/code/nanopore_code`</a>  
- Raw plots and tables are stored in <a target="_blank" href = "https://github.com/felixgrunberger/pyrococcus_annotation/tree/master/figures/nanopore_figures">`data/figures/nanopore_figures`</a>.  

## Illumina d(ifferential) RNA sequencing  
Output from ANNOgesic, R scripts and figures are stored with a `rnaseq_`.. prefix.  



