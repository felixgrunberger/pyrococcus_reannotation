# Next generation DNA-Seq and differential RNA-Seq allow re-annotation of the *Pyrococcus furiosus* DSM 3638 genome and provide insights into archaeal antisense transcription  

#### Authors:  
Felix Grünberger, Robert Reichelt, Boyke Bunk, Cathrin Spröer, Jörg Overmann, Reinhard Rachel, Dina Grohmann & Winfried Hausner  

#### Abstract:  
*Pyrococcus furiosus* DSM 3638 is a model organism for hyperthermophilic archaea with an optimal growth temperature near 100∘C. The genome was sequenced about 18 years ago. However, some publications suggest that in contrast to other Pyrococcus species, the genome of *P. furiosus* DSM 3638 is prone to genomic rearrangements. Therefore, we re-sequenced the genome using third generation sequencing techniques. The new de novo assembled genome is 1,889,914 bp in size and exhibits high sequence identity to the published sequence. However, two major deviations were detected: (1) The genome is 18,342 bp smaller than the NCBI reference genome due to a recently described deletion. (2) The region between PF0349 and PF0388 is inverted most likely due an assembly problem for the original sequence. In addition, numerous minor variations, ranging from single nucleotide exchanges, deletions or insertions were identified. The total number of insertion sequence (IS) elements is also reduced from 30 to 24 in the new sequence. Re-sequencing of a 2-year-old “lab culture” using Nanopore sequencing confirmed the overall stability of the *P. furiosus* DSM 3638 genome even under normal lab conditions without taking any special care. To improve genome annotation, the updated DNA sequence was combined with an RNA sequencing approach. Here, RNAs from eight different growth conditions were pooled to increase the number of detected transcripts. Furthermore, a differential RNA-Seq approach was employed for the identification of transcription start sites (TSSs). In total, 2515 TSSs were detected and classified into 834 primary (pTSS), 797 antisense (aTSS), 739 internal and 145 secondary TSSs. Our analysis of the upstream regions revealed a well conserved archaeal promoter structure. Interrogation of the distances between pTSSs and aTSSs revealed a significant number of antisense transcripts, which are a result of bidirectional transcription from the same TATA box. This mechanism of antisense transcript production could be further confirmed by in vitro transcription experiments. We assume that bidirectional transcription gives rise to non-functional antisense RNAs and that this is a widespread phenomenon in archaea due to the architecture of the TATA element and the symmetric structure of the TATA-binding protein.

## Information about this repository  
This Github repository contains code, figures and documentation for reproducing our work and may be used for future similar projects. This work is published in <a href= "https://www.frontiersin.org/articles/10.3389/fmicb.2019.01603/full">Frontiers is Microbiolgy</a>.  

## Structure  
The repository is structured in the folders `code` (R code and Rmd documents), `data` (genomic data and output from ANNOgesic/ISESfinder, etc.. ) and `figures` (raw R plots as PDFs).    


## *Pyrococcus assemblies*   
Script and figures for the comparison of different *Pyrococcus furiosus* assemblies have a `dnaseq_`.. prefix.  

## Nanopore sequencing  
This repository contains all necessary code for analyzing the data, starting with the output from <a target="_blank" href = "https://canu.readthedocs.io/en/latest/">`canu`</a> and <a target ="_blank" href = "https://nanopolish.readthedocs.io/en/latest/">`nanopolish`</a>:  
- *.fasta* data are stored in <a target="_blank" href = "https://github.com/felixgrunberger/pyrococcus_annotation/tree/master/data/genome_data">`data/genome_data`</a>  
- Documentation of *scripts* can be found at <a target="_blank" href= "https://github.com/felixgrunberger/pyrococcus_annotation/tree/master/code/nanopore_code">`data/code/nanopore_code`</a>  
- Raw plots and tables are stored in <a target="_blank" href = "https://github.com/felixgrunberger/pyrococcus_annotation/tree/master/figures/nanopore_figures">`data/figures/nanopore_figures`</a>.  

## Illumina d(ifferential) RNA sequencing  
Output from ANNOgesic, R scripts and figures are stored with a `rnaseq_`.. prefix.  


This project is under the general MIT License.
