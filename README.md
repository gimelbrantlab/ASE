# Allelic Imbalance Pipeline - obsolete
This is an OLD sandbox version of the pipeline for allele-specific analysis. 
It has been split into two GitHub repositories:
1. [ASEReadCounter*](https://github.com/gimelbrantlab/ASEReadCounter_star) - going from FASTQ to a table of allelic counts  
and
2. [QCumber](https://github.com/gimelbrantlab/QCumber) - starting with the allelic count table, estimate allelic imbalance and overdispersion

The overall flow of analysis is as follows:
![scheme](https://github.com/gimelbrantlab/ASE/blob/master/ASEReadsCounterstar_QCumber_flowchart.svg)

This pipeline is created to analyze allelic imbalance for F1 crosses of inbred mouse lines. It constructs individual paternal and maternal genomes, then maps the reads from RNA-seq experiments to these genomes and counts the number of reads which map to either the reference or alternate allele at each heterozygous SNP. Next it estimates allelic imbalance for individual genes summarizing information from SNPs and constructs confidence intervals if technical replicates are available.


## File preparation 
You will need to run _prepare_reference.py_ once for every F1 genome and also download a bunch of reference files (for instructions, please see [GenomePreparation.md](https://github.com/gimelbrantlab/ASE/blob/master/markdown/GenomePreparation.md)). 

## Allelic counts
For each RNA-seq experiment, you need to get allelic counts for individual SNPs or genes. Please refer to [Allelic_counts.md](https://github.com/gimelbrantlab/ASE/blob/master/markdown/Allelic_counts.md)) for more information.

## Allelic Imbalance estimation and differential allelic expression

As a result of completing the previous step, you should have a file "*_processed_gene_v3.txt" containing information about number of maternal and paternal counts per gene. Next step is to estimate allelic imbalances for each gene and condition, calculate confidence intervals using technical replicates and then finally detect genes demonstrating differential allelic imbalance. To do this, please refer to the [manual](https://github.com/gimelbrantlab/ASE/blob/master/markdown/manual.md).

## Credits
TODO: Write credits

## License
TODO: Write license


