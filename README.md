# Allelic Imbalance Pipeline

This pipeline is created to analyze allelic imbalance for F1 crosses of inbred mouse lines. It constructs individual paternal and maternal genomes, then maps the reads from RNA-seq experiments to these genomes and counts the number of reads which map to either the reference or alternate allele at each heterozygous SNP. Next it estimates allelic imbalance for individual genes summarizing information from SNPs and constructs confidence intervals if technical replicates are available.

## Files preparation

You will need to run genome_and_files_preparation.sh once for every F1 genome (for more information, please see [instructions](https://github.com/gimelbrantlab/ASE/blob/master/GenomePreparation.md)). Next, run getASE.sh for all RNA-seq experiments to get allelic imbalance for individual SNPs. The method for constructing confidence intervals is in preparation and coming here soon.

## Credits
TODO: Write credits

## License
TODO: Write license


