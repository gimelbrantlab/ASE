# Files needed for the analysis

1. RNA-seq data for F1-cross, fastq file (single end or paired end)

> For example, paired-end `SRR1106781_1.fastq.gz`,`SRR1106781_2.fastq.gz` and `SRR1106786_1.fastq.gz`,`SRR1106786_2.fastq.gz` replicates. 
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/001/SRR1106781/SRR1106781_1.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/001/SRR1106781/SRR1106781_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/006/SRR1106786/SRR1106786_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/006/SRR1106786/SRR1106786_2.fastq.gz
```

2. Reference genome

> For example, `GRCm38_68.fa`.
```
wget ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa
```

3. gtf annotation for the reference genome

> For example, `Mus_musculus.GRCm38.68.gtf.gz` for corresponding reference genome version.
```
wget ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz
```

4. Either:
* One/Two vcf files for maternal and paternal imbred line (should be "compatible" with reference genome)
* Joint vcf file for multiple species, where two lines are presented (should be "compatible" with reference genome)
* Individual vcf file or joint individuals vcf file (should be "compatible" with reference genome)
(see correponding section about reference preparation)
    
> For example, `mgp.v5.merged.snps_all.dbSNP142.vcf.gz` for multiple mice lines.
```
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
```

> *Note: remember about vcf calling by mutect, varscan, etc, if you heve no pre-existing vcf*


# Input preprocessing before anything:

1. Reference genome fasta should be indexed (`samtools faidx`).

2. vcf-files should be compressed (`bgzip`) and indexed (`samtools tabix`).


# Reference preparation:

One script to rule them all:

```
python3 /full/path/to/ASE/python/prepare_reference_tmp.py --PSEUDOREF True --HETVCF True \
  --pseudoref_dir /full/path/to/dir/for/pseudo/ref/out/ \
  --vcf_dir /full/path/todir/for/vcf/outputs/ \
  --ref /full/path/to/GRCm38_68.fa \
  --name_mat 129S1_SvImJ --name_pat CAST_EiJ \
  --vcf_joint mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
  --gtf /n/scratch2/sv111/ASE/Mus_musculus.GRCm38.68.gtf
```
For help: 
```
python3 /home/am717/ASE/python/prepare_reference_tmp.py --help
```


## Pseudoreference fasta creation:
> `--PSEUDOREF True`

* Input:

|  | Inbred lines | Inbred lines | Individual | Individual | 
| --- | --- | --- | --- | --- |
|  | Joint lines vcf | Separate line(s) vcf | Joint individuals vcf | Separate individual vcf |
| FASTA Reference genome | --ref | --ref | --ref | --ref |
| VCF Variant file(s)[1]    | --vcf_joint | --vcf_mat, --vcf_pat | --vcf_joind | --vcf_ind |
| Name(s)                | --name_mat, --name_pat | --name_mat, --name_pat | --name_ind | --name_ind |
| FASTA Output directory | --pseudo_dir | --pseudo_dir | --pseudo_dir | --pseudo_dir |
| VCF Output directory   | --vcf_dir | --vcf_dir | --vcf_dir | --vcf_dir |

[1] If one or the alleles in case of inbred lines is reference, then everything should be provided as mat or pat only, consistently.

* Output:
  * Pseudoreference genome fastas with own directories.
  * Support vcf or bed files (if needed).


## Heterozygous(parental) VCF creation:
> `--HETVCF True`

* Input:

|  | Inbred lines | Inbred lines | Individual | Individual | 
| --- | --- | --- | --- | --- |
|  | Joint lines vcf | Separate line(s) vcf | Joint individuals vcf | Separate individual vcf |
| FASTA Reference genome | --ref | --ref | --ref | --ref |
| nothing or GTF or BED Selected regions annotation [2] | --gtf or --bed | --gtf or --bed | --gtf or --bed | --gtf or --bed |
| VCF Variant file(s)[1]    | --vcf_joint | --vcf_mat, --vcf_pat | --vcf_joind | --vcf_ind |
| Name(s)                | --name_mat, --name_pat | --name_mat, --name_pat | --name_ind | --name_ind |
| VCF Output directory   | --vcf_dir | --vcf_dir | --vcf_dir | --vcf_dir |

[1] If one or the alleles in case of inbred lines is reference, then everything should be provided as mat or pat only, consistently.
[2] Bed will be considered as main; if gtf provided, automatically considered regions=exons and groups=genes; bed file should have 4 columns and prepared in advance: contig, start position, end position, group ID (no colnames).

* Output:
  * VCF with heterozygous positions, one allele as reference and the second as alternative.
  * Support vcf files (if needed).


# RNA-seq preparation:

## Alignment (STAR) on parental genomes:

* Make shure that your fasta files are ready for the STAR alignment step, each of them should be [indexed with STAR](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf) (`STAR --runMode genomeGenerate`).
```
pseudoRefDirs=/full/path/to/dir/for/pseudo/ref/out/
STAR --runMode genomeGenerate --genomeDir $pseudoRefDirs/129S1_SvImJ/ --genomeFastaFiles $pseudoDir/129S1_SvImJ/129S1_SvImJ_pseudo.fa
STAR --runMode genomeGenerate --genomeDir $pseudoRefDirs/CAST_EiJ/ --genomeFastaFiles $pseudoDir/CAST_EiJ/CAST_EiJ_pseudo.fa
```
* Make shure that you have unziped (`gunzip`) sample fasta/fastq files.

Run alignment:
```
pseudoRefDirs=/full/path/to/dir/for/pseudo/ref/out/

STAR --readFilesIn /full/path/to/SRR1106781_1.fastq.gz /full/path/to/SRR1106781_2.fastq.gz \
     --outFileNamePrefix /full/path/to/alignment/output/SRR1106781_on129S1. \
     --runThreadN 4 --outSAMtype SAM \
     --outSAMattrRGline ID:mat \
     --genomeDir $pseudoRefDirs/129S1_SvImJ/ \
     --outFilterMultimapNmax 1 --sjdbGTFfile /full/path/to/Mus_musculus.GRCm38.68.gtf
STAR --readFilesIn /full/path/to/SRR1106781_1.fastq.gz /full/path/to/SRR1106781_2.fastq.gz \
     --outFileNamePrefix /full/path/to/alignment/output/SRR1106781_onCAST. \
     --runThreadN 4 --outSAMtype SAM \
     --outSAMattrRGline ID:pat \
     --genomeDir $pseudoRefDirs/CAST_EiJ/ \
     --outFilterMultimapNmax 1 --sjdbGTFfile /full/path/to/Mus_musculus.GRCm38.68.gtf
     
STAR --readFilesIn /full/path/to/SRR1106786_1.fastq.gz /full/path/to/SRR1106786_2.fastq.gz \
     --outFileNamePrefix /full/path/to/alignment/output/SRR1106786_on129S1. \
     --runThreadN 4 --outSAMtype SAM \
     --outSAMattrRGline ID:mat \
     --genomeDir $pseudoRefDirs/129S1_SvImJ/ \
     --outFilterMultimapNmax 1 --sjdbGTFfile /full/path/to/Mus_musculus.GRCm38.68.gtf
STAR --readFilesIn /full/path/to/SRR1106786_1.fastq.gz /full/path/to/SRR1106786_2.fastq.gz \
     --outFileNamePrefix /full/path/to/alignment/output/SRR1106786_onCAST. \
     --runThreadN 4 --outSAMtype SAM \
     --outSAMattrRGline ID:pat \
     --genomeDir $pseudoRefDirs/CAST_EiJ/ \
     --outFilterMultimapNmax 1 --sjdbGTFfile /full/path/to/Mus_musculus.GRCm38.68.gtf
```

* Output: aligned reads sam-files for each replicate and parental genome.

## Allele distributing (merge):

Sorted files by read names (`samtools sort -n `):
```
[TODO: concrete]
samtools sort -n -O sam -o /full/path/to/SRR1106781_on129S1.Nsorted.sam -@ 4 /full/path/to/SRR1106781_on129S1.Aligned.out.sam
amtools sort -n -O sam -o /full/path/to/SRR1106781_onCAST.Nsorted.sam -@ 4 /full/path/to/SRR1106781_on129S1.Aligned.out.sam
amtools sort -n -O sam -o /full/path/to/SRR1106786_on129S1.Nsorted.sam -@ 4 /full/path/to/SRR1106786_on129S1.Aligned.out.sam
amtools sort -n -O sam -o /full/path/to/SRR1106786_onCAST.Nsorted.sam -@ 4 /full/path/to/SRR1106786_onCAST.Aligned.out.sam
```
Then merge:
```
python /full/path/to/ASE/python/alleleseq_merge_stream_v2.py \ 
       --mat_sam /full/path/to/SRR1106781_on129S1.Nsorted.sam \
       --pat_sam /full/path/to/SRR1106781_onCAST.Nsorted.sam \
       --o /full/path/to/SRR1106781_merged.sam \
       --paired 1
python /full/path/to/ASE/python/alleleseq_merge_stream_v2.py \ 
       --mat_sam /full/path/to/SRR1106786_on129S1.Nsorted.sam \
       --pat_sam /full/path/to/SRR1106786_onCAST.Nsorted.sam \
       --o /full/path/to/SRR1106786_merged.sam \
       --paired 1
```
Output: one sam file with mat and pat readgroups per replicate.

## Reads sampling for mutual analysis:

All the sam files in the analysis should be sampled to the same lib size (for example, min(sizes)), see paper for reasoning.

Sort merged files by read names (`samtools sort -n `):
```
samtools sort -n -O sam -o /full/path/to/SRR1106781_merged.Nsorted.sam -@ 4 /full/path/to/SRR1106781_merged.sam
samtools sort -n -O sam -o /full/path/to/SRR1106786_merged.Nsorted.sam -@ 4 /full/path/to/SRR1106786_merged.sam
```
Then sample (and repeat as many times as you nead, then just process separatelly), for paired end (the case of the example):
* first calculate sizes (`samtools view -c`):
```
for sam in /full/path/to/SRR1106781_merged.Nsorted.sam /full/path/to/SRR1106786_merged.Nsorted.sam
do
  echo $sam'\t'`samtools view -c $sam` >> /path/to/samsizes.tsv
done
```
* take minimum: 
```
minsize=$(cut -f2 /path/to/samsizes.tsv | sort -V | head -1)
```
* and sample all files to that number of reads, in paired-end case, for example:
```
for sam in /full/path/to/SRR1106781_merged.Nsorted.sam /full/path/to/SRR1106786_merged.Nsorted.sam
do
  grep "^@" $sam > $sam".sample"$minsize"reads.sam"
  grep -v "^@" $sam | sed '$!N;s/\n/ IHOPETHATNEVERWOULDAPPERINSAMFILE /' | shuf -n $(( $minsize/2 )) | \
       sed 's/ IHOPETHATNEVERWOULDAPPERINSAMFILE /\n/' >> $sam".sample"$(($minsize/2))"Preads.sam"
done

```
(for single end, even simplier: pipe of `grep -v "^@"` and `shuf -n $minsize`)

Output: one sampled sam file per replicate.

## SNP allele coverage counting:

Convert sam to sorted bam (`samtools sort`):
```
samtools sort -o /full/path/to/SRR1106781_merged_sample26302221Preads.sorted.bam /full/path/to/SRR1106781_merged.Nsorted.sam.sample26302221Preads.sam
samtools sort -o /full/path/to/SRR1106786_merged_sample26302221Preads.sorted.bam /full/path/to/SRR1106786_merged.Nsorted.sam.sample26302221Preads.sam
```

Obtain table with SNP allele counts:
```
python /home/am717/scripts/allelecounter.py --vcf /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.exons.vcf.gz \
       --bam /full/path/to/SRR1106781_merged_sample26302221Preads.sorted.bam \
       --ref $pseudoDir/129S1_SvImJ/129S1_SvImJ_pseudo.fa \
       --sample F1 --min_cov 0 --min_baseq 2 --min_mapq 10 \
       --o /full/path/to/SRR1106781_merged_sample26302221Preads.stat_0.txt
python /home/am717/scripts/allelecounter.py --vcf /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.exons.vcf.gz \
       --bam /full/path/to/SRR1106786_merged_sample26302221Preads.sorted.bam\
       --ref $pseudoDir/129S1_SvImJ/129S1_SvImJ_pseudo.fa \
       --sample F1 --min_cov 0 --min_baseq 2 --min_mapq 10 \
       --o /full/path/to/SRR1106786_merged_sample26302221Preads.stat_0.txt
```

Output: one table per replicate.

## Creating SNP / Grouped SNPs tables:

Make sure that you have: 
* bed file (with four columns: contig, start and end positions, group ID; no column names) with selected regions, for example, exon regions grouped by genes:
```
awk '$3=="exon" && ($1 ~ /^[1-9]/ || $1 == "X" || $1 == "Y")' /full/path/to/Mus_musculus.GRCm38.68.gtf | cut -f1,4,5,9 | awk -F $'\t' 'BEGIN {OFS = FS} {split($4, a, "gene_id \""); split(a[2], b, "\""); print $1, $2-1, $3, b[1]}' > /full/path/to/output/Mus_musculus.GRCm38.68.EXONS.bed
```
* snp table, which can be obtained from vcf with heterozigous positions created above, as 5 first columns, for example:
```
grep "^#CHROM" /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.exons.vcf | cut -f1-5 > /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.snp_table.txt
grep -v "^#" /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.exons.vcf | sort -V | cut -f1-5 > /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.snp_table.txt
```

```
Rscript --vanilla /home/am717/scripts/counts_to_snp_genes.R \ 
        -d /full/path/to/dir/with/stat/files/ \
        -n SRR1106781_merged_sample26302221Preads,SRR1106786_merged_sample26302221Preads \
        -r Gendrel_81_85 \
        -o /full/path/to/dir/with/stat/files/ \
        -v /full/path/to/Het_Allelic_129S1_SvImJ_CAST_EiJ.snp_table.txt \
        -b /full/path/to/output/Mus_musculus.GRCm38.68.EXONS.bed 
```

Output: SNP table and Grouped SNP table (for example, genes) per set of replicates.


## Correction constant, CI(AI), and differential analysis :
(see manual)

