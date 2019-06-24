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
python3 /home/am717/ASE/python/prepare_reference_tmp.py --PSEUDOREF True --HETVCF True \
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

Make shure that your fasta files are ready for the STAR alignment step, each of them should be [indexed with STAR](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf) (`STAR --runMode genomeGenerate`).

```
STAR --runMode genomeGenerate 
```

```
STAR --readFilesIn $item1 $item2 \
     --outFileNamePrefix $D$item"_on129S1." \
     --runThreadN 4 --outSAMtype SAM \
     --outSAMattrRGline ID:mat \
     --genomeDir /n/scratch2/sv111/ASE/129S1SvImJ_pseudo/ \
     --outFilterMultimapNmax 1 --sjdbGTFfile /n/scratch2/sv111/ASE/Mus_musculus.GRCm38.68.gtf
```

## Allele distributing (merge):

Make shure that your file is sorted by read names:
```
samtools sort -n -O sam -o $D$item"_on129S1.sam" -@ 4 $D$item"_on129S1.Aligned.out.sam"
```
Then merge:
```
python /home/am717/scripts/alleleseq_merge_stream_v2.py \ 
       --pat_sam $D$item"_onCAST.sam" \
       --mat_sam $D$item"_on129S1.sam" \
       --o $D$item"_merged.sam" \
       --paired 1
```
Output: one sam file with mat and pat readgroups per replicate.

## Reads sampling for mutual analysis:

All the sam files in the analysis should be sampled to the same lib size (for example, min(sizes)), see paper for reasoning.

Make shure that your file is sorted by read names:
```
samtools sort -n -O sam -o $samfile -@ 4 $samfile0
```
Then sample (and repeat as many times as you nead, then just process separatelly):

for single end:
```
grep "^@" $samfile > $sampledsam
grep -v "^@" $samfile | shuf -n $(( $minsize )) >> $sampledsam
```

for paired end:
```
grep "^@" $samfile > $sampledsam
grep -v "^@" $samfile | sed '$!N;s/\n/ IHOPETHATNEVERWOULDAPPERINSAMFILE /' | shuf -n $(( $minsize/2 )) | \
      sed 's/ IHOPETHATNEVERWOULDAPPERINSAMFILE /\n/' >> $sampledsam
```

Output: one sampled sam file (x trials) per replicate.

## SNP allele coverage counting:

Convert sam to bam:
```
samtools sort -o $sampledbam $sampledsam
```
Obtain table with SNP allele counts:
```
python /home/am717/scripts/allelecounter.py --vcf $VCF \
       --sample F1 --bam $sampledbam \
       --ref $ref129S1 --min_cov 0 --min_baseq 2 --min_mapq 10 \
       --o $D"input_data/"$sample".sample"$minsize"."$i".stat_0.txt"
```

Output: one table (x trials) per replicate.

## Creating SNP / Grouped SNPs tables:

```
Rscript --vanilla /home/am717/scripts/counts_to_snp_genes.R \ 
        -d $D/input_data \
        -n $prefixes \
        -r $pr_name \
        -o $D \
        -v $snpf1infoexons \
        -b $exonsbed 
        
```

Output: SNP table and Grouped SNP table (for example, genes) (x trials) per set of replicates.


## Correction constant:

```
Rscript /home/am717/scripts/cc_anytab.R -d $D -t $tab -e Geuvadis_$pr_name -n 7
```

Output: RData with Correction Constants.

## CI(AI) and differential analysis:
