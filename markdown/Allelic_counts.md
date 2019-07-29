## How to get Allelic counts per SNP or per gene starting with Fastq files

### Step 1: Download a test RNA-seq data or use your own data.

We will be using RNA-seq data from [Gendrel et al. 2016](https://www.ncbi.nlm.nih.gov/pubmed/27101886). They had two technical replicates for each of 2 NPC clones: 

`SRR1106781_1.fastq.gz` and `SRR1106781_2.fastq.gz` for clone 1 and `SRR1106786_1.fastq.gz` and `SRR1106786_2.fastq.gz` for clone 2.

Download the data:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/001/SRR1106781/SRR1106781_1.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/001/SRR1106781/SRR1106781_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/006/SRR1106786/SRR1106786_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/006/SRR1106786/SRR1106786_2.fastq.gz
```
### Step 2: Run alignment:

You need to align your reads to so-called preudogenomes, i.e. reference genome with SNP from corresponding maternal and paternal genomes. So if you don't have pseudogenomes ready yet, please refer to [instructions](https://github.com/gimelbrantlab/ASE/blob/master/GenomePreparation.md)) to generate them.

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

### Step 3: Merge two files with reads aligned to maternal and to paternal genomes

Sorted files by read names (`samtools sort -n `):
```
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

### Step 4: Reads sampling

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
  echo -e $sam'\t'`samtools view -c $sam` >> /path/to/samsizes.tsv
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
  grep "^@" $sam > $sam".sample"$(( minsize/2 ))"reads.sam"
  grep -v "^@" $sam | sed '$!N;s/\n/ IHOPETHATNEVERWOULDAPPERINSAMFILE /' | shuf -n $(( $minsize/2 )) | \
       sed 's/ IHOPETHATNEVERWOULDAPPERINSAMFILE /\n/' >> $sam".sample"$(( $minsize/2 ))"Preads.sam"
done
```

(for single end, even simplier: pipe of `grep -v "^@"` and `shuf -n $minsize`)

Output: one sampled sam file per replicate.

### Step 5: Extract SNP coverage information from the alignements

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

### Step 6: Get allelic counts per SNP and per gene

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


As a result of the last step, you will get several files in your output folder, including Gendrel_81_85_processed_gene_extended2.txt. This file will be used for the downstream analysis, see description [here](https://github.com/gimelbrantlab/ASE/blob/master/markdown/manual.md).
