# Files needed for the analysis

1. RNA-seq data for F1-cross, fastq file (single end or paired end)

> For example, `rep_i.fastq` for each of i={1..N} replicates. 

2. Reference genome

> For example, `GRCm38_68.fa`.
```
wget ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa
```

3. gtf annotation for the reference genome

> For example, `Mus_musculus.GRCm38.68.gtf.gz` for corresponding reference genome version.
```
wget ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz
gunzip Mus_musculus.GRCm38.68.gtf.gz
```

4. Two vcf files for maternal and paternal imbred line (should be "compatible" with reference genome)
  
> For example, `129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz` and `CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz` for 129S1 and CAST mice lines.
```
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz
```
or

Joint vcf file for multiple species, where two lines are presented (should be "compatible" with reference genome)
  
> For example, `mgp.v5.merged.snps_all.dbSNP142.vcf.gz` for multiple mice lines.
```
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
```

> *Note: remember about vcf calling by mutect, varscan, etc*

# Input preprocessing:

1. Reference fasta should be indexed.

2. vcf-files should be bgzip-ed and indexed with tabix.

So the full list of files is: 

# Reference preparation:

```
prepare_reference.sh --ref GRCm38_68.fa --gtf Mus_musculus.GRCm38.68.gtf.gz \
  --name_mat 129S1_SvImJ --name_pat CAST_EiJ --vcf_mat 129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz --vcf_pat CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz \
  --o_dir *output_directory_path/*
prepare_reference.sh --ref GRCm38_68.fa --gtf Mus_musculus.GRCm38.68.gtf.gz \
  --name_mat 129S1_SvImJ --name_pat CAST_EiJ --vcf mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
  --o_dir *output_directory_path/*
```

> *Note: If directories are not writable, will create all supporting files at `output_directory_path/refereces/`.
If something is not provided, it will be also created at `output_directory_path/refereces/` (i.e. will eat some extra place!).*

* Creating pseudogenomes from reference genomes and vcf file(s)

* F1-cross vcf files

                                 
                      joint.vcf =======================> separate.vcf (x2)
                         ||      (GATK: SelectVariants)      //   ||
                         ||                ==================     ||
                  (GATK: ||              // (bcftools: merge)     || (GATK:
          SelectVariant) ||    ==========                         || SelectVariant)
                         \/  \/                                   \/
                       pair.vcf                         separate.snp.vcf (x2)
                  (GATK: ||
          SelectVariant) ||
                         \/
                     pair.snp.vcf
                         ||
         (pair_to_f1.py) ||
                         \/
                     F1.snp.vcf
                                                                                          

* Gene-Transcript-Exon Annotations
