export WD=/home/ubuntu
export NTHREADS=32
export picard=$WD/picard/build/libs/picard-2.23.1-2-g3f22b39-SNAPSHOT-all.jar
export gatk=$WD/gatk/gatk
export snpeff=$WD/snpEff/snpEff.jar
export freebayes=$WD/freebayes/build/freebayes

# GATK pipeline
# Call SNPs and indels (Haplotype Caller)
#for folder in $(find $WD/lab_dlbcl/hSNVFiles -name "*_RNA-Seq_gatk")
#do
#   curr=$(basename "$folder" _RNA-Seq_gatk)
#   cd $WD/lab_dlbcl/hSNVFiles/"$curr"_RNA-Seq_gatk
#   echo $curr
#   $gatk HaplotypeCaller -R $WD/GRCh38.p13.genome.fa  -I "$curr"_cutadapt_RG_MD_SC_BQ.bam \
#   -O "$curr"_snp.g.vcf -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation --native-pair-hmm-threads 32 #Call SNP/$
#   cd $WD
#   printf ""$curr"_snp\t$WD/lab_dlbcl/hSNVFiles/"$curr"_RNA-Seq_gatk/"$curr"_snp.g.vcf.gz\n" >> cohort.sample_map.txt # create sample map file
#done

# Create intervals list file
#awk '{print $1}' gencode.v34.annotation.gtf | uniq -d > intervals.list
#mkdir -p $WD/lab_dlbcl_germ_db

# Import single-sample GVCF to GenomicsDB
#$gatk GenomicsDBImport  --genomicsdb-workspace-path  $WD/lab_dlbcl_germ_db  --batch-size 50  -L intervals.list  --sample-name-map cohort.sample_map.txt --reader-threads $NTHREADS

# Joint Genotyping
#$gatk GenotypeGVCFs -R $WD/GRCh38.p13.genome.fa  -V gendb://$WD/lab_dlbcl_germ_db  -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel.vcf.gz -G StandardAnnotation \
#-G AS_StandardAnnotation -G StandardHCAnnotation

# Variant Filtration
# Build the model
#gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz .
#gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi .
#gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz .
#gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi .
#gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz .
#gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi .
#gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz .
#gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi .
#gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz .
#gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi .
# Variant Filtration for SNPs
cd 
#$gatk VariantRecalibrator -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel.vcf.gz \
#   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
#   --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
#   --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
#   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 All_20180418.vcf.gz \
#   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp.recal \
#   --tranches-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp.tranches \
#   --rscript-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp.plots.R

# Variant Filtration for Indels
#$gatk VariantRecalibrator -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel.vcf.gz \
#   --resource:known_indels,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
#   --resource:omni,known=false,training=true,truth=false,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
#   --resource:1000G,known=false,training=true,truth=false,prior=10.0 Homo_sapiens_assembly38.known_indels.vcf.gz \
#   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 All_20180418.vcf.gz \
#   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode INDEL -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel.recal \
#   --tranches-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel.tranches \
#   --rscript-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel.plots.R
#
## Apply VQSR
#$gatk ApplyVQSR -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel.vcf.gz -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp.vcf.gz \
#--truth-sensitivity-filter-level 99.0 --tranches-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp.tranches --recal-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp.recal \
#-mode SNP
#
#$gatk ApplyVQSR -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel.vcf.gz -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel.vcf.gz \
#--truth-sensitivity-filter-level 99.0 --tranches-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel.tranches --recal-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel.recal \
#-mode INDEL 

# Genotype refinement
#$gatk SelectVariants -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp.vcf.gz -select-type SNP -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snps.vcf.gz
#$gatk SelectVariants -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel.vcf.gz -select-type INDEL -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indels.vcf.gz
#
#$gatk CalculateGenotypePosteriors -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snps.vcf.gz -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_GQ_prior.vcf.gz
#$gatk CalculateGenotypePosteriors -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indels.vcf.gz  -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_GQ_prior.vcf.gz 
#
#
#$gatk VariantFiltration -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_GQ_prior.vcf.gz \
#-O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ.vcf.gz \
#-filter "GQ<20" --filter-name "GQ"
#$gatk VariantFiltration -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_GQ_prior.vcf.gz \
#-O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ.vcf.gz \
#-filter "GQ<20" --filter-name "GQ"
##Side note: delete files which are not needed for downstream analysis
#cd $WD/lab_dlbcl/hSNVFiles
#rm $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel.recal $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel.recal.idx $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp.recal \
#$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp.recal.idx lab_dlbcl_filtered_snp.vcf.gz lab_dlbcl_filtered_snp.vcf.gz.tbi lab_dlbcl_filtered_indel.vcf.gz \
#lab_dlbcl_filtered_snp.vcf.gz.tbi lab_dlbcl_snp_indel.vcf.gz lab_dlbcl_snp_indel.vcf.gz.tbi
#cd

#vcftools --gzvcf $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ.vcf.gz --remove-filtered-all --recode --recode-INFO-all \
#--out $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_pass
#vcftools --gzvcf $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ.vcf.gz --remove-filtered-all --recode --recode-INFO-all \
#--out $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_pass

#java -jar $picard SortVcf I=$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_pass.recode.vcf O=$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_pass_sorted.recode.vcf
#java -jar $picard SortVcf I=$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_pass.recode.vcf O=$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_pass_sorted.recode.vcf

# Repeat GATK pipeline with STAR aligner
# Call SNPs and indels (Haplotype Caller)
#for folder in $(find $WD/lab_dlbcl/hSNVFiles -name "*_RNA-Seq_gatk_star")
#do
#   curr=$(basename "$folder" _RNA-Seq_gatk_star)
#   cd $WD/lab_dlbcl/hSNVFiles/"$curr"_RNA-Seq_gatk_star
#   echo $curr
#   $gatk HaplotypeCaller -R $WD/GRCh38.p13.genome.fa  -I "$curr"_cutadapt_RG_MD_SC_BQ_star.bam \
#   -O "$curr"_snp_star.g.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation --native-pair-hmm-threads 32 #Call SNP/$
#   cd $WD
#   printf ""$curr"_snp\t$WD/lab_dlbcl/hSNVFiles/"$curr"_RNA-Seq_gatk_star/"$curr"_snp_star.g.vcf.gz\n" >> cohort.sample_map_star.txt # create sample map file
#done

# Import single-sample GVCF to GenomicsDB
#$gatk GenomicsDBImport --genomicsdb-workspace-path  $WD/lab_dlbcl_germ_db_star  --batch-size 50  -L intervals.list  --sample-name-map cohort.sample_map_star.txt \
#--reader-threads $NTHREADS

# Joint Genotyping
#$gatk GenotypeGVCFs -R $WD/GRCh38.p13.genome.fa  -V gendb://$WD/lab_dlbcl_germ_db_star  -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel_star.vcf.gz -G StandardAnnotation \
#-G AS_StandardAnnotation -G StandardHCAnnotation
#
#cd 
#$gatk VariantRecalibrator -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel_star.vcf.gz \
#   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
#   --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
#   --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
#   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 All_20180418.vcf.gz \
#   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_star.recal \
#   --tranches-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_star.tranches \
#   --rscript-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_star.plots.R

# Variant Filtration for Indels
#$gatk VariantRecalibrator -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel_star.vcf.gz \
#   --resource:known_indels,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
#   --resource:omni,known=false,training=true,truth=false,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
#   --resource:1000G,known=false,training=true,truth=false,prior=10.0 Homo_sapiens_assembly38.known_indels.vcf.gz \
#   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 All_20180418.vcf.gz \
#   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode INDEL -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel_star.recal \
#   --tranches-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel_star.tranches \
#   --rscript-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel_star.plots.R

# Apply VQSR
#$gatk ApplyVQSR -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel_star.vcf.gz \
#-O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_star.vcf.gz \
#--truth-sensitivity-filter-level 99.0 --tranches-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_star.tranches \
#--recal-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_star.recal \
#-mode SNP
#
#$gatk ApplyVQSR -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel_star.vcf.gz \
#-O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_star.vcf.gz \
#--truth-sensitivity-filter-level 99.0 --tranches-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel_star.tranches \
#--recal-file $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indel_star.recal \
#-mode INDEL 

# Genotype refinement
#$gatk SelectVariants -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_star.vcf.gz -select-type SNP -O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snps_star.vcf.gz
#$gatk SelectVariants -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_star.vcf.gz -select-type INDEL \
#-O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indels_star.vcf.gz
#
#$gatk CalculateGenotypePosteriors -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snps_star.vcf.gz \
#-O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_GQ_prior_star.vcf.gz
#$gatk CalculateGenotypePosteriors -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indels_star.vcf.gz \
#-O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_GQ_prior_star.vcf.gz
#
#
#$gatk VariantFiltration -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_GQ_prior_star.vcf.gz \
#-O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_star.vcf.gz \
#-filter "GQ<20" --filter-name "GQ"
#$gatk VariantFiltration -R $WD/GRCh38.p13.genome.fa -V $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_GQ_prior_star.vcf.gz \
#-O $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_star.vcf.gz \
#-filter "GQ<20" --filter-name "GQ"a

#vcftools --gzvcf $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_star.vcf.gz --remove-filtered-all --recode --recode-INFO-all  \
#--out $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_star_pass
#vcftools --gzvcf $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_star.vcf.gz --remove-filtered-all --recode --recode-INFO-all \
#--out $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_star_pass

#java -jar $picard SortVcf I=$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_star_pass.recode.vcf O=$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_star_pass_sorted.recode.vcf
#java -jar $picard SortVcf I=$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_star_pass.recode.vcf O=$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_star_pass_sorted.recode.vcf


# Call SNPs and indels (Freebayes)
#cd
#for folder in $(find $WD/lab_dlbcl/hSNVFiles -name "*_RNA-Seq_gatk")
#do
#   curr=$(basename "$folder" _RNA-Seq_gatk)
#   echo $folder/"$curr"_cutadapt_RG_MD_SC_BQ.bam >> bam.list
#done
#$freebayes --fasta-reference $WD/GRCh38.p13.genome.fa --bam-list bam.list > $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snp_indel_freebayes.vcf

# Separate SNPs and indels
#cd $WD/lab_dlbcl/hSNVFiles
#vcftools --vcf lab_dlbcl_snp_indel_freebayes.vcf  --keep-only-indels --recode --recode-INFO-all --out lab_dlbcl_indels-only_fb.vcf
#vcftools --vcf lab_dlbcl_snp_indel_freebayes.vcf  --remove-indels  --minQ 30  --max-missing 1  --recode --recode-INFO-all --out lab_dlbcl_snps-only_fb.vcf

#cd $WD/lab_dlbcl/hSNVFiles
#shopt -s extglob
#rm -v !(*.plots.R|*.pdf|*.recal|*.idx|*.tranches|*.recode.vcf|*_nolowGQ.vcf.gz|*_nolowGQ.vcf.gz.tbi|*_nolowGQ_star.vcf.gz|*_nolowGQ_star.vcf.gz.tbi)

# Statistics of the vcf files
#gzip -d $WD/lab_dlbcl/hSNVFiles/*.vcf.gz #unzip files

#for file in $(find $WD/lab_dlbcl/hSNVFiles -name "*.vcf")
#do
#   number=$(grep -v "#" $file | wc -l)
#   printf "$file\t$number\n" >> $WD/lab_dlbcl/hSNVFiles/stats.txt
#done

#Find overlapping germline variants between vcf files
#bedtools intersect -a $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_pass_sorted.recode.vcf -b $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_star_pass_sorted.recode.vcf \
#$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snps-only_fb.vcf.recode.vcf -wa -wb -header -filenames > $WD/lab_dlbcl/hSNVFiles/common_snps.vcf
bedtools intersect -a $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_snp_nolowGQ_star_pass_sorted.recode.vcf -b $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_snps-only_fb.vcf.recode.vcf \
> $WD/lab_dlbcl/hSNVFiles/star_fb_overlap.vcf
grep -v "#" $WD/lab_dlbcl/hSNVFiles/star_fb_overlap.vcf | wc -l
rm $WD/lab_dlbcl/hSNVFiles/star_fb_overlap.vcf

#bedtools intersect -a $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_pass_sorted.recode.vcf -b $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_star_pass_sorted.recode.vcf \
#$WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indels-only_fb.vcf.recode.vcf -wa -wb  -header -filenames > $WD/lab_dlbcl/hSNVFiles/common_indels.vcf
bedtools intersect -a $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_filtered_indel_nolowGQ_star_pass_sorted.recode.vcf -b $WD/lab_dlbcl/hSNVFiles/lab_dlbcl_indels-only_fb.vcf.recode.vcf \
> $WD/lab_dlbcl/hSNVFiles/star_fb_overlap_indel.vcf
grep -v "#" $WD/lab_dlbcl/hSNVFiles/star_fb_overlap_indel.vcf | wc -l
rm $WD/lab_dlbcl/hSNVFiles/star_fb_overlap_indel.vcf


