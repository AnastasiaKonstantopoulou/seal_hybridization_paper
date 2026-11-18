The main dataset used here included whole-genome sequences from 119 grey
seal (*Halichoerus grypus*) individuals.  
Read data from 17 Arctic ringed seal (*Pusa hispida*) individuals from
Greenland were added to serve as an outgroup in the *f*-statistics
estimation.

All analyses were performed in the CSC Computing Environment – IT Center
for Science, Finland.

## Setup

For simplicity, it is assumed that all steps are executed from the
current directory and that all tools are available in the `tools`
directory.

``` bash
work=$(pwd)
tools=$work/bin
rmasker=$tools/RepeatMasker
seqbility=$tools/seqbility
bedtools=$tools/bedtools
busco=$tools/busco
miniprot=$tools/miniprot
bwa_mem2=$tools/bwa-mem2
samtools=$tools/samtools
gatk3=$tools/gatk3
gatk4=$tools/gatk4
bcftools=$tools/bcftools
vcftools=$tools/vcftools
beagle=$tools/beagle
genmulti=$tools/generate_multihetsep.py
msmc2=$tools/msmc2
combinecc=$tools/combineCrossCoal.py
msmc_im=$tools/MSMC_IM.py
plink=$tools/plink
smartpca=$tools/eigensoft/bin/smartpca
convertf=$tools/eigensoft/bin/convertf
admixture=$tools/admixture
king=$tools/king
deeptools=$tools/deeptools
minimap2=$tools/minimap2
nextflow=$tools/nextflow
swiftover=$tools/swiftover
blast=$tools/blast
snpeff=$tools/snpEff*
bigWigToBedGraph=$tools/ucsc-tools/bigWigToBedGraph
```

`seqbility` contains `splitfa`, `gen_raw_mask.pl`, `gen_mask` and
`makeMappabilityMask.py`, `blast` contains `makeblastdb` and `blastp`
and `deeptools` contains `bamCoverage`.

# 1. Positive mask and functional categories

## 1.1 Positive mask

The repetitive elements in the reference genome were identified with
RepearMasker and dog repeat libraries, creating the output
${refname}.fa.out.gff:

``` bash
$rmasker -xsmall -gff --engine rmblast -species "Canis familiaris" ${refname}.fa
```

The positive mask was created with Heng Li’s SNPable, following the
instructions of Stephan Schiffels. The mappable regions were identified
by splitting the reference genome into 75 bp fragments and realigning
them to the reference genome. Based on that, a FASTA file with mappable
regions was created with gen_raw_mask.pl and gen_mask, and a BED file
was created with makeMappabilityMask.py.

``` bash
$seqbility/splitfa ${refname}.fa 75 | split -l 20000000

for i in $(ls x??); do bwa aln -R 1000000 -O 3 -E 3 ${refname}.fa $i > $i.sai; done

for i in $(ls x??); do bwa samse ${refname}.fa $i.sai $i > $i.sam; done

cat x??.sam | $seqbility/gen_raw_maks.pl > rawMask.fa

$seqbility/gen_mask -l 35 -r 0.5 rawMask.fa > mask.fa

python makeMappabilityMask.py

for i in $(ls posmask); do zcat posmask/$i | bgzip -c >> posmask.bed.gz; done
```

The repeat elements (negative mask) and the mappable regions (positive
mask) were combined to a joint mask:

``` bash
$bedtools subtract -a posmask.bed -b ${refname}.fa.out.gff > joint_mask.bed
```

## 1.2 Functional categories

**BUSCO analysis**

The BUSCO analysis was performed using the lineage “carnivora”.

``` bash
$busco -i ${refname}.fa -m genome -l carnivora 
```

Gene models in the ribbon seal reference were inferred with miniprot
using the elephant seal protein sequences.

``` bash
$miniprot --gff ${refname}.fa elephant_seal_protein.faa > elephant_seal_genes_tmp.gff
```

BUSCO generated an additional set of gene models, which were extracted
from the BUSCO output.

``` bash
cp BUSCO_${refname}.fa/logs/miniprot_align_carnivora_odb12_out.log > busco_genes_tmp.gff
```

GFF files were sorted and filtered by keeping only gene models with 80%+
identity.

``` bash
grep -v "##" ${gffset}_tmp.gff | sort -k1,1 -k4,5n > ${gffset}.gff

grep mRNA ${gffset}.gff | cut -f9 | sed 's/[;=]/ /g' | awk '$6>=0.8{print $2}' > ${gffset}_keep.txt

grep -f ${gffset}_keep.txt ${gffset}.gff > ${gffset}_keep.gff
```

Matches to BUSCO CDS were identified and extracted.

``` bash
grep CDS elephant_seal_genes_keep.gff | gff2bed | cut -f-3 | bedtools sort -i - | uniq \
> elephant_seal_keep_CDS.bed

grep CDS busco_genes_keep.gff | gff2bed | cut -f-3 | bedtools sort -i - > busco_genes_keep_CDS.bed


$bedtools intersect -a <(grep CDS elephant_seal_genes_keep.gff) -b busco_genes_keep_CDS.bed \
-f 0.9 -F 0.9 | gff2bed | bedtools merge -i - -c 10 -o first > first.bed

cut -f4 first.bed | sed 's/[;=]/ /g' | cut -d" " -f2 | uniq | sort | uniq > keepIDs.txt

grep -f keepIDs.txt elephant_seal_genes_keep.gff > elephant_seal_BUSCO_genes.gff
```

To obtain corresponding Ensembl gene IDs and gene symbols for downstream
analyses, the elephant seal proteins were aligned against ferret protein
sequences (MusPutFur1.0.pep.all; Ensembl version 115) using BLAST.

``` bash
# Build a database of ferret protein sequences

$blast/makeblastdb -in Mustela_putorius_furo.MusPutFur1.0.pep.all.fa -dbtype prot \
-parse_seqids -out ferret_protein


# Blast and keep first match with at least 75% identity. 

$blast/blastp -db ferret_protein -query elephant_seal_protein.faa \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -evalue 1e-7 -num_threads 8 -out blast_out.tsv
  
less blast_out.tsv | awk '{if($3>75 && done[$1]!=1){print $0;done[$1]=1}}' \
> blast_out_75.tsv
```

# 2. WGS data processing

## 2.1 Mapping and variant calling

These steps use BWA-MEM2 for mapping the data in FASTQ format to the
reference genome, SAMtools for processing the output files and GATK for
realignment around indels and variant calling.

``` bash
$bwa_mem2 mem -R${readgroup} -K 100000000 -Y ${refname}.fa ${fastq1}${fastq2}
| $samtools view -h - | $samtools fixmate -m -  ${sample_id}_mapped.bam 

$samtools sort ${sample_id}_mapped.bam -Obam -o ${sample_id}_aligned.bam

$gatk3 -T RealignerTargetCreator \
-R ${refname}.fa \
-I ${sample_id}_aligned.bam \
-o targets.intervals
    
$gatk3 -T IndelRealigner \
-R ${refname}.fa \
-I ${sample_id}_aligned.bam \
-targetIntervals targets.intervals \
-o ${sample_id}_realigned.bam

$samtools markdup ${sample_id}_realigned.bam ${sample_id}_markdup.bam

$gatk4 HaplotypeCaller \
-R ${refname}.fa \
-I ${sample_id}_markdup.bam \
-O ${sample_id}.gvcf.gz \
-ERC GVCF
```

## 2.2 Joint calling

These steps use GATK.

**A map file listing GVCFs for *all* samples was created**

``` bash
for gvcf in *.gvcf.gz; do
sample_id=$(basename "$gvcf" .gvcf.gz)
echo -e "$sample_id\t$(realpath "$gvcf")"
done > samples.map
```

**Single-sample GVCFs were imported into GenomicsDB and jointly
genotyped**

In practice, joint genotyping was done in batches for computational
efficiency, by breaking up the genome in smaller intervals and
specifying `-L $batch` for each.

``` bash
$gatk4 GenomicsDBImport \
    --genomicsdb-workspace-path database \
    --sample-name-map samples.map

$gatk4 GenotypeGVCFs \
    -R ${refname}.fa \
    -V gendb://database \
    -O grey_full.vcf.gz
```

**In addition, the combined Greenland ringed seal GVCFs were jointly
genotyped, including only sites present in the filtered grey seal VCF**

Note: the filtering process is described below.

``` bash
$gatk4 GenotypeGVCFs \
-R ${refname}.fa \
-V ringed.gvcf.gz \
-L grey_final.vcf.gz \
-O ringed.vcf.gz \
--force-output-intervals
```

**And merged with the grey seal VCF**

``` bash
$bcftools merge grey_final.vcf.gz ringed.vcf.gz -Oz -o merged_grey_ringed.vcf.gz
```

## 2.3 Filtering

**Sites fixed for the same ALT allele across samples were removed**

The resulting VCF was used in the SnpEff analysis.

``` bash
$bcftools view -q 0.0001:minor grey_full.vcf.gz -Oz -o grey_seg.vcf.gz
```

**Only regions included in the positive mask were retained**

For this, the positive mask described above was used.

``` bash
$bcftools view -T joint_mask.bed grey_seg.vcf.gz -Oz -o grey_seg_positive.vcf.gz
```

**Biallelic SNPs were kept**

``` bash
$bcftools view -v snps -m2 -M2 grey_seg_positive.vcf.gz \
-Oz -o grey_seg_positive_snp.vcf.gz
```

**Filtering based on sequencing depth**

The maximum and minimum values were set based on the the total mapping
coverage (across all samples) distribution.

``` bash
$bcftools view -i 'INFO/DP>1000 && INFO/DP<3500' \
-Oz -o grey_final.vcf.gz grey_seg_positive_snp.vcf.gz
```

A version of the VCF with stricter filtering was also created to use in
the Identity-by-state (IBS) distance estimation.

``` bash
$bcftools view -i 'INFO/DP>1400 && INFO/DP<3500' \
-Oz -o grey_final_strict.vcf.gz grey_seg_positive_snp.vcf.gz
```

## 2.4 Sequencing depth estimation in 50bp windows

Sequencing depth was estimated along the largest contig in 50bp windows.

``` bash
cat sample_dir.txt | while read smp dir; do 
    $samtools depth -a -r Hfa001 $dir/${smp}_markdup.bam \
    | awk '{s+=$3;if($2%50==0){OFS="\t";print $1,$2,s/50;s=0}}' > $smp.tsv; 
done
```

## 2.5 Assessment of potential asymmetric reference-mapping bias

Quality metrics were compared between the grey seal subspecies and the
ringed seals used as an outgroup for the *f*-statistics estimation.

**Mean depth and missingness per site, for each group**

For this, VCFtools was used.

``` bash
vcf=merged_grey_ringed.vcf.gz

atl=Atlantic_no_HG517.tsv
bal=Baltic_no_HG301.tsv
rin=ringed_smp.tsv

for smp in $atl $bal $rin; do
    base=$(basename $smp | sed 's/_.*\.tsv$//')
    out=stats/vcf_stats_${base}
    $vcftools --gzvcf $vcf --keep $smp --site-mean-depth --out $out
    $vcftools --gzvcf $vcf --keep $smp --missing-site --out $out
done
```

**Mapping quality**

Mapping quality (MQ) was extracted from BAM files with SAMtools. For
each individual, both mean MQ across sites, and the percentage of sites
with MQ\<30 were estimated.

``` bash
$samtools stats $bam | grep ^MAPQ | awk -v smp="$smp" 'BEGIN {OFS="\t"}
{sum += $2 * $3; total += $3; if ($2 < 30) low += $3}
END {print smp, (low/total)*100, sum/total}'> MQ_below_30_mean.tsv
```

**Reference-allele balance**

Reference-allele balance for each individual was calculated as the
proportion of REF reads/total reads at heterozygous sites.

``` bash
$bcftools query -f '%CHROM\t%POS[\t%SAMPLE\t%GT\t%AD]\n' merged_grey_ringed.vcf.gz \
| awk 'BEGIN {FS=OFS="\t"}
{for (i=3; i<=NF; i+=3) {sample=$(i); gt=$(i+1); split($(i+2), ad, ",")
        if (gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") {ref=ad[1]; alt=0
            for (j=2; j<=length(ad); j++) alt+=ad[j]
            total=ref+alt
            if (total>0) print $1,$2,sample,ref/total}}}' > ref_balance.tsv
```

# 3. Demographic analyses with MSMC2

Grey seal data were phased with beagle, separately for Baltic and
Atlantic grey seals.

``` bash
dir=msmc/ne_phasing   

$beagle nthreads=5 ne=100000 gl=${file}.vcf.gz out=$dir/${file}_imputed.vcf.gz
$beagle nthreads=5 ne=100000 gt=$dir/${file}_imputed.vcf.gz out=$dir/${file}_phased.vcf
```

The two phased VCF files were merged into `grey_phased.vcf.gz` with
`bcftools`.  
Then, the combined data were converted to multihetsep.

``` bash
for chr in {001..499}; do
    CHR=Hfa$chr
    len=`grep -w $CHR $FAI | cut -f2`
    if [ $len -gt 100000 ]; then 
        list=""
        for SMP in $SAMPLES; do
        $bcftools view -s $SMP -r $CHR $VCF -Oz -o ctgs/$SMP.$chr.vcf.gz
        list=$list" "msmc/ctgs/$SMP.$chr.vcf.gz
        done
        $genmulti --mask mask/$CHR.mask.bed.gz $list > multi/chr${chr}.multihetsep.txt \
        2> multi/log.txt
        rm $list
    fi
done
```

For the MSMC2 analyses, 6 samples of each subspecies were used.  

For the demographic inference, MSMC2 within populations was run with
four individuals per analysis.  
e.g., one run for the Atlantic:

``` bash
$msmc2 -p 27*1+1*2+1*3 -I 0,1,2,3,4,5,6,7 -o output4_atl/ABCD multi/chr*.multihetsep.txt \
&> output4_atl/ABCD.runlog
```

MSMC2 within populations was also run for two individuals per analysis,
to combine with the across-population runs.

``` bash
$msmc2 -p 27*1+1*2+1*3 -I 0,1,2,3 -o output2_atl/AB multi/chr*.multihetsep.txt &> \
output2_atl/AB.runlog
$msmc2 -p 27*1+1*2+1*3 -I 12,13,14,15 -o output2_bal/GH multi/chr*.multihetsep.txt &> \
output2_bal/GH.runlog
```

MSMC2 across populations was run with two individuals from each
population.  
e.g.:

``` bash
$msmc2 -p 27*1+1*2+1*3 \
-I 12-0,12-1,12-2,12-3,13-0,13-1,13-2,13-3,14-0,14-1,14-2,14-3,15-0,15-1,15-2,15-3 \
-o output4_bal_atl/GAHB multi/chr*.multihetsep.txt &> output4_bal_atl/GAHB.runlog
```

MSMC2 rCCR analysis output files were combined with the corresponding
two-sample MSMC2 analyses with combineCrossCoal.

``` bash
$combinecc output4_bal_atl/GAHB.final.txt output2_bal/GH.final.txt \
output2_atl/AB.final.txt > output4_bal_atl/BG_GAHB.final.txt
```

MSMC_IM were then run on the output:

``` bash
$msmc_im -beta 1e-8,1e-6 -mu 1.826e-8 -o output4_ccr/BG_GAHB output4_ccr/BG_GAHB.final.txt
```

# 4. Genetic variation and population structure

**Create input data for PCA and admixture analysis**

For the principal component analysis (PCA) and the admixture analysis a
thinned version of the grey seal VCF was used. Here, thinning was done
so that the distance between the SNPs is at least 1kb.

``` bash
$plink --vcf grey_final.vcf.gz --allow-extra-chr --bp-space 1000 \
--recode vcf bgz --out grey_final_thin
```

The thinned VCF data were converted to BED format.

``` bash
$plink --vcf grey_final_thin.vcf.gz --allow-extra-chr --allow-no-sex \
--make-bed --out grey_final_thin
```

The FAM and BIM files were slightly modified.

``` bash
#Rename the first column of the FAM file. The smp_pop.tsv used here contains a 
#sample ID and a corresponding population column (Baltic/Atlantic).

mv grey_final_thin.fam grey_final_thin.fam_

awk 'NR==FNR {pop[$1]=$2; next} ($1 in pop) {$1=pop[$1]} {print}' smp_pop.tsv \
grey_final_thin.fam_ > grey_final_thin.fam


#Make the first column of the BIM file numeric 

mv grey_final_thin.bim grey_final_thin.bim_

awk 'BEGIN{OFS="\t"} {ctg = substr($1,4,6) + 0;print ctg,$2,$3,$4,$5,$6}' \
grey_final_thin.bim_ > grey_final_thin.bim
```

## 4.1 PCA

Data were converted to Eigenstrat format and a PCA was performed.

``` bash
base=grey_final_thin

cat > convert.txt << EOF
genotypename:    $base.bed
snpname:         $base.bim 
indivname:       $base.fam 
outputformat:    EIGENSTRAT
genotypeoutname: $base.geno
snpoutname:      $base.snp
indivoutname:    $base.ind
EOF

$convertf -p convert.txt 


mv grey_final_thin.ind grey_final_thin.ind_
sed 's/:/ /' grey_final_thin.ind_ | awk '{print $2,$3,$1}' > grey_final_thin.ind


cat > smartpca.txt << EOF
genotypename:    $base.geno
snpname:         $base.snp
indivname:       $base.ind
evecoutname:     $base.pca.evec.txt
evaloutname:     $base.pca.eval.txt
usenorm:         NO
EOF

$smartpca -p smartpca.txt  > smartpca.out
```

## 4.2 IBS distance estimation

The IBS distance estimation was done with PLINK and used the VCF file
with stricter depth filtering.

``` bash
$plink --vcf grey_final_strict.vcf.gz --distance 1-ibs square gz \
--allow-extra-chr --out grey_final_strict
```

## 4.3 Admixture analysis

This analysis was done with ADMIXTURE, here for *K*=2-5.

``` bash
base=grey_final_thin

for K in {2..5}; do
$admixture --cv $base.bed $K | tee $base.${K}.log
done
```

## 4.4 Relatedness estimation

Relationship inference was done with KING. Prior to the estimation,
`grey_final.vcf.gz` was converted into BED format, similarly to before.
Family IDs were the same for all individuals.

``` bash
$king -b grey_final.bed --related --prefix grey_rel
```

# 5. Assessment of hybrid status of individual HG301

## 5.1 *f*<sub>4</sub>-statistic and *f*<sub>4</sub>-ratio test

*f*-statistics were computed with the R package ADMIXTOOLS2. Here, the
setups used in Figure 2C are shown, as an example.

Prior to this, the merged grey-ringed VCF was converted to BED format
with PLINK and the files were modified as in 4.

``` r
library(tidyverse)
library(admixtools)

# Define individuals and populations

codes_smp <- read.table("merged_grey_ringed.fam")

codes_smp <- codes_smp[, c(1,2)]

colnames(codes_smp) <- c("code", "smp")


## Keep HG301 as its own "population", but leave other DEN_ROD as "DEN_ROD"

codes_smp[codes_smp$code == "DEN_ROD" & codes_smp$smp == "HG301", "code"] <- "DEN_ROD_301"


## Pool populations

codes_smp[startsWith(codes_smp$code, "ICE"), "code"] <- "ICE"

codes_smp[startsWith(codes_smp$code, "UK"), "code"] <- "UK"


## SWE_STA + FIN + EST (Baltic populations that did not undergo extinction; OLD_BAL) 

codes_smp[startsWith(codes_smp$code, "FIN") | startsWith(codes_smp$code, "EST"), "code"] <- "OLD_BAL"

codes_smp[codes_smp$code == "SWE_STA", "code"] <- "OLD_BAL"


extract_f2(pref = "merged_grey_ringed", outdir = "f2_grey_ringed",
             inds = codes_smp$smp, pops= codes_smp$code)


# Compute f4 statistics 

test_ind <- "DEN_ROD_301"

## From the Atlantic

f4.res_atl <- f4(f2_blocks, "OLD_BAL", test_pop_ind[1], "UK","GRE_RIN")

## From the Baltic

f4.res_bal <- f4(f2_blocks, "UK", test_pop_ind[1], "OLD_BAL","GRE_RIN")


# Perform f4-ratio test 

## alpha for Atlantic origin. Population order: A, O, X, C, B

f4r.atl <- qpf4ratio(f2_blocks,as.matrix(tibble("ICE","GRE_RIN",test_ind, "OLD_BAL", "UK")))

## alpha for Baltic origin. Population order: A, O, X, C, B

f4r.bal <- qpf4ratio(f2_blocks,as.matrix(tibble("OLD_BAL","GRE_RIN",test_ind, "UK", "SWE_FAL")))
```

## 5.2 Genome-wide heterozygosity compared to parental populations

First, genetic diversity (*π*) per individual per site was estimated
with VCFtools.

``` bash
for smp in $(cut -f1 metadata.tsv); do
    $bcftools view -s $smp grey_final.vcf.gz | $vcftools --vcf - --site-pi --out $smp
done
```

Heterozygosity (*H*) in 250-kbp windows was computed in R following
Olkkonen and Löytynoja (2023).

window_function.R used here is provided separately.

``` r
source("window_function.R") # Estimates coverage in 250-kbp windows, 
# correcting the window width with the positive mask and nan-sites. 

library(ggplot2) 
library(dplyr)


pi_est_files <- list.files(pattern = "^HG.*\\.sites\\.pi$")

pi_est_windows_list <- list()


for (i in seq_along(pi_est_files)) {
  filename <- pi_est_files[i]
  data <- read.table(filename, sep="\t", header=TRUE)
  pi_est_windows_list[[filename]] <- windows(data)
  print(paste("Processed", filename))
}


all_windows <- bind_rows(pi_est_windows_list, .id = "ind")

all_windows$ind <- gsub("\\.sites\\.pi$", "", all_windows$ind)


# Load metadata

metadata <- read.table("metadata_with.tsv")
colnames(metadata) <- c("ind", "species", "location", "country", "group", "code", "sex")


# Assign HG301 as "Hybrid"

metadata$group[metadata$ind == "HG301"] <- "Hybrid"


# Exclude migrant (HG517)

metadata <- metadata %>%
  filter(!(ind == "HG517"))

all_windows <- all_windows %>%
  filter(!(ind == "HG517"))


# Add group column

ind_to_group <- setNames(metadata$group, metadata$ind)

all_windows$group <- ind_to_group[all_windows$ind]


# Categorize in right-closed bins of width 0.0001, and calculate the percentage of windows 
# falling into each bin

bins.plot.df <- function(dat, ind_name, group_name) {
  x <- dat$pi[dat$nsites > 100000]  # filter windows
  x <- x[!is.na(x)]
  breaks <- c(-0.001, seq(0, max(x) + 1e-5, 0.0001), 0.02)
  h <- hist(x, breaks = breaks, plot = FALSE)$counts
  percent <- h / sum(h) * 100
  bin_mid <- (head(breaks, -1) + tail(breaks, -1)) / 2
  data.frame(ind = ind_name, group = group_name, bin = bin_mid, percent = percent)
}


individuals <- unique(all_windows$ind)

H <- lapply(individuals, function(ind_id) {
  dat_ind <- subset(all_windows, ind == ind_id)
  bins.plot.df(dat_ind, ind_name = ind_id, group_name = dat_ind$group[1])
}) %>% bind_rows()


# Estimate mean percentage per bin, for each group 
#[Baltic, Atlantic, or Hybrid (only one individual)]

H_means <- H_long %>%
  group_by(group, bin) %>%         
  summarise(mean_percent = mean(percent), .groups = "drop")  


# Estimate average heterozygosity per group across individuals and windows

group_means <- all_windows %>%
  group_by(group) %>%
  summarise(mean_H = mean(pi), .groups = "drop")
```

## 5.3 Triangle plots

Triangle plots were generated using the R package *triangulaR*.

The hybrid (HG301) and the migrant (HG517) samples were set as “Test”
individuals, and all the Atlantic and Baltic individuals as parental
populations 1 and 2 (P1 and P2), respectively.

``` bash
less metadata.tsv | cut -f1,5 | awk '$2 == "Baltic" {$2 = "P1"}
$2 == "Atlantic" {$2 = "P2"}; $1 == "HG301" || $1 == "HG517" {$2 = "Test"} {print}' > pop_map.tsv
```

First, allele frequencies for P1 and P2 and absolute allele frequency
differences were calculated as in 8.1.

The selection of SNPs used for computing hybrid indices and interclass
heterozygosity was based on different allele-frequency-difference
(*Δ**A**F*) thresholds (0.7, 0.65, 0.60, and 0.55). Then, the VCF file
was filtered to contain only sites that meet each of these thresholds.
e.g., for *Δ**A**F* ≥ 0.55:

``` bash
$bcftools view -R pos_55.txt -Oz -o grey_55.vcf.gz grey_final.vcf.gz
```

The rest of the analysis was done in R. e.g., for SNPs meeting
*Δ**A**F* ≥ 0.55:

``` r
library(vcfR)
library(triangulaR)
library(ggplot2)
library(wesanderson)
library(ggrepel)


# Load pre-computed diagnostic snps 

vcf_55 <- read.vcfR("grey_55.vcf.gz")


# Estimate hybrid index and heterozygosity 

pop_map <- read.table("pop_map.tsv", header = TRUE)

hi.het_55 <- hybridIndex(vcfR = vcf_55, pm = pop_map, p1 = "P1", p2 = "P2")

hi.het_55$pop[hi.het_55$pop == "P1"] <- "Baltic"
hi.het_55$pop[hi.het_55$pop == "P2"] <- "Atlantic"


# Make triangle plot

p_55 <- triangle.plot(hi.het_55, colors = cols) +
  geom_text_repel(
    aes(hybrid.index, heterozygosity, label = id),
    data = subset(hi.het_55, id %in% c("HG301", "HG517")),
    size = 3,
    nudge_y = -0.02,
    min.segment.length = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggtitle("delta = 0.55, 459 SNPs")
```

## 5.4 Heterozygosity at sites differentiated between subspecies

To distinguish a first-generation hybrid more directly from
later-generation classes, we computed the expected proportion of
heterozygous genotypes in HG301 at the 442 autosomal sites with
*Δ**A**F* ≥ 0.55. This was then compared with the observed proportion
(also see 7.3).

``` r
library(dplyr)

freq <- read.table("Bal_Atl_freq_55_auto.tsv", row.names = NULL, header = TRUE)


dat <- diag_freq %>%
  transmute(CHROM, POS,
            p = AF1B,# Baltic frequency of allele 1
            q = AF1A, # Atlantic frequency of allele 1
            d = p - q,
            m = (p + q) / 2)

het_exp <- function(p, q) {
  r <- (p + q) / 2 
  c(F1= mean(p*(1-q) + (1-p)*q),
    F2 = mean(2 * r * (1-r)),  
    BC_B = mean(r*(1-p) + (1-r)*p),  
    BC_A = mean(r*(1-q) + (1-r)*q))
}

tab <- sapply(c(0.55, 0.60, 0.65, 0.70), function(cut) {
  s <- dat[abs(dat$d) >= cut, ]
  c(n = nrow(s), het_exp(s$p, s$q))
})

colnames(tab) <- paste0("dAF>=", c(0.55, 0.60, 0.65, 0.70))
print(round(tab, 3))
```

# 6. Sex linked contigs and sex identification

For detecting the X contigs, three complementary approaches were used:
coverage ratio comparison between males and females (m:f),
heterozygosity comparison and homology to the dog X chromosome and the
Saimaa ringed seal X-linked scaffolds.

## 6.1 Coverage comparison

This approach follows Catalán et al.(2025).

Coverage was estimated in 10-kbp windows with deepTools for 5 samples of
each morphologically determined sex.

``` bash
cat selected_sample_dir.txt | while read smp dir; do 
    $deeptools/bamCoverage -b $dir/${smp}_markdup.bam -bs 10000 --normalizeUsing RPKM --verbose \
    --outFileFormat bedgraph -o ${smp}_10kb_RPKM.bedgraph
done
```

The bedgraphs were combined into a single BED file using
`bedtools unionbedg`.

Then, computations were then done in R.

``` r
library(dplyr)
library(readr)
library(tidyr)

cov <- read.table("selected_samples_10kb_RPKM.bed", header=TRUE)


# Sex info for the selected individuals

sexinfo <- read_tsv("smp.sex.tsv", col_names = c("ind", "sex"))


# After determining the the five largest contigs showing no notable coverage differences 
# between samples, we used their average coverage (across windows of the five contigs) 
# per individual to normalize the coverage of each window for that individual. 

biggest_auto_chrs <- c("Hfa001", "Hfa002", "Hfa003", "Hfa004", "Hfa005")

auto_mean <- cov %>%
  filter(chrom %in% biggest_auto_chrs) %>%
  summarise(across(-c(chrom, start, end), mean, na.rm = TRUE))

cov_norm <- cov
for (s in colnames(cov_norm)[-(1:3)]) {
  cov_norm[[s]] <- cov_norm[[s]] / as.numeric(auto_mean[[s]])
}


# Add sex info and perform a non-parametric Wilcoxon Rank Sum Test to test for significant 
# differences in contig coverage values between sexes, applying a Bonferroni multiple 
# test correction, as in Catalán et al. (2025).


## Add preliminary sex info to the full cov file

cov_norm_sex <- cov_norm %>%
  pivot_longer(-c(chrom, start, end),
               names_to = "ind",
               values_to = "coverage") %>%
  left_join(sexinfo, by = "ind")


## Perform test

results <- cov_norm_sex %>%
  group_by(chrom) %>%
  summarise(
    p = wilcox.test(coverage ~ sex)$p.value),
    male_mean = mean(coverage[sex == "m"], na.rm = TRUE),
    female_mean = mean(coverage[sex == "f"], na.rm = TRUE)
  ) %>%
  mutate(p_adj = p.adjust(p, method = "bonferroni"))


results <- results %>%
  mutate(mf_ratio = male_mean / female_mean)


# From the contigs that show significant differences, the ones with a m:f 
# coverage ratio of 0.5±0.1 were selected

res_x_contigs <- results %>%
  filter(p_adj < 0.05 & mf_ratio >= 0.4 & mf_ratio <= 0.6)

# To verify the assigned sex,the coverage values of single individuals (not only
# the mean of m and f) were also inspected.
```

## 6.2 Sex identification

One of the candidate X contigs identified with this approach, Hfa209,
was observed to contain *TSPYL2*, which in dogs and humans is located on
the X chromosome. Males and females in the full dataset were separated
by comparing mean coverage values (averaged across windows) for all
individuals on this contig.

For this, coverage was estimated in 10-kbp windows for contig Hfa209 in
all individuals, as well as for the autosomal contigs Hfa001-Hfa005 (for
normalization). Coverage estimates were generated separately for each
contig and were then combined into a single file (per individual).

``` r
regions="Hfa001 Hfa002 Hfa003 Hfa004 Hfa005 Hfa209"

cat sample_dir.txt | while read smp dir; do
    out=${smp}_10kb_RPKM.bedgraph
    > $out
    for r in $regions; do
        bamCoverage -b ${dir}/${smp}_markdup.bam -bs 10000 --normalizeUsing RPKM \
        --region $r --outFileFormat bedgraph -o tmp_${smp}_${r}.bedgraph
        
        grep -v "^track" tmp_${smp}_${r}.bedgraph >> $out #Remove header and append
        rm tmp_${smp}_${r}.bedgraph # Remove temporary file
    done
done
```

Again, the bedgraphs were combined into a single BED file using
`bedtools unionbedg`.

Sex identification was done in R.

``` r
cov_total <- read.table("all_samples_10kb_RPKM.bed", header=TRUE)


# Normalize, as before

cov_norm_total <- cov_total
for (s in colnames(cov_norm_total)[-(1:3)]) {
  cov_norm_total[[s]] <- cov_norm_total[[s]] / as.numeric(auto_mean[[s]])
}


# Take average across Hfa209 windows

cov_norm_mean_Hfa209 <- cov_norm_total %>%
  filter(chrom== "Hfa209") %>%
  summarise(across(-c(start, end), mean, na.rm = TRUE)) 


cov_norm_mean_Hfa209_long <- cov_norm_mean_Hfa209 %>%
  pivot_longer(-chrom, names_to = "ind", values_to = "coverage") 

# This can be plotted to separate males and females.
```

## 6.3 Heterozygosity comparison

Heterozygosity (*H*) in 250-kbp windows was estimated as in R 5.2, and
contigs with difference in *H* between the previously identified males
and females were detected.

``` r
# Heterozygosity in windows for each individual was normalized by the mean H of the 
# five largest autosomal contigs.

auto_mean <- all_windows %>% #contains windowed H estimates for all individuals; see 5.2
  filter(ctg %in% biggest_auto_chrs) %>%
  group_by(ind) %>%
  summarise(mean_pi_ind = mean(pi, na.rm = TRUE), .groups = "drop")


all_windows_norm <- all_windows %>%
  left_join(auto_mean, by = "ind") %>%
  mutate(pi_norm = pi / mean_pi_ind)


# Add sex info for all individuals

sexinfo_all <- read_tsv("smp_all.sex.tsv", col_names = c("ind", "sex"))

all_windows_norm_sex <- all_windows_norm %>%
  left_join(sexinfo_all, by = "ind")


# Perform Wilcoxon Rank Sum Test, as before

results_H <- all_windows_norm_sex %>%
  group_by(ctg) %>%
  summarise(
    p = wilcox.test(pi_norm ~ sex)$p.value)
    male_mean = mean(pi_norm[sex == "m"], na.rm = TRUE),
    female_mean = mean(pi_norm[sex == "f"], na.rm = TRUE)
  ) %>%
  mutate(p_adj = p.adjust(p, method = "bonferroni"))


res_x_contigs_H <- results_H %>%
  filter(p_adj < 0.05)

# The contigs with significant differences between sexes and near zero H in males 
# were compared to the X contigs determined by the coverage approach.
```

## 6.4 Alignment approach

As a third line of evidence, the ribbon seal reference genome was mapped
against a dog (ROS_Cfam_1.0.dna_sm.toplevel; Ensembl version 115) and a
Saimaa ringed seal genome (mPhs449_Pusa_saimensis_1; Grethlein et
al. 2026), using minimap2.

For example, for the Saimaa ringed seal:

``` bash
minimap2 saimaa.fa ${refname}.fa -o ribbon_vs_saimaa.paf

# Add PID column (BLAST-like alignment identity)

awk '{pid=($10/$11)*100; print$0, pid}' ribbon_vs_saimaa.paf > ribbon_vs_saimaa_pid.paf

# Focus on regions mapping to the Saimaa X chromosome (CM183163.1)

awk '$6=="CM183163.1"' ribbon_vs_saimaa_pid.paf > ribbon_vs_saimaa_X.paf

# Inspect based on PID and alignment (number of matching bases) length 
# (e.g., pid>60%, length>10kbp)

column -t ribbon_vs_saimaa_X.paf | awk '$10>10000 && $19>60 {print}' \
> ribbon_vs_saimaa_X_10kbp_pid60.paf
```

Based on the three lines of evidence, a final set of X contigs was
compiled.

**Calculate the total length of X contigs**

``` bash
while read contig; do
    grep -w $contig ${refname}.fa.fai | cut -f1,2
done < X_contigs_final.tsv > X_contigs_final_length.tsv

less X_contigs_final_length.tsv | awk '{sum += $2} END {print sum}' 
```

# 7. Hybrid’s parental populations

To identify the hybrid’s maternal populations, a PCA using only X-linked
contigs and only male samples was performed.

## 7.1 Haploid variant calling of X-contig regions in males

X-contig regions were extracted from BAM files of male samples and
variants were called with the option `ploidy 1`.

``` bash
$samtools view -b -M ${sample_id}_markdup.bam -L X_contigs_final.bed -h -o ${sample_id}_X.bam

$gatk4 HaplotypeCaller \
-R ${refname}.fa \
-I ${sample_id}_X.bam \
-O ${sample_id}_X.gvcf.gz \
-L X_contigs_final.list \
-ploidy 1 \
-ERC GVCF
```

Joint calling was performed as in 2.2, which produced
`X_ctg_male.vcf.gz`.

The VCF was filtered to include only biallelic SNPs that are within
positively masked regions and segregate among the samples, as in 2.3.
Then, depth filtering was done as follows, and was based on the mapping
coverage distribution across 65 males. Additionally, any 10-kbp window
that deviated from the expected 0.5±0.1 m:f coverage ratio was masked.

``` bash
$bcftools view -i 'INFO/DP>100 && INFO/DP<800' \
-Oz -o X_ctg_male_dpt.vcf.gz  X_ctg_male_seg_positive_snp.vcf.gz

$bcftools view -T ^X_mask_regions.bed X_ctg_male_dpt.vcf.gz \
-Oz -o X_ctg_male_final.vcf.gz
```

## 7.2 PCA

Input data were generated and the output FAM and BIM files were modified
as in 4.

``` bash
$plink --vcf X_ctg_male_final.vcf.gz --allow-extra-chr --allow-no-sex --make-bed \
--chr $x1 $x2 ... \
--remove females.txt \
--out males_x
```

The BIM file was further modified so that the contigs artificially start
from number 1.

``` bash
mv males_x.bim males_x.bim_

less males_x.bim_ | awk 'BEGIN{ctg=0;OFS="\t"}{if($1!=pc){ctg+=1}pc=$1;$1=ctg;print $0}' \
> males_x.bim
```

Next, data were converted to Eigenstrat format and a PCA was performed,
as in 4.1.

## 7.3 Baltic-like and Atlantic-like alleles carried by HG301

To estimate allele frequencies on X-linked contigs, the VCF file for
males was first merged with the VCF file for females at their common
sites.

``` bash
# From the full VCF, extract sites included in the VCF of X contigs for males
#keeping only females

$bcftools view -T $Xsites_males -S $female_samples grey_final.vcf.gz -Oz -o X_ctg_females.vcf.gz

# Keep only sites of the male VCF that are included in the female VCF, and merge

$bcftools view -T $Xsites_females X_ctg_male_final.vcf.gz -Oz -o X_ctg_male_female_sites.vcf.gz

# Merge

$bcftools merge X_ctg_females.vcf.gz X_ctg_male_female_sites.vcf.gz -Oz -o X_ctg_females_males.vcf.gz
```

Allele frequencies were calculated as in 8.1, for both the full and the
X-contig VCF. e.g., for the X-contig VCF:

``` bash
$vcftools --gzvcf X_ctg_females_males.vcf.gz --keep Baltic_no_HG301.tsv --freq \
--out Baltic_X
$vcftools --gzvcf X_ctg_females_males.vcf.gz --keep keep Atlantic_no_HG517.tsv --freq \
--out Atlantic_X
```

Allele “frequencies” were also calculated for the hybrid, for both the X
contigs and the full VCF.

``` bash
$vcftools --gzvcf X_ctg_females_males.vcf.gz --keep HG301.txt --freq --out HG301_X
$vcftools --gzvcf grey_final.vcf.gz --keep HG301.txt --freq --out HG301
```

Then, the data were imported in R as in 8.1, and allele frequency
differences were estimated. Alleles at sites with *Δ**A**F* ≥ 0.55 were
classified as either Baltic-like or Atlantic-like, for both autosomal
and X contigs, depending on their frequency in each population. Then, it
was determined which allele the hybrid carries at each site.

``` r
X <- c("Hfa025", "Hfa047", "Hfa049", "Hfa063", "Hfa102", "Hfa139", "Hfa142", 
       "Hfa147", "Hfa160", "Hfa173", "Hfa174", "Hfa177", "Hfa178", "Hfa205", 
       "Hfa209", "Hfa225", "Hfa239", "Hfa241", "Hfa253", "Hfa256", "Hfa260", 
       "Hfa269", "Hfa274", "Hfa326", "Hfa340", "Hfa358", "Hfa363", "Hfa367", 
       "Hfa370", "Hfa388", "Hfa389", "Hfa392", "Hfa395", "Hfa404", "Hfa414", 
       "Hfa427", "Hfa428", "Hfa432", "Hfa440", "Hfa441", "Hfa442", "Hfa447", 
       "Hfa451", "Hfa460", "Hfa467", "Hfa469", "Hfa470", "Hfa488", "Hfa489", 
       "Hfa491") #X contigs 

# Autosomes

Bal_Atl_freq$AF_diff <- Bal_Atl_freq$AF1B - Bal_Atl_freq$AF1A # not absolute difference


## Determine which population has the REF allele in higher frequency at sites 
## with ΔAF >= 0.55

Bal_Atl_freq_55 <- Bal_Atl_freq %>%
  filter(abs(AF_diff) >= 0.55) %>%
  mutate(Higher_AF = if_else(AF_diff < 0, "Atlantic", "Baltic")) %>%
  filter(!(CHROM %in% X)) # Keep only autosomal contigs

## Merge with frequencies for hybrid

Hybrid <- read.table("HG301.frq", row.names = NULL)
colnames(Hybrid) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "AF1", "AF2")

Hybrid$AF1 <- as.numeric(sub(".*:", "", Hybrid$AF1))
Hybrid$AF2 <- as.numeric(sub(".*:", "", Hybrid$AF2))


Bal_Atl_Hybrid_55 <- Bal_Atl_freq_55 %>%
  left_join(Hybrid %>% select(CHROM, POS, AF1), by = c("CHROM", "POS"))

## Determine if the hybrid carries the Baltic-like or Atlantic-like allele in 
## each position. e.g., if the REF allele is more common in the Baltic, and the 
## hybrid is homozygous for it -> Baltic

Bal_Atl_Hybrid_55 <- Bal_Atl_Hybrid_55 %>%
  mutate(Hybrid_allele_pop = case_when(
      AF1 == 0.5 ~ "Heterozygous",
      AF1 == 1 & AF_diff > 0 ~ "Baltic",
      AF1 == 0 & AF_diff < 0 ~ "Baltic",
      AF1 == 1 & AF_diff < 0 ~ "Atlantic",
      AF1 == 0 & AF_diff > 0 ~ "Atlantic"))


# X contigs: similar steps to those for the autosomes, but with no "Heterozygous" 
# category (the hybrid is a male, and therefore haploid for the X chromosome)

Bal_Atl_Hybrid_X_55 <- Bal_Atl_Hybrid_X_55 %>% #This was also done for sites with ΔAF >= 0.50
  mutate(Hybrid_allele_pop = case_when(
      AF1 == 1 & AF_diff > 0 ~ "Baltic",
      AF1 == 0 & AF_diff < 0 ~ "Baltic",
      AF1 == 1 & AF_diff < 0 ~ "Atlantic",
      AF1 == 0 & AF_diff > 0 ~ "Atlantic"))


# Combine the data for the X contigs and the autosomes

HG301_alleles <- bind_rows(Bal_Atl_Hybrid_X_55 %>% mutate(Chromosome = "X"),
  Bal_Atl_Hybrid_55 %>% mutate(Chromosome = "Autosomes"))


# Count Baltic-like and Atlantic like alleles, and estimate percentages for both 
# X-contig and autosomal sites
# For heterozygous sites, 1/2 of the count was assigned to each population

HG301_counts <- HG301_alleles %>%
  group_by(Chromosome) %>%
  summarise(Baltic = sum(Hybrid_allele_pop == "Baltic") 
            + 0.5*sum(Hybrid_allele_pop == "Heterozygous"),
    Atlantic = sum(Hybrid_allele_pop == "Atlantic")
    + 0.5*sum(Hybrid_allele_pop == "Heterozygous"),
    .groups = "drop") %>%
  mutate(Baltic_percent = Baltic/(Baltic + Atlantic)*100,
    Atlantic_percent = Atlantic/(Baltic + Atlantic)*100)
```

# 8. Allele frequency estimation and gene functions

## 8.1 Allele frequency differences (ΔAF)

#### 8.1.1 ΔAF estimation and visualization in chromosomes

  

First, allele frequencies per site were estimated separately for the
Baltic and the Atlantic individuals, using VCFtools.  
The hybrid (HG301) was excluded from the Baltic samples and the migrant
(HG517) from the Atlantic ones.

``` bash
$vcftools --gzvcf grey_final.vcf.gz --keep Baltic_no_HG301.tsv --freq --out Baltic
$vcftools --gzvcf grey_final.vcf.gz --keep Atlantic_no_HG517.tsv --freq --out Atlantic
```

Then, absolute allele frequency differences were computed in R.

``` r
# Read data (e.g., for the Baltic)

Baltic_frq <- read.table("Baltic.frq", row.names = NULL)

colnames(Baltic_frq) <- c("CHROM", "POS", "N_ALLELES", "N_CHRB", "AF1B", "AF2B")
#AF1 corresponds to the REF and AF2 to the ALT allele.

# Make freq numeric

Baltic_frq$AF1B <- as.numeric(sub(".*:", "", Baltic_frq$AF1B))
Baltic_frq$AF2B <- as.numeric(sub(".*:", "", Baltic_frq$AF2B))

# Merge with the data for the Atlantic and compute absolute ΔAF 

Bal_Atl_freq <- merge(Baltic_frq[, -c(3,4)],Atlantic_frq[, -c(3,4)],
                      by=c("CHROM", "POS"), sort=FALSE)

Bal_Atl_freq$AF_diff <- abs(Bal_Atl_freq$AF1B - Bal_Atl_freq$AF1A)


write.table(Bal_Atl_freq, "Bal_Atl_AF_diff.tsv", quote=FALSE, col.names=TRUE,
            row.names = FALSE, sep="\t")
```

**Visualization in chromosomes**

To visualize ΔAF in chromosomes, the values were lifted over to a dog
genome.

For this, a liftover chain between the ribbon seal and the dog genome
was generated using the Nextflow LiftOver pipeline. The pipeline is
automatically fetched from GitHub by `nextflow`, while helper tools had
been installed locally.

``` bash
#Get dog genome

wget https://ftp.ensembl.org/pub/release-115/fasta/canis_lupus_familiaris/dna/\
Canis_lupus_familiaris.ROS_Cfam_1.0.dna_sm.toplevel.fa.gz

#Create seal-to-dog chain file

gunzip Canis_lupus_familiaris.ROS_Cfam_1.0.dna_sm.toplevel.fa.gz

dog=Canis_lupus_familiaris.ROS_Cfam_1.0.dna_sm.toplevel

$nextflow run evotools/nf-LO --source ${refname}.fa --target ${dog}.fa --distance far \
--aligner minimap2 --outdir seal_to_dog --no_netsynt
```

Next, the *Δ**A**F* values were transferred to the corresponding
coordinates of the dog genome using swiftover
(<https://github.com/blachlylab/swiftover>).

``` bash
# Convert the ΔAF TSV file into a BED file 

less Bal_Atl_AF_diff.tsv | awk -F"\t" 'BEGIN {OFS="\t"} \
{print $1, $2-1, $2, $3"&"$4"&"$5"&"$6"&"$7}' > Bal_Atl_AF_diff.bed


# Transfer values

$swiftover -t bed -c $seal_to_dog_chain -i Bal_Atl_AF_diff.bed \ 
-u unmatched_seal_to_dog.bed > dog_AF_diff.bed

mv dog_AF_diff.bed dog_AF_diff.bed_

less dog_AF_diff.bed_ | sed 's/&/\t/g' > dog_AF_diff.bed
```

#### 8.1.2 Bootstrap support of outlier SNPs and peaks

  

The bootstrap support of outlier SNPs (meeting the cutoff of
*Δ**A**F* ≥ 0.55) was evaluated by resampling individuals with
replacement within the Baltic and Atlantic populations, estimating
allele frequencies and *Δ**A**F* values for 1000 replicates, and
assessing how frequently each initially identified outlier SNP met the
*Δ**A**F* cutoff across replicates.

First, a genotype matrix was generated with VCFtools. The VCF file used
did not include the hybrid (HG301) and migrant (HG517) individuals.

``` bash
$vcftools --gzvcf $vcf --012 --out grey_no_hyb_mig
```

The estimation of allele frequencies and *Δ**A**F* values, and the rest
of the analysis was done in R, using the following script.

``` r
library(dplyr)
library(data.table)

gt <- data.frame(t(read.table("grey_no_hyb_mig.012")))
gt <- gt[-1,]
gt <- gt %>% mutate(across(everything(), ~ na_if(., -1)))


pop <- read.table("ind_pop.tsv", header=TRUE) # has 3 columns: sample, site code, 
# and subspecies (population)

baltic <- which(pop$pop == "Baltic")
atlantic <- which(pop$pop == "Atlantic")

n.baltic <- length(baltic)
n.atlantic <- length(atlantic)


# Group individuals in major subpopulation groups, within each subspecies 

strata_map <- c(DEN_ROD="BAL_ROD", DEN_CHR="BAL_CHR",
                SWE_FAL="BAL_NE", SWE_STA="BAL_NE", EST_GOR="BAL_NE",
                FIN_BOB="BAL_NE", FIN_GOF="BAL_NE",
                DEN_THY="ATL_NS", UK_ILM="ATL_UK",
                ICE_BRE="ATL_ICE", ICE_PAT="ATL_ICE",
                ICE_STR="ATL_ICE", ICE_DJU="ATL_ICE",
                NOR_FRO="ATL_NR", RUS_KOL="ATL_NR")


pop$stratum <- strata_map[pop$loc]

bal.strata <- split(baltic, pop$stratum[baltic])
atl.strata <- split(atlantic, pop$stratum[atlantic])


# Function to estimate ΔAF

af <- function(M) {
  n_obs <- rowSums(!is.na(M))
  ifelse(n_obs > 0, rowSums(M, na.rm = TRUE) / (2 * n_obs), NA_real_)
}


# Estimate ΔAF for original sample list 

org.baltic <- af(gt[, baltic, drop = FALSE])
org.atlantic <- af(gt[, atlantic, drop = FALSE])

org.af <- abs(org.baltic - org.atlantic)


# Identify SNPs above threshold

threshold <- 0.55

selected.snps <- which(org.af >= threshold)

s1 <- numeric(nrow(gt))
s2 <- numeric(nrow(gt))


# Bootstrapping

set.seed(777)

n.boot <- 1000

boot.count <- numeric(nrow(gt))

for (i in 1:n.boot) {
  
  # Resample within each Baltic site
  baltic.smp <- unlist(
    lapply(bal.strata, function(site) {
      sample(site, length(site), replace = TRUE) 
    }), use.names = FALSE)
  
  # Resample within each Atlantic site
  atlantic.smp <- unlist(
    lapply(atl.strata, function(site) {
      sample(site, length(site), replace = TRUE)
    }), use.names = FALSE)
  
  # Calculate bootstrap allele frequencies
  boot.baltic <- af(gt[, baltic.smp, drop = FALSE])
  boot.atlantic <- af(gt[, atlantic.smp, drop = FALSE])
  
  # Calculate ΔAF
  boot.af <- abs(boot.baltic - boot.atlantic)
  
  s1 <- s1 + ifelse(is.na(boot.af), 0, boot.af)
  s2 <- s2 + ifelse(is.na(boot.af), 0, boot.af^2)
  
  # Count how often each SNP meets threshold
  boot.count <- boot.count + (!is.na(boot.af) & boot.af >= threshold)
}


# Get mean and SE

boot.mean <- s1/n.boot
boot.se   <- sqrt(pmax(0, (s2 - n.boot*boot.mean^2)/(n.boot - 1)))


# Make dataframe of bootstrapping results

results <- data.frame(
  SNP = selected.snps,
  original_AFdiff = org.af[selected.snps],
  bootstrap_support = 100 *
    boot.count[selected.snps] / n.boot,
  boot_se = boot.se[selected.snps] )


# Add position information

pos <- read.table("grey_no_hyb_mig.012.pos", header=FALSE)

pos.selected <- pos[selected.snps, ]

pos.selected$SNP <- rownames(pos.selected)

colnames(pos.selected) <- c("ctg", "pos", "SNP")

results_pos <- results %>%
  mutate(SNP = as.character(SNP)) %>%
  left_join(
    pos.selected %>% mutate(SNP = as.character(SNP)),
    by = "SNP"
  ) %>% select (
    ctg, pos, original_AFdiff, bootstrap_support, boot_se)


saveRDS(list(count = boot.count, s1 = s1, s2 = s2), "boot_all.rds")


write.table(results, "AF_boot_55.tsv", sep = "\t", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)
write.table(results_pos, "AF_pos_boot_55.tsv", sep = "\t", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)
```

**Group SNPs in peaks and estimate bootstrap statistics for peaks**

Consecutive outlier SNPs within each contig were grouped into the same
ΔAF peak if the distance between them was ≤ 500 kbp.

``` r
library(dplyr)

# Get file with ΔAF values and keep positions with ΔAF >=0.55.

pos_diff <- read.table("Bal_Atl_AF_diff.tsv", row.names = NULL, header=TRUE)

pos_diff_55 <- subset(pos_diff_new, AF_diff>=0.55)

# Group in peaks

pos_diff_55_peaks <- pos_diff_55 %>%
  arrange(CHROM, POS) %>%
  group_by(CHROM) %>%
  mutate(gap = POS - lag(POS, default = dplyr::first(POS)),
    peak = cumsum(gap > 500000)) 


# Get bootstrap results

AF_diff_boot_55 <- read.table("AF_pos_boot_55.tsv")

colnames(AF_diff_boot_55) <- c("CHROM", "POS", "original_AFdiff", 
                               "bootstrap_support", "boot_se")


# Add peak information

AF_diff_boot_55 <- AF_diff_boot_55 %>%
  left_join(pos_diff_55_peaks, by = c("CHROM", "POS")) %>%
  select(CHROM, POS, original_AFdiff, bootstrap_support, boot_se, gap, peak)


# Compute statistics for peaks

## Bootstrap support of SNP with highest ΔAF per peak

AF_diff_boot_55_max <- AF_diff_boot_55 %>%
  group_by(CHROM, POS) %>%
  mutate(n_snps = n()) %>%
  slice_max(original_AFdiff, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(ctg_peak = paste(CHROM, peak, sep = "_")) %>%
  arrange(bootstrap_support)

## Median bootstrap support across SNPs

AF_diff_boot_55_median <- AF_diff_boot_55 %>%
  group_by(CHROM, peak) %>%
  summarise(n_snps = n(), median = median(bootstrap_support), mean = mean(bootstrap_support)) %>%
  mutate(ctg_peak = paste(CHROM, peak, sep = "_"))

## Combine

AF_diff_boot_55_sum <- left_join(AF_diff_boot_55_max, AF_diff_boot_55_median, by = "ctg_peak") 

AF_diff_boot_55_sum <- AF_diff_boot_55_sum %>%
  mutate(peak = ctg_peak, ctg = CHROM.x, pos_max = POS, AF_max = original_AFdiff, 
         boot_max = bootstrap_support, se_max = boot_se, n_snps = n_snps.x) %>% 
  select(peak, ctg, pos, n_snps, median, mean, pos_max, AF_max, boot_max, se_max) %>%
  arrange(desc(median), desc(boot_max))  %>%
  mutate(peak = factor(peak, levels = unique(peak)))
```

## 8.2 Genes and functions

#### 8.2.1 Putative functionally-relevant genes surrounding high-ΔAF peaks

  

To inspect genes in proximity to the areas of high differentiation, we
used the functional region annotation described above.

Genes within 1 Mbp of regions of high allele frequency difference
(*Δ**A**F* ≥ 0.55) were extracted and then filtered based on Gene
Ontology (GO) terms, retaining those with annotations related ion and
osmotic regulation or lipid metabolism. For this, GO annotations for
ferret gene IDs (MusPutFur1.0) were retrieved from Ensembl using the R
package biomaRt.

``` r
library(biomaRt)
library(dplyr)
library(tidyr)


# Select database and dataset

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "mpfuro_gene_ensembl"
)


# Ensembl gene IDs to get GO terms for

genes_t <- read.table("genes_AF_diff_55_1Mbp.tsv", quote = "", sep = "\t", 
                      col.names=c("contig","start","end","protein_ncbi","description",
                                  "ensembl_protein","ensembl_gene_id_version","gene_symbol"))


genes <- genes_t[,7] # Get column with Ensembl IDs


# Build query and get GO IDs and names 

go_list <- getBM(attributes=c("ensembl_gene_id_version", "go_id", "name_1006"),
                 filters = "ensembl_gene_id_version",
                 values = genes,
                 mart=ensembl)


# Remove rows with no GO names and get all terms per gene ID

go_list_clean <- go_list %>%
  mutate(name_1006 = na_if(name_1006, ""))

go_list_clean <- go_list_clean %>%
  filter(!is.na(name_1006))


go_per_gene <- go_list_clean %>%
  group_by(ensembl_gene_id_version) %>%
  summarise(
    go_id = paste(unique(go_id), collapse = ";"),
    go_name = paste(unique(name_1006), collapse = ";"))


# Merge with gene info

genes_with_go <- genes_t %>%
  left_join(go_per_gene, by = "ensembl_gene_id_version")


#Filter for specific names

pattern <- paste(c("temperature", "thermo", "lipo", "lipid", "sterol", 
                   "salinity" "osmo", "sodium", "potassium", "chloride", 
                   "\\banion\\b", "\\bcation\\b", "ion transport", "ion channel",
                   "monoatomic ion"), collapse = "|")

go_filtered <- genes_with_go %>%
  filter(grepl(pattern, go_name, ignore.case = TRUE))
```

#### 8.2.2 Test for overrepresention of relevant functional categories in high-ΔAF regions

 

Beforehand, different sets of peak windows were defined based on their
bootstrap support, regarded as the bootstrap support of the of the SNP
with the highest *Δ**A**F*, and by extending the main peak region with
genomic windows of varying lengths on either side (1Mbp, 500kbp or
100kbp). e.g., to get peaks with ≥ 50% bootstrap support and their 1Mbp
surrounding region:

``` r
library(dplyr)

pos_diff_55_50boot_1000 <- AF_diff_boot_55_sup %>% # Contains peaks, contigs, 
  # positions and bootstrap supports for SNP with highest ΔAF per peak
  filter(boot_max>=50) %>%
  group_by(CHROM, peak) %>%
  summarise(start = pmax(min(POS) - 1000000, 0), end = max(POS) + 1000000,
    .groups = "drop") %>%
  select(CHROM, start, end)

write.table(pos_diff_55_50boot_1000, "pos_diff_55_50boot_1000.tsv", 
            quote=FALSE, col.names=FALSE, row.names = FALSE, sep="\t")
```

**Genomic randomization approach**

We assessed whether the proportion of functionally relevant genes within
the peak windows exceeded that expected by chance.

Genes were considered as relevant for ion/osmo regulation or
thermoregulation if they were annotated with GO terms containing the
following salinity- or temperature-related keywords, respectively.

``` bash
sal='salinity|osmo|sodium|potassium|chloride|\banion\b|\bcation\b|ion transport \
|ion channel|monoatomic ion'

thermo='temperature|thermo|lipo|lipid|sterol'
```

Function to count all unique Ensembl gene IDs, genes with
functionally-relevant GO terms, and the proportion of functionally
relevant to total genes

``` bash
GENE_COL=5 #Column that contains Ensembl ferret gene IDs
GO_COL=8 #Column that contains GO annotations

count_genes() {
    local file="$1" pattern="$2"

    awk -F'\t' -v id="$GENE_COL" -v go="$GO_COL" -v pat="$pattern" '
    {gene=$id
        if (gene=="") next
        total[gene]=1
        if (tolower($go) ~ tolower(pat)) functional[gene]=1} #Each gene ID
        #is counted once
    END {for (g in total) t++
        for (g in functional) f++
        printf "%d\t%d\t%.8f\n", t, f, (t>0 ? f/t : 0)
    }' "$file"}
```

The observed proportion of salinity-related and temperature-related
genes in the different sets of peak windows was calculated. e.g., for
the set defined by peaks with ≥ 50% bootstrap support and their 1Mbp
surrounding region:

``` bash
# Extract genes found in the corresponding peak regions from the annotation file
# (containing ferret Ensembl IDs and GO annotations)

$bedtools intersect \
    -a elephant_seal_BUSCO_geneIDs_function_ensembl_GO.tsv \
    -b pos_diff_55_50boot_1000.tsv \
    -wa -u > observed_genes.bed


# Count functionally-relevant and total genes, and their proportion

## Salinity
count_genes observed_genes.bed "$sal" > observed_sal.tsv

## Temperature
count_genes observed_genes.bed "$thermo" > observed_temp.tsv
```

Genomic windows of the same number and length as the different sets of
peak windows were sampled 1000 times. Next, these windows were
intersected with the annotation file.

``` bash
peaks=pos_diff_55_50boot_1000.tsv #example

# Generate random windows

for i in $(seq 1 1000); do
    $bedtools shuffle \
        -i $peaks \
        -g reference/genome.txt \
        -excl $peaks \
          -noOverlapping \
        -seed "$i" \
        > permutations_50_1000/random_peaks_${i}.bed
done


# Intersect with annotation file and
# count functionally-relevant and total genes, and their proportion

for i in $(seq 1 1000); do
    $bedtools intersect \
        -a elephant_seal_BUSCO_geneIDs_function_ensembl_GO.tsv \
        -b permutations_50_1000/random_peaks_${i}.bed \
        -wa -u > permutations_50_1000/genes_${i}.bed
        
    ## Salinity
    sal_stats=$(count_genes permutations_50_1000/genes_${i}.bed "$sal")
    echo -e "$i\t$sal_stats" >> permutation_50_1000_results_sal.tsv

    ## Temperature
    thermo_stats=$(count_genes permutations_50_1000/genes_${i}.bed "$thermo")

    echo -e "$i\t$thermo_stats" >> permutation_50_1000_results_temp.tsv

done
```

# 9. Variant effect prediction and conservation scores

Before predicting the functional effect of variants segregating between
the subspecies, the ancestral state of alleles was inferred and the VCF
was repolarized. The version of the VCF repolarized was not masked and
included both biallelic SNPs and biallelic indels.

The ancestral state inference and repolarization of the VCF were done
using custom scripts.

**Script to add AA information:**

``` bash
cat > scripts/addAA.sh << 'EOF'

vcf=$1
temp=$2

bam=mirang_sorted.bam
ref=Ribbon_PR_18_6_2024.fa

$gatk4 HaplotypeCaller -I $bam -R $ref -L $vcf.vcf.gz -O temp/$temp.gvcf.gz \
--output-mode EMIT_ALL_CONFIDENT_SITES -ERC BP_RESOLUTION 

$gatk4 GenotypeGVCFs -R $ref -L $vcf.vcf.gz  --call-genotypes true \
--force-output-intervals $vcf.vcf.gz -V temp/$temp.gvcf.gz -O temp/$temp.vcf.gz

$bcftools query -i 'INFO/DP>4' -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' temp/$temp.vcf.gz \
| awk '{OFS="\t";printf "%s.%09d\t",$1,$2; print $0}' > temp/$temp.1.sites.tsv
$bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $vcf.vcf.gz \
| awk '{OFS="\t";printf "%s.%09d\t",$1,$2; print $0}' >  temp/$temp.2.sites.tsv

join -j1 temp/$temp.2.sites.tsv temp/$temp.1.sites.tsv \
| awk '$4==$8 && ($5==$9 || $9=="."){
    OFS="\t";
    if($10=="0/0" || $10=="0|0") {
    print $2,$3,$4,$5,$4
    } else if($10=="1/1" || $10=="1|1") {
    print $2,$3,$4,$5,$5
    } else {
    print $2,$3,$4,$5,"."
    }
}' >  temp/$temp.AA.sites.tsv

bgzip temp/$temp.AA.sites.tsv
tabix -s1 -b2 -e2 temp/$temp.AA.sites.tsv.gz

echo '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">' \
> temp/$temp.hdr.txt

$bcftools annotate -a temp/$temp.AA.sites.tsv.gz -c CHROM,POS,REF,ALT,INFO/AA -h \
temp/$temp.hdr.txt -Oz -o ${vcf}_AA.vcf.gz $vcf.vcf.gz

$bcftools view -m2 -M2 ${vcf}_AA.vcf.gz | perl scripts/polarizeVCF.pl \
| bcftools view - -Oz -o ${vcf}_bial_pol.vcf.gz
$bcftools index ${vcf}_bial_pol.vcf.gz
EOF
```

**Script to polarize VCF data:**

``` bash
cat > scripts/polarizeVCF.pl << 'EOF'

while(<>) {
    if(/^#/) {
        print $_;
        next;
    }
    chomp;
    my($chr,$site,$rs,$ref,$alt,$info) = (split(/\t/,$_))[0,1,2,3,4,7];
    my $ancestral_allele;
    $ancestral_allele = (split(/AA=/,$info))[1];
    $ancestral_allele = (split(/[\|;]/,$ancestral_allele))[0];
    $ancestral_allele = uc($ancestral_allele);

    if($ancestral_allele eq $ref) { #ok as it is
        print "$_\n";
        next;
    } else { #Switch alleles
        my @tmp = split(/;/,$info);
        foreach my $t ( @tmp ) {
            next if($t !~ /AF/);
            my($key,$val) = split(/=/,$t);
            my $new_AF = 1-$val;
            $_ =~ s/$key=$val/$key=$new_AF/;
        }
        my($chr,$site,$rs,$ref,$alt,$qual,$filter,$info,$format,@tmp) = split(/\t/,$_);
        my $swp = $ref;
        $ref = $alt; 
        $alt = $swp;
        my $str = join "\t",$chr,$site,$rs,$ref,$alt,$qual,$filter,$info,$format,"";
        foreach my $t ( @tmp ) {
            if($t =~ /^[01][\/\|][01]/) {
                my $sep='/';
                if($t =~ /\|/) { $sep='|'; }
                #($t) =~ s/\|/\//;
                my($a,$b,$c,$d,$e) = (split(/[\/\|:,]/,$t));
                if($a == 0)    {    $a = 1;    }
                elsif($a == 1) {    $a = 0;    }
                if($b == 0)    {    $b = 1;    }
                elsif($b == 1) {    $b = 0;    }

                if($a == 1 && $b == 0) { $a = 0; $b = 1; }
                
                $t = "$a$sep$b:$d,$c:$e";
            }
            $str .= "$t\t";
        }
        $str =~ s/\t$//;
        print "$str\n";
    }
}
EOF
```

These scripts were utilized to add AA and re-polarize sites where
‘INFO/AA=ALT’.

``` bash
bash scripts/addAA.sh grey_seg.vcf.gz mirang_seals
```

As multiallelic sites could not be polarized, the polarized *biallelic*
sites were subsequently combined with the non-polarized multiallelic
sites.

``` bash
$bcftools isec -p isec -c all -Oz grey_seg_bial_pol.vcf.gz grey_seg.vcf.gz 

$bcftools concat -a isec/0002.vcf.gz isec/0001.vcf.gz -Ou | $bcftools sort --temp-dir \
tmp_merge -Oz -o full_with_pol.vcf.gz

$bcftools index full_with_pol.vcf.gz
```

Prediction of variant effects was done with SnpEff. To this end, the
program files were copied into the local working directory and the input
data were structured according to the program’s requirements.

``` bash
mkdir -p SnpEff/data/ribbonseal

cp $snpeff SnpEff

cp ${refname}.fa genes.gff SnpEff/data/ribbonseal

cd SnpEff/data/ribbonseal 

mv ${refname}.fa sequences.fa
```

Next, the config file was edited and an annotation database was created.

``` bash
echo "ribbonseal.genome : ribbonseal" >> snpEff.config

java -jar snpEff.jar build -gff3 -v ribbonseal -noCheckCds  -noCheckProtein
```

And the VCF file was annotated.

``` bash
file_ann=grey_full_with_pol_ann

( java -jar snpEff.jar -v -stats $vcf/${file_ann}.html ribbonseal $vcf/full_with_pol.vcf.gz \
| $bcftools view -Oz -o $vcf/${file_ann}.vcf.gz ) &> $vcf/${file_ann}.log 
```

The VCF was additionally annotated with **Genomic Evolutionary Rate
Profiling (GERP) scores** obtained from Ensembl datasets for ferret.

First, the GERP scores were retrieved and converted to BED format.

``` bash
wget https://ftp.ensembl.org/pub/release-115/compara/conservation_scores/\
92_mammals.gerp_conservation_score/gerp_conservation_scores.mustela_putorius_furo.MusPutFur1.0.bw

$bigWigToBedGraph gerp_conservation_scores.mustela_putorius_furo.MusPutFur1.0.bw \
ferret_conservation_scores.bed
```

Then, a ferret-to-ribbon seal liftover chain was created and the GERP
scores were lifted over to the ribbon seal genome, similarly to 8.1.

``` bash
# Get ferret genome and create chain file

wget https://ftp.ensembl.org/pub/release-115/fasta/mustela_putorius_furo/dna/\
Mustela_putorius_furo.MusPutFur1.0.dna_sm.toplevel.fa.gz

gunzip Mustela_putorius_furo.MusPutFur1.0.dna_sm.toplevel.fa.gz

ferret=Mustela_putorius_furo.MusPutFur1.0.dna_sm.toplevel

$nextflow run evotools/nf-LO --source ${ferret}.fa --target ${refname}.fa --distance far \
--aligner minimap2 --outdir ferret_to_seal --no_netsynt

# Transfer scores

$swiftover -t bed -c $ferret_to_seal_chain -i ferret_conservation_scores.bed \
-u unmatched_ferret_to_seal.bed > ribbon_conservation_scores.bed
```

With this, the VCF could be annotated.

``` bash
echo '##INFO=<ID=GERP,Number=1,Type=Float,Description="GERP conservation score">' \
> gerp_header.txt

$bcftools annotate -a ribbon_conservation_scores.bed -h gerp_header.txt \
-c CHROM,FROM,TO,INFO/GERP -Oz -o grey_full_with_pol_ann_GERP.vcf.gz \
grey_full_with_pol_ann.vcf.gz
```

# References

Catalán A, Gygax D, Candolin U, Tusso S, Duchen P, Höhna S. 2025.
Sex-biased Migration and Demographic History of the Big European Firefly
Lampyris noctiluca. Molecular Biology and Evolution. 42:msaf12342.
10.1093/molbev/msaf123.

Olkkonen E, Löytynoja A. 2023. Analysis of population structure and
genetic diversity in low-variance Saimaa ringed seals using low-coverage
whole-genome sequence data. STAR Protoc. 4:1025674.
10.1016/j.xpro.2023.102567.

Grethlein M, Fekete Z, Goffart S, Kiebler A, Kunnasranta M, Niemi M, Santoro DF, Wehrenberg G, Winter S, Prost S, et al, unpublished data, https://www.biorxiv.org/content/10.64898/2026.08.13.744633v1, last accessed November 11, 2026
