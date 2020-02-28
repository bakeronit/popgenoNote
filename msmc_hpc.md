## Estimate change of effective population size using MSMC

1. Prepare VCF file using samtools. (I used all samples to increase the accuracy of phasing step)

**For each scaffold:**
```
samtools mpileup -q 20 -Q 20 -C 50 -u -r $chr -f $reference $sample1.bam $sample2.bam ... |bcftools call -c -V indels | bcftools view --max-alleles 2 |bgzip > vcf_files/${chr}.vcf.gz
```

2.Read aware phasing with Shapeit2

2.1 Prepare bam list file for every scaffold.

```{bash}
sample1 bam1    chr1
sample2 bam2    chr1
```
2.2 Extract phase informative reads (by scaffold).

```{bash}
extractPIRs --bam ${chr}.bamlist --vcf vcf_files/${chr}.vcf.gz --out PIRs/${chr}.PIRsList --base-quality 20 --read-quality 20
```
2.3 Phasing

```{bash}
shapeit -assemble --input-vcf vcf_files/${chr}.vcf.gz --input-pir PIRs/${chr}.PIRsList -O haplotypeData/$chr

shapeit -convert --input-haps haplotypeData/$chr --output-vcf phasedVCF/${chr}.phased.vcf.gz
```
3.Merge phased vcf and unphased vcf

```{bash}
bcftools merge --force-samples vcf_files/${chr}.vcf.gz phasedVCF/${chr}.phased.vcf.gz |
awk 'BEGIN {OFS="\t"}'
    $0 ~/^##/ {print $0}
    $0 ~/^#CHROM/ {for(i=1;i<84;i++) printf "%s"OFS, i; print $84}
    $0 !~/^#/ {for(i=10;i<=84;i++) $i=$(i+75); for(i=1;i<84;j++) printf "%s"OFS,i; print $84} |
 bcftools view -Oz > finalVCF/${chr}.vcf.gz
 tabix -p vcf finalVCF/${chr}.vcf.gz
```
I picked 7 individuals to do the following msmc analysis.

4. Combine vcf for all samples I picked.

4.1 generate input vcf and mask bed file

```{bash}
#calculate mean coverage for each sample on each scaffold
mean_cov=$(samtools depth -r $chr $bam| awk '{{sum += $3}} END {{print sum/NR}}')

#get masked vcf
bcftools view -s $sample finalVCF/${chr}.vcf.gz | \
bamCaller.py $mean_cov ${sample}/${chr}.mask.bed.gz |gzip -c > ${sample}/${chr}.vcf.gz
```
This created 7\*970 files...

4.2.Create mappability file

```{bash}
./bin/seqbility-20091110/splitfa $reference 35 | split -l 20000000

ls x* | while read id; do bwa aln -t 20 -R 1000000 -O 3 -E 3 $reference $id > ${id}.sai;done

ls x??| while read id; do bwa samse -f ${id}.sam $reference ${id}.sai ${id} && gzip ${id}.sam;done

wait

gzip -dc x??.sam.gz | ./bin/seqbility-20091110/gen_raw_mask.pl > $raw_fasta

./bin/seqbility-20091110/gen_mask -l 35 -r 0.5 $raw_fasta > $masked_fasta

#edit makeMappabilityMask.py first
python2 ./bin/msmc-tools/makeMappabilityMask.py
```

4.3 Combine all

```{bash}
generate_multihetsep.py --chr $chr \
--mask samples/AI_1_022/${chr}_mask.bed.gz \
--mask samples/BR_4_091/${chr}_mask.bed.gz \
--mask samples/RS1_S_314/${chr}_mask.bed.gz \
--mask samples/RS2_C11_784/${chr}_mask.bed.gz \
--mask samples/RS3_1_184/${chr}_mask.bed.gz \
--mask samples/AR_125_388/${chr}_mask.bed.gz \
--mask samples/AR_133_357/${chr}_mask.bed.gz \
--mask mappability/Adigi_${chr}.mask.bed.gz \
samples/AI_1_022/${chr}.vcf.gz \
samples/BR_4_091/${chr}.vcf.gz \
samples/RS1_S_314/${chr}.vcf.gz \
samples/RS2_C11_784/${chr}.vcf.gz \
samples/RS3_1_184/${chr}.vcf.gz \
samples/AR_125_388/${chr}.vcf.gz \
samples/AR_133_357/${chr}.vcf.gz > multihetsep/indv7.${chr}.multihetsep.txt
```

5.run MSMC

5.1 run msmc
```{bash}
msmc -t 20 -R -o indv7 multihetsep/indv7.*.multihetsep.txt
```
This takes long time, so I cancle the job and turn to msmc.

5.2 run msmc2
Run msmc2 for all individuals in every population
```{bash}
msmc2 -t 10 -s -I 0,1,2,3 -o inshore_4hap multihetsep/indv7.*.multihetsep.txt

msmc2 -t 10 -s -I 4,5,6,7,8,9 -o southoffshore_6hap multihetsep/indv7.*.multihetsep.txt

msmc2 -t 10 -s -I 10,11,12,13 -o northoffshore_4hap multihetsep/indv7.*.multihetsep.txt
```

