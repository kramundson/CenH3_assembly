Aim: Assemble CenH3 from Desiree short reads. We want to be able to figure out CenH3 copy
number and design allele-specific primers for two CRISPR gRNA targets, M4 and M14.

Datasets:
1. Sanger sequenceing reads from Sundaram
2. Illumina 150nt PE WGS of 37 T0 lines, approx 7M reads per individual.

Strategy:
1. Align reads to reference genome + alignment QC (all were done in a Snakemake workflow for a different project)
2. Estimate CenH3 gene copy number from sequencing depth
3. Extract reads that align to annotated CenH3 and their mates
4. Merge overlapping paired-end reads
5. De novo assembly of CenH3-mapped reads + Sanger reads

1. Align reads to reference genome + alignment QC

See ```/share/comailab/kramundson/ximena/tetraploid_SNPs``` for Snakemake workflow.
To give an idea of what was done, I'll include example commands from snakemake here:

Trim adapter and low-quality sequence from raw reads:

```
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 10 -m 40 -j 8 -o {output} {input.fq1} {input.fq2}
```

BWA mem alignment to reference genome:

```
"bwa mem -R {params.rg} -t {params.bwa_threads} {input.ref} {input.reads} | "
"samtools sort -@{params.sort_threads} -m {params.sort_mem} -o {output} -"
```

Note: reference used for alignment was ```/share/comailab/kramundson/ximena/ref/potato_dm_v404_all_pm_un_chloro_mito.fasta```
md5sum: 26694987d696f73939143ae3cc9dd64b

Remove PCR duplicates from each library using Picard MarkDuplicates:

```
"java {params.java_heap} -jar {params.jarpath} MarkDuplicates INPUT={input} "
"OUTPUT={output.bam} METRICS_FILE={output.metrics} {params.opt} 2> {log}"
```

Realign indels for each library using GATK RealignerTargetCreator and IndelRealigner

```
# create targets
gatk -T RealignerTargetCreator \
    -R {input.ref} \
    -I {input.bam} \
    -o {output} \
    2> {log}
    
# realign target indels
gatk -T IndelRealigner \
    -R {input.ref} \
    -I {input.bam} \
    -targetIntervals {input.intervals} \
    -o {output} \
    2> {log}
```

For overlapping mates, clip one of the two mates in the overlap region

```
bam clipOverlap --in {input} --out {output} --stats --unmapped 2> {log}
```

Finally, merge all samples from Cornell-derived Desiree

```
samtools merge Cornell_Desiree.bam Cornell_Desiree*.bam && samtools index && Cornell_Desiree.bam
```

2. Estimate gene copy number from sequencing depth

```
samtools depth -r chr01:73237104-73244176 Cornell_Desiree.bam > Cornell_Desiree_CenH3_depth.tsv
samtools depth -r chr01:73237104-73244176 -Q 30 Cornell_Desiree.bam > Cornell_Desiree_CenH3_Q30depth.tsv
```

Roughly estimating from samtools flagstat output, the genome-wide mean sequencing depth is ~60x

340268146 aligned reads (each 150nt)
340268146 * 150 / 844e6 = 60.4742

Note: this is a rough estimate that doesn't take things like overlapping mates into account.

![alt text](./images/IGV_screenshot_CenH3_Desiree.png)

Figure: IGV screenshot of CenH3. Aligned reads are from 37 T0 lines, each of which was
sequenced to low coverage with 150nt PE reads.

![alt text](./.images/prelim_CenH3_dp.png)

Figure: Green line shows roughly estimated genome-wide mean sequencing depth (60). In blue
is the per-position sequencing depth using reads of any mapping quality. In red is the
per-position sequencing depth using only with mapping quality Q30 or better. Y-axis breaks
are in increments of 15, which reflects the expected change in read depth for integer
copy number changes of CenH3 genes. For example, read depth of 75 suggests 5 CenH3 copies
in 4 monoploid genomes.

Overall, I'm having a hard time interpreting this figure.
There is definitely increased read depth over the CenH3 gene body compared to the average
genome-wide depth, but the depth increase is not uniform throughout the gene. Maybe this
could be due to extra copies of CenH3 genes and/or fragments of pseudogenes in Desiree.
Could this be cultivar-specific or potato-specific?

3. Extract reads aligning to CenH3 gene and mates.

```
samtools view Cornell_Desiree.bam chr01:73237104-73244176 | cut -f 1 | sort -u > 2019_0226_Des_CenH3_headers.txt

zgrep -F -A 3 -f 2019_0226_Des_CenH3_headers.txt ../reads/*SK*-1.fq.gz | \
    sed -e 's/..\/reads\/KFRAG_[0-9]\{5,\}_SK_*[0-9]\{1,\}-[0-9]\{1,\}.fq.gz://g' | \
    grep -v "^--$" > Des_CenH3_1.fastq 2> Des_CenH3_1.err &

zgrep -F -A 3 -f 2019_0226_Des_CenH3_headers.txt ../reads/*SK*-2.fq.gz | \
    sed -e 's/..\/reads\/KFRAG_[0-9]\{5,\}_SK_*[0-9]\{1,\}-[0-9]\{1,\}.fq.gz://g' | \
    grep -v "^--$" > Des_CenH3_2.fastq 2> Des_CenH3_2.err &
```

```
echo $(($(wc -l Des_CenH3_1.fastq | cut -d ' ' -f 1)/4))
echo $(($(wc -l Des_CenH3_2.fastq | cut -d ' ' -f 1)/4))
```

3,367 paired end reads were extracted

4. Merge overlapping paired-end reads

```
pear -f Des_CenH3_1.fastq -r Des_CenH3_2.fastq -o Des_CenH3_pear -j 2 &> Des_CenH3_pear.log
```

Output is a set of four fastq files:
 * .assembled.fastq extension for merged reads
 * .unassembled.1.fastq extension for unmerged forward mates
 * .unassembled.2.fastq extension for unmerged reverse mates
 * .discarded.fastq extension, no reads are in this file
 
5. Prep Sanger sequencing reads.

Note: Spades does not run error correction on Sanger, Oxford Nanopore, or PacBio reads.
For now, I'll provide untrimmed Sanger reads as a single multi-FASTA file.

```
# local
cat Sundaram_CenH3_Sanger_Sequencing/*CenH3*check*/Des*.seq > Des_CenH3_sanger.fasta
# upload Des_CenH3_sanger.fasta to buffy for assembly
```

6. SPAdes assembly

```
module load spades/3.13.0
spades.py -1 Des_CenH3_pear.unassembled.forward.fastq -2 Des_CenH3_pear.unassembled.reverse.fastq \
    --merged Des_CenH3_pear.assembled.fastq --sanger Des_CenH3_sanger.fasta -o 2_spades_with_pear_sanger
    
# copy to local, saved as ./spades_assemblies/2_spades_with_pear_sanger/
```

7. In BANDAGE package tBLASTn M4 and M14 target sites to scaffolds. Want to know which
scaffolds in the assembly have the gRNA target sites, whether there are mismatches at
targets, and whether there are mismatches near targets.

I had to adjust BLAST params for short query sequences.
Params used:
-word_size 7 -reward 1 -penalty -3
Filter 90% similarity over 85% length

From the graph output, it looks like there are locations in Sundaram's sanger reads where
I have four haplotigs. If I can overlay the CRISPR target site, it might be possible to
design allele specific primers.

TODO show clustal output. Mismatches in protospacer. PAM intact for all scaffolds.
Looks like M14 allele specific primers might be possible. M4 probably not without better
assembly. Need to design primers that are far enough away to amplify and Sanger sequence.

