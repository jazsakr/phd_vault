---
aliases:
  - minimap2
tags: []
type: package
status: maintained
---

Paper: [Minimap2: pairwise alignment for nucleotide sequences](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778)

GitHub: [lh3](https://github.com/lh3)/**[minimap2](https://github.com/lh3/minimap2)**

Documentation: [_Manual Reference Pages  -_ minimap2](https://lh3.github.io/minimap2/minimap2.html)

---

**Description**: Mapping and alignment between collections of DNA sequences.

---

# Install

```bash
conda install bioconda::minimap
# or
conda install -c bioconda minimap
```

---

# Usage

The following mapping parameters are based on `minimap` [recommendations](https://github.com/lh3/minimap2?tab=readme-ov-file#getting-started) for [[Nanopore]].

NOTE: Default output of `minimap` is a SAM file, which takes up a lot of space. It is recommend to pipe `minimap`’s output into `samtools` to sort then compress to a bam file and then index. See [[samtools]] page for more information.

## RNA

[Map long mRNA/cDNA reads](https://github.com/lh3/minimap2?tab=readme-ov-file#map-long-mrnacdna-reads)
**Note**: To prioritize annotated splice junctions when mapping, you need to provide a bed file of annotations. 

Use minimap2’s utility [paftools.js](https://github.com/lh3/minimap2/blob/master/misc/README.md#introduction) to convert the annotation file from gtf to bed file.
```bash
paftools.js gff2bed annotation.gff > annotation.bed
```

Then map:
```bash
# unsorted SAM file
minimap2 -t 24 -ax splice --secondary=no --MD --junc-bed annotation.bed reference_genome.fa sample.fastq.gz > output.sam 2> sample.log

# sorted BAM file
minimap2 -t 24 -ax splice --secondary=no --MD --junc-bed annotation.bed reference_genome.fa sample.fastq.gz 2> sample.log | samtools sort --threads 24 - > sample.bam && samtools index -@ 24 sample.bam
```

## map genomic DNA
[Map genomic reads](https://github.com/lh3/minimap2?tab=readme-ov-file#map-long-noisy-genomic-reads)

```bash
# unsorted SAM file
minimap2 -t 24 -ax map-ont --secondary=no reference_genome.fa sample.fastq.gz > sample.sam 2> sample.log

# sorted BAM file
minimap2 -t 24 -ax map-ont --secondary=no reference_genome.fa sample.fastq.gz 2> sample.log | samtools sort --threads 24 - > sample.bam && samtools index -@ 24 sample.bam
```

## Using bam files as input

For samples with base modifications, they are stored in bam files since fastq file format does not allow this information to be recorded. This is an issue when using `minimap` because it takes FASTQ files as input. There has been a [feature request for using sam/bam files as input](https://github.com/lh3/minimap2/issues/870#issuecomment-1033306393) and a current work around is [copying over the SAM flag through to the mapped output](https://github.com/igvteam/igv/issues/1435#issuecomment-1813318839).

First, `samtools fastq -T` will copy the SAM tags (`MM` for the base modification and `ML` for the modification probability) to `minimap` and then the `-y` `minimap` option will copy all tags through to the final output file.

```bash
# unsorted SAM file
samtools fastq --threads 24 -T MM,ML sample_unsorted.bam | minimap2 -t 24 -ax map-ont --secondary=no -y reference_genome.fa - > sample.sam

# sorted BAM file
samtools fastq --threads 24 -T MM,ML sample_unsorted.bam | minimap2 -t 24 -ax map-ont --secondary=no -y reference_genome.fa - | samtools sort --threads 24 - > sample.bam && samtools index -@ 24 sample.bam
```