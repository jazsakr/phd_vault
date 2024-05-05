---
aliases: 
tags: 
type: package
status: maintained
---

Documentation:
- [Documentation to all `samtools` commands](http://www.htslib.org/doc/samtools.html) (lists and links all `samtools` subcommand documentations)
- [The SAM Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf) (describes the format of a SAM/BAM file)
- [The Optional Fields Specification](https://samtools.github.io/hts-specs/SAMtags.pdf) (describes what the tags in the SAM/BAM file)

Resource:
- [Decoding SAM flags](https://broadinstitute.github.io/picard/explain-flags.html) (lets you enter SAM flag value and tell you what properties are associated with that value)
- [Tutorial:Piping With Samtools, Bwa And Bedtools](https://www.biostars.org/p/43677/) (examples of piping different samtools commands together)

---

# Install
## conda

```bash
conda install bioconda::samtools
# or
conda install -c bioconda samtools
```

---

# Usage

## `view`

```bash
# filter out any unmapped and not primary alignments
samtools view -F 260 -h -b --threads 12 sample.bam > sample_filtered.bam 

# filter bam file for specific reads
samtools view -h -b --threads 12 sample.bam -N read_ids.txt > sample_filtered.bam

# count primary reads
samtools view -F 260 -c sample.bam

# remove a specific tag
samtools view --threads 12 -b -h --remove-tag mv sample_tmp.bam > sample.bam

# filter reads then sort them (useful if the basecaller maps reads but not sort them)
samtools view -h -F 260 --threads 6 sample_unsorted.bam | samtools sort -@ 6 -o sample.bam - && samtools index -@ 12 sample.bam
```

## `merge`

```bash
# combine bam files while keeping track which reads came from which
samtools merge -o final.bam --threads 12 -r sample1.bam sample2.bam
```

## `cat`
```bash
# concat unmapped bam files from list of bam files
samtools cat --threads 12 -b bams_to_concat.txt -o sample.bam
```

## `reset`
```bash
# remove mapping information from bam file (helpful if you want to remap without rebasecalling)
samtools reset --threads 12 sample.bam | samtools view -b --threads 12 - > sample_unmapped.bam
```

## `stats`
```bash
# outputs statistics of your bam file
samtools stats --threads 12 sample.bam > sample_stats.txt
```

## `coverage`
```bash
# produces a histogram or table of coverage per chromosome.
samtools coverage sample.bam > sample_coverage.txt
```

---

## Piping `samtools` commands together

```bash
# sorting an unsorted bam file
samtools view -h sample_unsorted.sam | samtools sort -o sample.bam - && samtools index -@ 12 sample.bam
```

---

# Quick Reports

##  `stats` report

`samtools stats` gives a very comprehensive statistics report. Sometimes, I just need some basic information. For example, if I just want to know some mapping information about my samples, I would run some code to get output that would look something like this:

```
total reads:	7460016 (7.46 M)
reads mapped:	5747566 (5.75 M)
reads unmapped:	1712450 (1.71 M)
supplementary alignments:	1015394 (1.02 M)
total length:	63115221520 (63.12 Gb)
bases mapped:	57290692041 (57.29 Gb)
average length:	8460
maximum length:	1384696
```

The following code will:
1. Run `samtools stats`.
2. Look for “Summary numbers (SN)” section.
3. Find the following fields:  raw total sequences, reads mapped, reads unmapped, supplementary alignments, total length, bases mapped (cigar), average length, maximum length (refer to the documentation for a description of each).
4. Convert number reads to millions and add to end of the line in parentheses.
5. Convert bases to gigabases and add to the end of the line in parentheses.

```bash
# code as block for script
samtools stats --threads 12 sample.bam | grep ^SN | cut -f 2- | grep -e 'raw total sequences' -e 'reads mapped' -e 'reads unmapped' -e 'supplementary alignments' \
-e 'total length' -e 'bases mapped (cigar):' -e 'average length' -e 'maximum length' | awk '/raw total sequences:/ { total_reads=$4 } /reads mapped:/ { mapped_reads=$3 } \
/reads unmapped:/ { unmapped_reads=$3 } /supplementary alignments:/ { supp_alignments=$3 } /total length:/ { total_length=$3; total_gb=total_length/1000000000 } \
/cigar/ { mapped_bases=$4; mapped_gb=mapped_bases/1000000000 } /average length:/ { avg_length=$3 } /maximum length:/ { max_length=$3 } END \
{ printf "total reads:\t%s (%.2f M)\nreads mapped:\t%s (%.2f M)\nreads unmapped:\t%s (%.2f M)\nsupplementary alignments:\t%s (%.2f M)\ntotal length:\t%s \
(%.2f Gb)\nbases mapped:\t%s (%.2f Gb)\naverage length:\t%s\nmaximum length:\t%s\n", total_reads, total_reads/1000000, mapped_reads, mapped_reads/1000000, \
unmapped_reads, unmapped_reads/1000000, supp_alignments, supp_alignments/1000000, total_length, total_gb, mapped_bases, mapped_gb, avg_length, max_length }' \
> sample_stats.txt

# one liner for the command line
samtools stats --threads 12 sample.bam | grep ^SN | cut -f 2- | grep -e 'raw total sequences' -e 'reads mapped' -e 'reads unmapped' -e 'supplementary alignments' -e 'total length' -e 'bases mapped (cigar):' -e 'average length' -e 'maximum length' | awk '/raw total sequences:/ { total_reads=$4 } /reads mapped:/ { mapped_reads=$3 } /reads unmapped:/ { unmapped_reads=$3 } /supplementary alignments:/ { supp_alignments=$3 } /total length:/ { total_length=$3; total_gb=total_length/1000000000 } /cigar/ { mapped_bases=$4; mapped_gb=mapped_bases/1000000000 } /average length:/ { avg_length=$3 } /maximum length:/ { max_length=$3 } END { printf "total reads:\t%s (%.2f M)\nreads mapped:\t%s (%.2f M)\nreads unmapped:\t%s (%.2f M)\nsupplementary alignments:\t%s (%.2f M)\ntotal length:\t%s (%.2f Gb)\nbases mapped:\t%s (%.2f Gb)\naverage length:\t%s\nmaximum length:\t%s\n", total_reads, total_reads/1000000, mapped_reads, mapped_reads/1000000, unmapped_reads, unmapped_reads/1000000, supp_alignments, supp_alignments/1000000, total_length, total_gb, mapped_bases, mapped_gb, avg_length, max_length }' > sample_stats.txt
```

## `coverage` report

I want to filter the coverage report for chromosome, number of reads, percent of the chromosome covered and the depth coverage per chromosome as well as make it more human readable like this:

```
chrom numreads coverage meandepth
chr1	457785	99.99%	17.08x
chr2	364249	99.99%	15.12x
chr3	321696	99.85%	15.09x
```

The following code will:
1. Select the columns: \#rname, numreads, coverage, meandepth and save in a temp file.
2. Add column names to final report file.
3. Ignore the column names in the temp file and short the rows by `chr`.
4. Append sorted rows to the final report file.
5. Delete the temporary file.

```bash
# code as block for script
samtools coverage sample.bam | awk 'NR==1 {print $1, $4, $6, $7; next} {printf "%s\t%s\t%.2f%%\t%.2fx\n", $1, $4, $6, $7}' > sample_coverage_unsorted.tsv 
echo -e "chrom\tnumreads\tcoverage\tmeandepth" > sample_coverage.tsv && cat sample_coverage_unsorted.tsv | tail -n +2 | sort -k1,1V | tr ' ' '\t' >> sample_coverage.tsv
rm sample_coverage_unsorted.tsv

# one liner for the command line
samtools coverage sample.bam | awk 'NR==1 {print $1, $4, $6, $7; next} {printf "%s\t%s\t%.2f%%\t%.2fx\n", $1, $4, $6, $7}' > sample_coverage_unsorted.tsv && echo -e "chrom\tnumreads\tcoverage\tmeandepth" > sample_coverage.tsv && cat sample_coverage_unsorted.tsv | tail -n +2 | sort -k1,1V | tr ' ' '\t' >> sample_coverage.tsv && rm sample_coverage_unsorted.tsv
```

##  Update`stats` report

After creating my coverage report, I want to get the average coverage and add it to the statistics report I generated earlier. 

The following code will:
1. Omit the lines with header, chrM, chrX, and chrY
2. Get the fourth column with coverage and drop the “x” character
3. Calculated the average coverage
4. Add the average coverage to the statistics report generated previously


```bash
# update the stats.txt to include average coverage from coverage.tsv
awk 'NR > 1 && $1 != "chrM" && $1 != "chrX" && $1 != "chrY" { sub(/x$/, "", $4); sum += $4; count++ } END { printf "Average depth: %.2fx\n", sum / count }' sample_coverage.tsv >> sample_stats.txt
```

The updated report will look something like this:

```
total reads:	7460016 (7.46 M)
reads mapped:	5747566 (5.75 M)
reads unmapped:	1712450 (1.71 M)
supplementary alignments:	1015394 (1.02 M)
total length:	63115221520 (63.12 Gb)
bases mapped:	57290692041 (57.29 Gb)
average length:	8460
maximum length:	1384696
average coverage: 22.34x
```

---

# Errors

## no version information available
```
samtools: /home/user/miniconda3/bin/../lib/libtinfow.so.6: no version information available (required by samtools)
samtools: /home/user/miniconda3/bin/../lib/libncursesw.so.6: no version information available (required by samtools)
samtools: /home/user/miniconda3/bin/../lib/libncursesw.so.6: no version information available (required by samtools)
```

[libtinfo.so.6: no version information available message](https://stackoverflow.com/questions/72103046/libtinfo-so-6-no-version-information-available-message-using-conda-environment) solved by installing the following
```bash
conda install -c conda-forge ncurses
```