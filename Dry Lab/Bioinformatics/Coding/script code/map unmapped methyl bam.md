This script will take unmapped bam files with modification information and map using [[samtools]] and [[minimap|minimap2]]. Then statistics and coverage report will be generated. 

:RiInformationLine: **NOTE**: For detailed explanation on mapping bam files, see [[minimap#Using bam files as input]]. For detailed explanation on how the reports are created, see [[samtools#Quick Reports]].

> [!warning]
> This script requires `minimap` and `samtools` to be installed


```bash
file="sample_unmapped"
fileout="sample"

path="~/analysis/experiment/bams/"

ref="~/references/chm13.fa"

# creat a log file
touch ${path}${fileout}.log

# get minimap and samtools version
n=$(minimap2 --version) && echo "minimap2 $n" >> ${path}${fileout}.log && samtools --version | head -1 >> ${path}${fileout}.log

# map unmapped methyl bam file
samtools fastq --threads 64 -T MM,ML ${path}${file}.bam | minimap2 -t 64 -ax map-ont --secondary=no -y ${ref} - 2>> ${path}${fileout}.log | samtools sort --threads 64 - > ${path}${fileout}.bam && \
samtools index -@ 64 ${path}${fileout}.bam

# get mapping statistics
samtools stats --threads 64 ${path}${fileout}.bam | grep ^SN | cut -f 2- | grep -e 'raw total sequences' -e 'reads mapped' -e 'reads unmapped' -e 'supplementary alignments' \
-e 'total length' -e 'bases mapped (cigar):' -e 'average length' -e 'maximum length' | awk '/raw total sequences:/ { total_reads=$4 } /reads mapped:/ { mapped_reads=$3 } \
/reads unmapped:/ { unmapped_reads=$3 } /supplementary alignments:/ { supp_alignments=$3 } /total length:/ { total_length=$3; total_gb=total_length/1000000000 } \
/cigar/ { mapped_bases=$4; mapped_gb=mapped_bases/1000000000 } /average length:/ { avg_length=$3 } /maximum length:/ { max_length=$3 } END \
{ printf "total reads:\t%s (%.2f M)\nreads mapped:\t%s (%.2f M)\nreads unmapped:\t%s (%.2f M)\nsupplementary alignments:\t%s (%.2f M)\ntotal length:\t%s \
(%.2f Gb)\nbases mapped:\t%s (%.2f Gb)\naverage length:\t%s\nmaximum length:\t%s\n", total_reads, total_reads/1000000, mapped_reads, mapped_reads/1000000, \
unmapped_reads, unmapped_reads/1000000, supp_alignments, supp_alignments/1000000, total_length, total_gb, mapped_bases, mapped_gb, avg_length, max_length }' \
> ${path}${fileout}_stats.txt

# get coverage info
samtools coverage ${path}${fileout}.bam | awk 'NR==1 {print $1, $4, $6, $7; next} {printf "%s\t%s\t%.2f%%\t%.2fx\n", $1, $4, $6, $7}' > ${path}${fileout}_coverage_unsorted.tsv \
&& echo -e "chrom\tnumreads\tcoverage\tmeandepth" > ${path}${fileout}_coverage.tsv && cat ${path}${fileout}_coverage_unsorted.tsv | tail -n +2 | sort -k1,1V | tr ' ' '\t' >> \
${path}${fileout}_coverage.tsv && rm ${path}${fileout}_coverage_unsorted.tsv

# update the mapping statistics
awk 'NR > 1 && $1 != "chrM" && $1 != "chrX" && $1 != "chrY" { sub(/x$/, "", $4); sum += $4; count++ } END { printf "Average depth: %.2fx\n", sum / count }' ${path}${fileout}_coverage.tsv >> ${path}${fileout}_stats.txt
```