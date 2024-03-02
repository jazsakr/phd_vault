internal link: [[Nanopore]]

---

We wanted to convert our older bulk RNA-seq runs with raw data as fast5 files to [[pod5]] then re-basecall them. 

> [!warning] 
> Code on this page works with a SINGLE final pod5 file and requires the package [[pigz]] for compression.

First, we move all the old fastq files created by guppy into a folder (you can either decide to keep them or delete them).
```bash
mkdir guppy_fastq_files
```

---

# Convert to pod5

Use the package [[pod5]] to convert fast5 files to pod5. Again, this code is based on a single, final pod5 file.

In the following example, the command will look recursively for fast5 files, output a single pod5 file and save the screen output and progress bar to a log file that will be parsed later.

```bash
file="sample"
pod5 convert fast5 --recursive ./ --output ${file}.pod5 2> ${file}_pod5.log 
```

NOTE: You can check progress and estimated completion by running `tail ${file}_pod5.log`

---

# pod5 conversion report

Next, we want to generate a report from the pod5.log file.

The follwoing code wll:
1. Get `pod5` package version.
2. Count the number of raw reads in the new pod5 file.
3. Count the number of fast5 files that did not convert to pod5 due to errors.
4. Count the number of converted reads (based on default 4000 reads per file).
5. Get file size of the final pod5 file.
6. Calculate the checksum of the final pod5 file.

```bash
echo "$(pod5 --version)" > ${file}_pod5_log_summary.txt && echo "counting raw reads" && echo -e "\ntotal raw reads" >> ${file}_pod5_log_summary.txt && echo "$(pod5 view ${file}.pod5 | grep -v read_id | wc -l)" >> ${file}_pod5_log_summary.txt && echo "counting errors" && echo -e "\nunconverted files" >> ${file}_pod5_log_summary.txt && grep -o 'ERROR:pod5:Encountered an exception in .*' ${file}_pod5.log | awk '{print $5}' | rev | cut -d'/' -f1 | rev >> ${file}_pod5_log_summary.txt && echo "total unconverted files $(grep fast5 ${file}_pod5_log_summary.txt | wc -l)" >> ${file}_pod5_log_summary.txt && echo "total unconverted reads $(expr $(grep fast5 ${file}_pod5_log_summary.txt | wc -l) \* 4000)" >> ${file}_pod5_log_summary.txt && echo "checking file size" && echo -e "\nfile size" >> ${file}_pod5_log_summary.txt && echo "$(du -sh ${file}.pod5)" >> ${file}_pod5_log_summary.txt && echo "creating checksum" && echo -e "\nchecksum" >> ${file}_pod5_log_summary.txt && md5sum ${file}.pod5 >> ${file}_pod5_log_summary.txt
```

Output will look something like this:

```
Pod5 version: 0.3.2

total raw reads
99117638

unconverted files

total unconverted files 2
total unconverted reads 8000

file size
990G	sample.pod5

checksum
5f179ef785e5992e33ff222eaf8fc3b4  sample.pod5
```

---

# Re-basecall

Then we re-basecall using [[dorado]].

The following is an example:
```bash
dorado basecaller ~/packages/dorado/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.1.0 ${file}.pod5 --device "cuda:0" --min-qscore 10 --emit-fastq > ${file}.fastq
```

---

# fastq report

Next we want to generate a report with information about our fastq file.

The following code will:
1. Get the `dorado` version. 
2. Count the number of reads in fastq file.
3. Compress the fastq file. 
4. Get file size of the fastq file.
5. Calculate the checksum of the fastq file. 

```bash
dorado --version 2> ${file}_fastq.log && sed -i 's/\+.*//' ${file}_fastq.log && echo -e "dorado version\n$(cat ${file}_fastq.log)" > ${file}_fastq.log && echo "counting reads in fastq" && echo -e "\nreads in fastq" >> ${file}_fastq.log && echo $(cat ${file}.fastq|wc -l)/4|bc >> ${file}_fastq.log && echo "compressing fastq" && pigz -p 12 ${file}.fastq && echo "checking file size" && echo -e "\nfile size" >> ${file}_fastq.log && echo "$(du -sh ${file}.fastq.gz)" >> ${file}_fastq.log && echo "creating checksum" && echo -e "\nchecksum" >> ${file}_fastq.log && md5sum ${file}.fastq.gz >> ${file}_fastq.log
```

The output will look something like this:
```
dorado version
0.5.0

reads in fastq
86620316

file size
59G	sample.fastq.gz

checksum
fafa73fe0c3a3dc05f9e00d874eebf60  sample.fastq.gz
```