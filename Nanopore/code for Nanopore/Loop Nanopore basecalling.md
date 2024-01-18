internal link: [[Nanopore]]

---

One of our samples kept [aborting during basecalling](https://github.com/nanoporetech/dorado/issues/548), and a suggestion was to [loop through each file the pod5 folder and basecall each individually](https://github.com/nanoporetech/dorado/issues/548#issuecomment-1888051514). Here is the script I used to do that. 

**NOTE**: This example creates bams as output. For fastqs, replace all instances of `bam` with `fastq` and edit `dorado_options`  accordingly.

# Check variables

To check variables before running basecalling loop on all pod5 files, the the following on the command line. 
It will find just 3 pod5 files, and for each print:
1. pod5 file name
2. name of bam output file from that pod5 file
3. pod5 file and path
4. path of bam output folder
5. bam file and path

```bash
for pod5_file_check in $(find ./ -type f -name '*.pod5' | head -3); do pod5_name_check="$(basename "${pod5_file_check}")"; echo $pod5_name_check; bam_file_check="$(basename "${pod5_file_check%.pod5}.bam")"; echo $bam_file_check; echo $pod5_file_check; output_dir_check="$(dirname "${pod5_file_check}" | sed 's/pod5/bam/')"; echo $output_dir_check; output_bam_check="${output_dir_check}/${bam_file_check}"; echo -e "$output_bam_check\n"; done
```

---

# Loop Script

If the variables look correct, edit the dorado variables and run this script in the run folder. This script will 
1. Find all the pod5 files and create all the variables (it is the same code as above).
2. Create the bam output folder if it doesn’t exist (it will be created in the same folder where pod5 output folder is).
3. If the output bam file does not exist, the basecall the corresponding pod5 file with `dorado`
4. Save dorado progess messages in a log file called `basecalling.log`
5. Remove lines in the `basecalling.log` file that contain the string `info` WITHOUT the string `Finished`. The will help keep the file from getting too big (especially for larger runs). This will save pod5 file name, `Finished` messages, and all other messages like `error` and `warning`. 

**NOTE**: Run this script where you’d like the `basecalling.log` to be saved.

```bash
#!/bin/bash

# dorado parameters
dorado_command="$HOME/packages/dorado/dorado-0.5.1-linux-x64/bin/dorado basecaller"
dorado_options="--modified-bases 5mCG_5hmCG  --device "cuda:0" "
dorado_config="$HOME/packages/dorado/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0"

for pod5_file in $(find ./ -type f -name '*.pod5'); do
	# set files names and directories
	pod5_name="$(basename "${pod5_file}")"
	bam_file="$(basename "${pod5_file%.pod5}.bam")"
	output_dir="$(dirname "${pod5_file}" | sed 's/pod5/bam/')"
	output_bam="${output_dir}/${bam_file}"

	# make a directory for bam files if it does not exist yet
	[ ! -d "${output_dir}" ] && mkdir -p "${output_dir}"

	# if the bam file does not already exist:
	if [ ! -e "${output_bam}" ]; then
		echo "Basecalling $pod5_name" 2>&1 | tee -a basecalling.log
		# basecall pod5 file with dorado and save progess statements in log file
		${dorado_command} ${dorado_config} "${pod5_file}" ${dorado_options} > "${output_bam}" 2>> >(tee -a basecalling.log >&2)
		# clean-up log file to remove extra information
		sed -i '/info/ { /Finished/!d; }' basecalling.log
		echo "---" >> basecalling.log
	fi
done

```

Alternatively, here it is as a one liner:

```bash
dorado_command="$HOME/packages/dorado/dorado-0.5.1-linux-x64/bin/dorado basecaller"; dorado_options="--modified-bases 5mCG_5hmCG --device cuda:0"; dorado_config="$HOME/packages/dorado/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0"; for pod5_file in $(find ./ -type f -name '*.pod5'); do pod5_name="$(basename "${pod5_file}")"; bam_file="$(basename "${pod5_file%.pod5}.bam")"; output_dir="$(dirname "${pod5_file}" | sed 's/pod5/bams/')"; output_bam="${output_dir}/${bam_file}"; [ ! -d "${output_dir}" ] && mkdir -p "${output_dir}"; [ ! -e "${output_bam}" ] && { echo "Basecalling $pod5_name" 2>&1 | tee -a basecalling.log; ${dorado_command} ${dorado_config} "${pod5_file}" ${dorado_options} > "${output_bam}" 2>> >(tee -a basecalling.log >&2); sed -i '/info/ { /Finished/!d; }' basecalling.log && echo "---"  >> basecalling.log; }; done
```

---

# Quick Report

We want to know which files were basecalled successfully and which had errors. Instead of scrolling through the log, this line of code will do the following:
1. Find all the pod5 files that finished successfully and save file names in a temporary file called `basecalling_temp.log`.
2. Compare the basecalling log and the temporary file to find all pod5 files that gave an error.
3. Save the name of all the pod5 files that gave an error in a file called `reads_causing_errors.txt`.
4. Use the `basecalling_temp.log` file to create a successfully basecalled bam file list called `completed_bams.txt`. 
5. Report the number of total pod5 files, pod5 files that gave error and successfully completed bam files in a file called `basecalling_report.txt`.

**NOTE**: Run this line of code where the `basecalling.log` is saved.

```bash
grep -B 1 "Finished" basecalling.log > basecalling_temp.log && sed -i 's/--/---/g' basecalling_temp.log && diff basecalling.log basecalling_temp.log | grep pod | awk '{print $NF}' > reads_causing_errors.txt && while read -r line; do sed -n "/$line/,/---/ p" basecalling.log >> error.log; done < reads_causing_errors.txt && cat basecalling_temp.log | grep pod | awk '{print $NF}' | sed 's/pod5/bam/g' | sort -t_ -k4,4n > completed_bams.txt && rm basecalling_temp.log && echo "total pod5 files: $(grep pod5 basecalling.log | wc -l)" > basecalling_report.txt && echo "pod5 files with error: $(cat reads_causing_errors.txt | wc -l)" >> basecalling_report.txt && echo "successfully basecalled bam files: $(cat completed_bams.txt | wc -l)" >> basecalling_report.txt
```

---

# Concatenate Script

Next, we want to combine all the bam files from the same sample. If you work with large data and have many bam files, your system might give an error that looks like [this](https://www.biostars.org/p/10105/#31578):
```
[E::hts_open_format] Failed to open file "./PAS38658_8b84a00a_15c3d01c_245.bam" : Too many open files
samtools cat: fail to open file './PAS38658_8b84a00a_15c3d01c_245.bam': Too many open files
```
To overcome this, the following code will:
1. Take the `completed_bams.txt` and split it into files with 500 bam file names (`completed_bams_1.txt`, `completed_bams_2.txt`, …).
2. Add the corresponding path to each bam file in each split file with bam file names.
3. Use each split file as list of input for `samtools cat` and create larger bam files (for example, a sample with 1850 total bam files will create 4 split lists and therefore 4 bam files at this step).
4. Concatenate the larger bam files to create one bam file for the sample.

**NOTE**: This code can be run as a script or as a one liner, just make sure to change the variables accordingly. Run it where you’d like the final bam file to be saved.

```bash
sample="sample"
t=12

split -l 500 -d completed_bams.txt completed_bams_ --suffix-length=1 --additional-suffix=.txt

for cat_files in completed_bams_*.txt; do while read -r line; do path=$(find ./ -type f -name "$line") && sed -i "s|$line|$path|g" $cat_files; done < $cat_files; done

for bam_file_list in completed_bams_*.txt; do echo "Concatenating bams from $bam_file_list" && n=$(echo $bam_file_list | sed -n 's/.*_\([0-9]\+\)\.txt/\1/p') && samtools cat --threads ${t} -b $bam_file_list -o ${sample}_${n}.bam; done && rm completed_bams_*.txt

samtools cat --threads ${t} -o ${sample}.bam ${sample}_*.bam && rm ${sample}_*.bam
```

Alternatively, here it is as a one liner:

```bash
sample="sample" && t=12

split -l 500 -d completed_bams.txt completed_bams_ --suffix-length=1 --additional-suffix=.txt && for cat_files in completed_bams_*.txt; do while read -r line; do path=$(find ./ -type f -name "$line") && sed -i "s|$line|$path|g" $cat_files; done < $cat_files; done && for bam_file_list in completed_bams_*.txt; do echo "Concatenating bams from $bam_file_list" && n=$(echo $bam_file_list | sed -n 's/.*_\([0-9]\+\)\.txt/\1/p') && samtools cat --threads ${t} -b $bam_file_list -o ${sample}_${n}.bam; done && rm completed_bams_*.txt && samtools cat --threads ${t} -o ${sample}.bam ${sample}_*.bam && rm ${sample}_*.bam
```