internal link: [[Nanopore]]

---

One of our samples kept [aborting during basecalling](https://github.com/nanoporetech/dorado/issues/548), and a suggestion was to [loop through each file the pod5 folder and basecall each individually](https://github.com/nanoporetech/dorado/issues/548#issuecomment-1888051514). Here is the script I used to do that. 

NOTE: This example creates bams as output. For fastqs, replace all instances of `bam` with `fastq` and edit `dorado_options`  accordingly.

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

# Script

If the variables look correct, edit the dorado variables and run this script in the run folder. This script will 
1. Find all the pod5 files and create all the variables (the same as the code above).
2. Create the bam output folder if it doesn’t exist.
3. If the output bam file does not exist, the basecall the corresponding pod5 file with `dorado`
4. Save dorado progess messages in a log file called `basecalling.log`
5. Remove lines in the `basecalling.log` file that contain the string `info` WITHOUT the string `Finished`. The will help keep the file from getting too big (especially for larger runs). All other messages like `error` and `warning` will be saved. 

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
		echo "Processing $pod5_name" 2>&1 | tee -a basecalling.log
		# basecall pod5 file with dorado and save progess statments in log file
		${dorado_command} ${dorado_config} "${pod5_file}" ${dorado_options} > "${output_bam}" 2>> >(tee -a basecalling.log >&2)
		# clean-up log file to remove extra information
		sed -i '/info/ { /Finished/!d; }' basecalling.log
	fi
done

```

Alternatively, here it is as a one liner:

```bash
dorado_command="$HOME/packages/dorado/dorado-0.5.1-linux-x64/bin/dorado basecaller"; dorado_options="--modified-bases 5mCG_5hmCG --device cuda:0"; dorado_config="$HOME/packages/dorado/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0"; for pod5_file in $(find ./ -type f -name '*.pod5'); do pod5_name="$(basename "${pod5_file}")"; bam_file="$(basename "${pod5_file%.pod5}.bam")"; output_dir="$(dirname "${pod5_file}" | sed 's/pod5/bams/')"; output_bam="${output_dir}/${bam_file}"; [ ! -d "${output_dir}" ] && mkdir -p "${output_dir}"; [ ! -e "${output_bam}" ] && { echo "Processing $pod5_name" 2>&1 | tee -a basecalling.log; ${dorado_command} ${dorado_config} "${pod5_file}" ${dorado_options} > "${output_bam}" 2>> >(tee -a basecalling.log >&2); sed -i '/info/ { /Finished/!d; }' basecalling.log; }; done
```