internal link: [[Nanopore]]

---

When I write a bash script that works up a Nanopore run after sequencing, I don’t want to test it on a real run folder. But simply copying a run folder is not helpful because runs can be very large and testing code on them might take a while with large files.

Here is some code to create a code-testing run folder from a real run. 

> [!warning] 
> This code is written for live basecalled bulk or barcoded runs with `pod5` files and `fastq` file outputs, **NOT** `bam` outputs (yet).

This script will:

1. Create a new directory.
2. Copy over folders and files EXCEPT the `sequencing_summary` file (because it is a very large file) or files that end in `.pod5`, `.fastq`, or `.bam`
3. Copy over only the first 10 lines of the `sequencing_summary` file. 
4. Copy over just 3 `pod5` files for each pod5 folder
5. Copy over 100 reads from each `fastq` file
	-  Since runs with barcoded reads have a slightly different folder structure than bulk reads, it will check for barcode folders if it does not find files in the `fastq` folders

**NOTE**: Run this script above the run folder you want to use to copy and make sure to replace `source_dir` and `testing_dir` with your own folder names.

```bash
# replace with your own folder names
source_dir="exp23"
testing_dir="exp_testing"

# make the new directory for testing code
mkdir $testing_dir

# copy over folders and some files
rsync -av --progress --exclude={"*.pod5","*.bam","*.fastq","*.fastq.gz"} \
--exclude="*sequencing_summary_*txt" "./$source_dir/" "./$testing_dir/"

# copy over a few lines from sequencing_summary file
seq_sum=`find ./$source_dir -type f -name "sequencing_summary*txt"` && \
seq_sum_testing="$(echo $seq_sum | sed "s|$source_dir|$testing_dir|")" \
&& head seq_sum > seq_sum_testing

# copy over just 3 pod5 files for each pod5 folder
for pod5_dir in $(find "$source_dir" -type d -name 'pod5*'); do \
pod5_testing_dir="$(echo $pod5_dir | sed "s|$source_dir|$testing_dir|")"; \
pod5_files=($(find "$pod5_dir" -type f -name '*.pod5' | head -3)); \
rsync -av --progress "${pod5_files[@]}" "$pod5_testing_dir"; done

# copy over 100 reads in each fastq file
for fastq_dir in $(find "$source_dir" -type d -name 'fastq*'); do
	fastq_testing_dir="$(echo $fastq_dir | sed "s|$source_dir|$testing_dir|")"

	# Check for .fastq.gz files in the fastq folder
	fastq_files=($(find "$fastq_dir" -maxdepth 1 -type f -name '*.fastq.gz' | head -10))

	if [ ${#fastq_files[@]} -gt 0 ]; then
		# Process files directly in the fastq folder
		for file in "${fastq_files[@]}"; do
			echo "$file" | rev | cut -d'/' -f1 | rev
			testing_file="$fastq_testing_dir/$(basename "$file")"
			zcat "$file" | head -400 | gzip > "$testing_file"
		done
	else
		# Check for barcode folders in the fastq folder
		subdirs=($(find "$fastq_dir" -mindepth 1 -maxdepth 1 -type d))
		# Process files in each barcode folder
		for subdir in "${subdirs[@]}"; do
			sub_fastq_files=($(find "$subdir" -type f -name '*.fastq.gz' | head -10))
			for file in "${sub_fastq_files[@]}"; do
				echo "$file" | rev | cut -d'/' -f1 | rev
				sub_testing_dir="$(echo $subdir | sed "s|$source_dir|$testing_dir|")"
				testing_file="$sub_testing_dir/$(basename "$file")"
				zcat "$file" | head -400 | gzip > "$testing_file"
			done
		done
	fi
done
```