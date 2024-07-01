internal link: [[Nanopore]]

---

This script works up a Nanopore run with barcoded reads after sequencing.

> [!warning] 
> This code is written for a folder structure produced by live basecalling and `pod5` files. Basecalling in the terminal with `dorado basecaller` might produce a folder structure that is slightly different.

This script will:
1.  Get all of the run names by finding all the run folders
2. For each run:
	1. Get the parent directory and `fastq_pass` path
	2. Create a log file
	3. Count the total number of raw reads
	4. Look for `pod5_skip` folder, and count number of skipped raw reads if found
	5. For each barcode: 
		1. Concatenate all `fastq` files
		2. Count total number of reads in `fastq`
		3. Compress `fastq`
		4. Create checksum for  `fastq`
	6. Repeat Step 5 for unclassified reads

**NOTE**: Run this script in the experiment folder.

```bash
# get run folder names in experiment folder
runs=($(find . -maxdepth 1 -type d | grep -v '^\.$' | cut -d "/" -f 2))

# for each run
for run in "${runs[@]}"; do
	# get folder names 
	f=`find ./$run -type f -name "*pass*_0.fastq.gz" | head -1`
	pdir=`echo $f | cut -d'/' -f 1,2,3`
	pass_dir=`echo $f | cut -d'/' -f 1,2,3,4`

	# create a log file
	echo -e "===\n$run\n===\n" && echo -e "===\n$run\n===\n---\n" > ${pdir}/${run}.log

	# count number of raw reads
	echo -e "counting total raw reads\n" && echo -e "total raw reads: $(pod5 view --recursive ${pdir} | grep -v read_id | wc -l)\n" >> ${pdir}/${run}.log

	# look for skipped folder and count the number of skipped raw reads
	if [[ -d "${pdir}/pod5_skip" ]]; then
		echo -e "counting skipped raw reads\n" && echo -e "skipped raw reads: $(pod5 view ${pdir}/pod5_skip | grep -v read_id | wc -l)\n" >> ${pdir}/${run}.log
	fi

	# process each barcode folder
	for i in {01..24}; do
		echo -e "barcode $i\n---" && echo -e "barcode $i\n---" >> ${pdir}/${run}.log
		echo "concatenating reads" && zcat ${pass_dir}/barcode${i}/*.fastq.gz > ${pass_dir}/${run}_barcode${i}.fastq
		echo "counting reads" && echo "number of reads: $(echo $(cat ${pass_dir}/${run}_barcode${i}.fastq|wc -l)/4|bc)" >> ${pdir}/${run}.log
		echo "compressing fastq" && pigz -p 12 ${pass_dir}/${run}_barcode${i}.fastq
		echo "checking file size" && echo -e "file size: $(du -sh ${pass_dir}/${run}_barcode${i}.fastq.gz | cut -d $'\t' -f 1)" >> ${pdir}/${run}.log
		echo -e "creating checksum\n" && echo -e "checksum: $(md5sum ${pass_dir}/${run}_barcode${i}.fastq.gz | cut -d ' ' -f 1)\n" >> ${pdir}/${run}.log
	done

	# process the unclassified folder
	echo -e "unclassified\n---" && echo -e "unclassified\n---" >> ${pdir}/${run}.log
	echo "concatenating reads" && zcat ${pass_dir}/unclassified/*.fastq.gz > ${pass_dir}/${run}_unclassified.fastq
	echo "counting reads" && echo "number of reads: $(echo $(cat ${pass_dir}/${run}_unclassified.fastq|wc -l)/4|bc)" >> ${pdir}/${run}.log
	echo "compressing fastq" && pigz -p 12 ${pass_dir}/${run}_unclassified.fastq
	echo "checking file size" && echo -e "file size: $(du -sh ${pass_dir}/${run}_unclassified.fastq.gz | cut -d $'\t' -f 1)" >> ${pdir}/${run}.log
	echo -e "creating checksum\n" && echo -e "checksum: $(md5sum ${pass_dir}/${run}_unclassified.fastq.gz | cut -d ' ' -f 1)\n" >> ${pdir}/${run}.log
done
```