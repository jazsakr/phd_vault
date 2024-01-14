Sometimes I want to write a script that will manipulate a [[Nanopore]] run folder but I need a place to test it. One example is when [basecalling kept aborting for one of our samples](https://github.com/nanoporetech/dorado/issues/548) and a suggestion was to [basecall each [[pod5]] file separately](https://github.com/nanoporetech/dorado/issues/548#issuecomment-1888051514). I needed to write [a script](Loop%20Nanopore%20basecalling.md) to do that and I wanted folder to test it on.

Here is how I make a copy of an existing run to use as my code testing folder:

```bash
# go above directory you want to use to copy

mkdir ont_code_testing

source_dir="exp_name"
testing_dir="ont_code_testing"

# copy over all files except the big sized files 
rsync -av --progress --exclude={"*.pod5","*.bam","*.fastq","*.fastq.gz"} "./$source_dir/" "./$testing_dir/"

# copy over just 5 pod5 files for each pod5 folder
for pod5_dir in $(find "$source_dir" -type d -name 'pod5'); do pod5_testing_dir="$(echo $pod5_dir | sed "s|$source_dir|$testing_dir|")"; pod5_files=($(find "$pod5_dir" -type f -name '*.pod5' | head -5)); rsync -av --progress "${pod5_files[@]}" "$pod5_testing_dir"; done
```

Future direction would be to copy over some of the other file types that was excluded. Might look something like this:

```bash
source_dir="exp_name"
testing_dir="ont_code_testing"

file_types=("pod5","bam","fastq","fastq.gz")

for type in $file_types; do 
	# add other lines of code from block above
done
```