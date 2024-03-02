---
aliases: 
tags: nanopore
type: package
status: maintained
---

GitHub: [nanoporetech](https://github.com/nanoporetech)/**[dorado](https://github.com/nanoporetech/dorado)**

---

**Description**: Open source basecaller for [[Nanopore|Oxford Nanopore]] reads.

---

# Install

```bash
cd ~/packages/dorado
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.0-linux-x64.tar.gz
tar -xvzf dorado-0.5.0-linux-x64.tar.gz
cd ~/miniconda3/bin
ln -s ~/packages/dorado/dorado-0.5.0-linux-x64/bin/dorado dorado
```

Download [Basecalling models](https://github.com/nanoporetech/dorado#dna-models)

```bash
cd ~/packages/dorado/dorado_models
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.3.0
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.2.0
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.1.0

dorado download --model rna004_130bps_sup@v3.0.1
```

# Usage

The basic syntax for using `dorado` is the following:
```bash
dorado subcommand model raw_data options > output
```

Depending on the raw data, it will more specifically look something like this:
```bash
# point to raw data file
dorado subcommand model sample.pod5 options > output

# point to raw data folder
dorado subcommand model ./pod5 options > output

# look recursively and find raw data
dorado subcommand model ./ --recursive other_options > output
```

## Examples
```bash
# cDNA
dorado basecaller ~/packages/dorado/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 ./pod5 --device "cuda:0" --min-qscore 10 --emit-fastq > ${file}.fastq

# genomic
# with mapping 
dorado basecaller ~/packages/dorado/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 pod5/ --modified-bases 5mCG_5hmCG --reference ~/references/mm39.fa --secondary no > ${file}_unsorted.bam
# withOUT mapping 
dorado basecaller ~/packages/dorado/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 pod5/ --modified-bases 5mCG_5hmCG > ${file}_unmapped.bam

# barcoding
dorado basecaller ~/packages/dorado/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 ./pod5 --kit-name SQK-NBD114-24 --min-qscore 10 > ${file}.bam
```

---

# Flags

[From post:](https://github.com/nanoporetech/dorado/issues/110) The basecaller can `--emit-moves` which gives you the `mv:B:c` tag in the output which contains a coarse sequence to signal mapping.

[Examples](https://github.com/nanoporetech/dorado/issues/284) for using `--resume-from`.

---

# Errors

## Sample rate not compatible 

If you choose the model with the incorrect sampling rate, you will get an [error like this](https://github.com/nanoporetech/dorado/issues/357):
```
[2023-09-19 10:42:40.122] [error] Sample rate for model (5000) and data (4000) are not compatible.
```

Then you need to use `dna_r10.4.1_e8.2_400bps_sup@v4.1.0` instead of `dna_r10.4.1_e8.2_400bps_sup@v4.2.0`

## Cannot create index
After using dorado basecaller, if you try to index the bam file and you get [this error](https://bioinformatics.stackexchange.com/questions/14179/no-coor-reads-not-in-a-single-block-at-the-end-0-1) :

```
[E::hts_idx_push] NO_COOR reads not in a single block at the end 21 -1
[E::sam_index] Read '0633b8c3-5fb9-4966-9aa3-3fb89ddb8df8' with ref_name='chr5', ref_length=151758149, flags=16, pos=117768162 cannot be indexed
samtools index: failed to create index for "igvfm_260-21_lig-dna_p2_1_unsorted.bam"
```

It means that the reads are not sorted by coordinate and you will need to sort with [[samtools]].
