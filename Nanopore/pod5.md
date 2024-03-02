---
aliases: 
tags: nanopore
type: package
status: maintained
---

GitHub: [nanoporetech](https://github.com/nanoporetech)/[pod5-file-format](https://github.com/nanoporetech/pod5-file-format)

Documentation: [Pod5 File Format Documentation](https://pod5-file-format.readthedocs.io/en/latest/index.html)

---

**Description**: File format for storing raw [[nanopore]] data.

---

# Install

```bash
pip install pod5
```

# Usage

[Converting fast5 files into pod5 files](https://github.com/nanoporetech/pod5-file-format/blob/master/python/pod5/README.md#pod5-convert-fast5) #convert 

```bash
# convert each fast5 into its corresponding pod5
pod5 convert fast5 --threads 12 ./fast5_pass/*.fast5 --output pod5_pass/ --one-to-one ./fast5_pass/ 2> pod5.log

# convert all fast5 files into one pod5 file
pod5 convert fast5 --threads 12 --recursive ./ --output converted.pod5 2> pod5.log

# covert fast5 files in specific fast5 folders into one pod5 file
pod5 convert fast5 --threads 12 /path/to/folder1/fast5 /path/to/folder2/fast5 /path/to/folder3/fast5 --output converted.pod5 2> pod5.log
```

# Debug

```bash
# add POD5_DEBUG=1 infront command, example:
POD5_DEBUG=1 pod5 view ${file}.pod5
```

# Errors

## `pod5 view`

```
POD5 has encountered an error: 'Error while processing 'igvf004_13A-gc_lig-ss_p2_1.pod5''

For detailed information set POD5_DEBUG=1'
```

### `polars.exceptions.ColumnNotFoundError: not_set`

According to the post: [pod5 view fails with polars.exceptions.ColumnNotFoundError: not_set #104](https://github.com/nanoporetech/pod5-file-format/issues/104), “Please install polars 0.19”
```bash
pip install polars==0.19
```

## `pod5 subset`

### `'enable_string_cache() missing 1 required positional argument: 'enable''`
```
POD5 has encountered an error: 'enable_string_cache() missing 1 required positional argument: 'enable''

For detailed information set POD5_DEBUG=1'
```

According to the post: # [pod5 filter/ subset failed with TypeError: enable_string_cache() missing 1 required positional argument: 'enable' #105](https://github.com/nanoporetech/pod5-file-format/issues/105)
```bash
pip install polars==0.19.7
```