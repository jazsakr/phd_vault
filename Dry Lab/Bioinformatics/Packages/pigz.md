---
aliases: 
tags: 
type: package
status: maintained
---

Homepage: [pigz - Parallel gzip](https://zlib.net/pigz/)

GitHub: [madler](https://github.com/madler)/[pigz](https://github.com/madler/pigz)

Documentation: [Manual](https://zlib.net/pigz/pigz.pdf)

---

[**Description**](https://zlib.net/pigz/): _pigz_, which stands for **p**arallel **i**mplementation of **gz**ip, is a fully functional replacement for gzip that exploits multiple processors and multiple cores to the hilt when compressing data.

---

# Install

```bash
conda install -c conda-forge pigz
```

# Usage

```bash
pigz -p 12 file.fastq

# compress multiple files
pigz -p 12 *.fastq
```