# Installation 

This [page](http://hgdownload.soe.ucsc.edu/admin/exe/) has installation instructions and link
```bash
# http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
# save to your base conda env
cd ~/miniconda/bin
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bigBedToBed ./
```

# Usage
```bash
bigBedToBed input.bb output.bed
```