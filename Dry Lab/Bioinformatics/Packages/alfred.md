---
aliases: 
tags:
  - statistics
  - counts
  - alignment
type: package
status: maintained
---

Paper: [Alfred: interactive multi-sample BAM alignment statistics, feature counting and feature annotation for long- and short-read sequencing](https://academic.oup.com/bioinformatics/article/35/14/2489/5232224)

GitHub: [Alfred: BAM alignment statistics, feature counting and feature annotation](https://github.com/tobiasrausch/alfred)

Documentation: [Alfred: BAM Statistics, Feature Counting and Feature Annotation](https://www.gear-genomics.com/docs/alfred/)

---

# Installation

```bash
cd ~/packages
git clone --recursive https://github.com/tobiasrausch/alfred.git 
cd alfred/ 
make all 
make install

# symlink to base environment
ln -s ~/packages/alfred/bin/alfred ~/miniconda3/bin/alfred
```