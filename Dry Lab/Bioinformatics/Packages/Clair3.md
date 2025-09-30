---
aliases: 
tags: 
type: package
status: maintained
---

Paper: [Symphonizing pileup and full-alignment for deep learning-based long-read variant calling](https://www.nature.com/articles/s43588-022-00387-x.pdf)

GitHub: [HKU-BAL](https://github.com/HKU-BAL)/[Clair3](https://github.com/HKU-BAL/Clair3)

Documentation: []()

---

**Description**:  “Clair3 is a germline small variant caller for long-reads.”

---

# Installation

[Option 4. Build an anaconda virtual environment](https://github.com/HKU-BAL/Clair3?tab=readme-ov-file#option-4-build-an-anaconda-virtual-environment) worked on local lab server.

```bash
conda create -c conda-forge -c bioconda -n clair3_env python tensorflow=2.15.0 whatshap samtools parallel -y
source activate clair3_env

# install other packages to run C in environment
conda install -c conda-forge xz zlib bzip2 automake curl pigz cffi make gcc gxx -y 

# clone Clair3 and compile longphase and cffi library for c implement
cd ~/packages/
git clone https://github.com/HKU-BAL/Clair3.git
cd Clair3

make PREFIX=${CONDA_PREFIX}

# install pypy
cd ${CONDA_PREFIX}/bin # ~/miniconda3/envs/clair3_env/bin
wget https://downloads.python.org/pypy/pypy3.10-v7.3.19-linux64.tar.bz2 && tar -jxvf ${CONDA_PREFIX}/bin/pypy3.10-v7.3.19-linux64.tar.bz2
ln -sf pypy3.10-v7.3.19-linux64/bin/pypy3 ${CONDA_PREFIX}/bin/pypy3

# download pre-trained models
cd ~/packages/Clair3 && mkdir models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz && tar -zxvf clair3_models.tar.gz -C ./models
```

# Usage

- [ONT Variant Calling Quick Demo](https://github.com/HKU-BAL/Clair3/blob/main/docs/quick_demo/ont_quick_demo.md#ont-variant-calling-quick-demo)

```bash
conda activate clair3_env

modeln='/home/user/packages/Clair3/models/r1041_e82_400bps_sup_v500'
bam='/nanopore/data/merged_bams/sample.bam'
ref='/home/user/references/chm13_v2.0-t2t.fa'
output='/nanopore/data/phasing/'

~/packages/Clair3/run_clair3.sh --bam_fn=${bam} --ref_fn=${ref} --threads=12 --platform="ont" --model_path=${modeln} --use_whatshap_for_final_output_phasing --enable_long_indel --output=${output}
```