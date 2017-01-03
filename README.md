# q2-deblur

Wrapper for Deblur

# Install

```
conda install -c bioconda VSEARCH MAFFT SortMeRNA==2.0 biom-format
pip install deblur
```

# Example

```
$ qiime deblur workflow --i-demultiplexed-seqs data.qza --o-table deblurred_table --o-representative-sequences deblurred_sequences
```
