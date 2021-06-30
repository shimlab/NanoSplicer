# NanoSplicer
NanoSplicer is a program for accurately identify splice sites using Nanopore data (both basecalls and squiggles).


## Keywords:
Oxford Nanopore sequencing, Transcriptomics

# Abstract
Alternative splicing is an essential mechanism that enables a single gene to produce multiple mRNA products (called isoforms). Oxford Nanopore sequencing produces long reads that have natural advantages for characterising isoforms. Alternative splicing can not only add or skip entire exons, but can also vary exon boundaries by selecting different nucleotides as splice sites. However, accurately identifying the latter is challenging for nanopore sequencing, as exon boundaries are often unclear due to the high error rate. One existing solution is polishing nanopore reads with short reads. While feasible, this approach requires both short and long reads, which adds considerable expensive. Furthermore, isoform distinguishing short reads are not always available (e.g. in 10X scRNAseq). Therefore, a method that could accurately identify exon boundaries solely from nanopore reads would have numerous advantages. We are developing a method to characterise exon boundaries using only nanopore sequencing data. Nanopore sequencing records changes in electrical current when a DNA or RNA strand is traversing through a pore. This raw signal (known as a squiggle) is then basecalled by machine learning algorithms, but this process is error-prone. Instead of looking at the basecalled sequences only, we are also using the squiggle to characterise the exon boundaries. We tested our method using synthetic mRNAs with known exon boundaries, demonstrating our initial implementation can correctly assign % of the squiggles overlapping splice junctions to the correct exon boundaries. We conclude using raw squiggle data is a promising approach for accurately identifying exon boundaries in nanopore sequencing reads.

# Overview
`JWR_checker.py`: Find junctions within reads (JWRs) from a spliced mapping result (BAM).

`JWR_subset.py`: Subset the result from the JWR_checker. Based on our performance assessment, JWRs with good junction alignment quality (JAQ) have usually an accurate mapping to the splice junctions. A subset of JWR can be obtained by only select JWRs with low JAQ. By default, `JWR_subset.py` selects the all JWRs with JAQ <= 0.9.

`NanoSplicer.py`: Run the identifications on the `JWR_checker.py` (or `JWR_subset.py` if applied) output. 

# Requirements
NanoSplicer has been currently test on python 3.6 and 3.7. Everything should work for python3.X. 

## Package Dependency
For `JWR_checker.py` and `JWR_subset.py`:
* `pandas`
* `pysam`
* `numpy`
* `tqdm`
* `h5py`

Additional requirements for `NanoSplicer`:
* `tombo`
* `ont_fast5_api`
* `matplotlib`
* `fcntl`
* `intervaltree`
* `scipy`
* `skimage`

# Install

```
git clone https://github.com/shimlab/NanoSplicer.git
```

# Container (referred)
The environment required for running NanoSplicer has been package in container and can be accessed using `singularity`, which is supported by most hpc:

```
singularity pull NanoSplicer_container.sif docker://youyupei/nanosplicer:v1
```

## For poeple not familiar with container:
You can run linux command within the container by using `singularity shell` or `singularity exec`. These command autimatically bind your home directory to the home directory inside the container, which means everything under `~/` will be accessible without extra step (including the sub-directory). If your data are save in a different directory, you'll need to bind the directory with `-B <local path>:<path in container>` when running `singularity shell` or `singularity exec`. For example, your data are in `/data`, you need to add `-B /data:/folder_in_container`, everything in `/data` can then be accessible in `/folder_in_container`. You may use the same name for convenience (e.g. `-B /data:/data`). More formal introduction can be found at https://sylabs.io/singularity/





# Example 
```
python3 JWR_checker.py 
```









## Contributing authors:
Yupei You,
Mike Clark*,
Heejung Shim*

## Author affiliations:
School of Mathematics and Statistics and Melbourne Integrative Genomics, University of Melbourne, Melbourne, Australia
Anatomy and Neuroscience, School of Biomedical Sciences, University of Melbourne, Melbourne, Australia
School of Mathematics and Statistics and Melbourne Integrative Genomics, University of Melbourne, Melbourne, Australia

