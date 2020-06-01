# NanoSplicer
NanoSplicer is a program for accurately identify splice sites using Nanopore data (both basecalls and squiggles).
[Documentation](https://youyupei.github.io/NanoSplicer)

# Abstract
Alternative splicing is an essential mechanism that enables a single gene to produce multiple mRNA products (called isoforms). Oxford Nanopore sequencing produces long reads that have natural advantages for characterising isoforms. Alternative splicing can not only add or skip entire exons, but can also vary exon boundaries by selecting different nucleotides as splice sites. However, accurately identifying the latter is challenging for nanopore sequencing, as exon boundaries are often unclear due to the high error rate. One existing solution is polishing nanopore reads with short reads. While feasible, this approach requires both short and long reads, which adds considerable expensive. Furthermore, isoform distinguishing short reads are not always available (e.g. in 10X scRNAseq). Therefore, a method that could accurately identify exon boundaries solely from nanopore reads would have numerous advantages. We are developing a method to characterise exon boundaries using only nanopore sequencing data. Nanopore sequencing records changes in electrical current when a DNA or RNA strand is traversing through a pore. This raw signal (known as a squiggle) is then basecalled by machine learning algorithms, but this process is error-prone. Instead of looking at the basecalled sequences only, we are also using the squiggle to characterise the exon boundaries. We tested our method using synthetic mRNAs with known exon boundaries, demonstrating our initial implementation can correctly assign % of the squiggles overlapping splice junctions to the correct exon boundaries. We conclude using raw squiggle data is a promising approach for accurately identifying exon boundaries in nanopore sequencing reads.

## Keywords:
Oxford Nanopore sequencing, Transcriptomics

## Contributing authors:
Yupei You,
Mike Clark*,
Heejung Shim*

## Author affiliations:
School of Mathematics and Statistics and Melbourne Integrative Genomics, University of Melbourne, Melbourne, Australia
Anatomy and Neuroscience, School of Biomedical Sciences, University of Melbourne, Melbourne, Australia
School of Mathematics and Statistics and Melbourne Integrative Genomics, University of Melbourne, Melbourne, Australia

