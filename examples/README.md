## Files:

1.**example.bam**: a bam file contain 10 alignments

2.**sequins_genome.fa**: the reference file used in the alignment

3.**fast5s/**: folder contains all the squiggles that corresponding to the reads in example.bam

4.**script.sh**: bash script containing the example commands. 

## Expected output of this example (after running `script.sh`):
1. Output from `JWR_checker.py`:
    **'all_JWR.h5'** and **'all_JWR.h5.csv'**: All JWRs extracted (43 in total) from the input BAM file ('example.bam')
2. Output from `JWR_subset.py`:
    **'JWR_subset.h5'** and **'JWR_subset.h5.csv'**: Subseletied JWRs (9 JWRs) with JAQ <= 0.8.
3. Output from `NanoSplicer.py`:
    **'output_prob_table.tsv'**: Probability table for 9 JWRs
    **'output_jwr.bed'**: NaoSplicer identification in BED12 format (7 JWRs)
    **'output_error_summary.csv'**: JWRs that failed to run by NanoSplicer (0 JWRs)
    *see https://github.com/shimlab/NanoSplicer#output for more detailed explaination of 'output_prob_table.tsv' and 'output_jwr.bed'
    


## Known issue:
1. There will be a performance warning when running JWR_checker and JWR_subset. The warning can be ignored. The warning is triggered when data is saved into the HDF5 file. It relates to python objects that can not be directly converted to c-type, causing non-optimal performance of the HDF5 file.
2. The progress bar in NanoSplicer.py is based on the number of fast5 files that have been processed. There is only 1 fast5 file in `example/fast5/`, the progress bar will be at 0% (0/1 completed) until the run has finished.
