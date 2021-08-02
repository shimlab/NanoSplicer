# get all the junction within read (JWR) from the BAM file
python3 ../JWR_checker.py --output_csv example.bam all_JWR.h5

# subselect the JWR with JAQ<=0.8
python3 ../JWR_subset.py --best_JAQ=0.8 --output_csv all_JWR.h5 JWR_subset.h5

# identify the splice junction for the subselected JWR
python3 ../NanoSplicer.py -i example.bam -f fast5s/ -r sequins_genome.fa JWR_subset.h5
