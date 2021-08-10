# BED12 files in NanoSplicer output for easier visualisation
Except for the assignment probability table. NanoSplicer output 2 BED12 files in addition: 

## BED version of identification for each JWR)
A BED file (`NanoSplicer output name].bed`) will be output: each BED entry contains the identification for a single JWR. JWRs identified with low assignment probability (Default: less than 0.8) are presented in `grey` while the other JWRs are in `black`.The 'score' (5th column) contains the probability of the identified splice junction being the true one. 

**Note:** JWRs with low Squiggle Information quality is excluded in this BED file. They usually come from the reads from relatively low-quality squiggle and are error-prone.



## Number of support for each splice junction
A BED file `junction_suppot.bed` will also be output. Different from the previous BED file. In `junction_suppot.bed`, each BED entry is a uniq splice junction with the `score` field indicates the number of JWRs are supporting it. All the JWRs with high probability included in the previous BED file will also be included here. Besides, any JWR that is not included because of the JAQ threshold in `JWR_subset` will also be counted.