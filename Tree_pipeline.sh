#!/bin/bash
unset -f module
. /tgac/software/testing/lmod/6.1/x86_64/lmod/lmod/init/profile
ml biopython


$OUTPUT="PST130_concatinatelan.phy"

# Create a directory on the working dir where the script was started from
mkdir sorted

# Find all the Isolates fasta files and order them by gene
FILES="../Isolates2015/*.fa"
shopt -s nullglob
for f in $FILES;
do
    python ./sort_fasta.py -f $f > ./sorted/$(basename $f).sorted
done

# With all the files ordered by gene and placed in the sorted folder, filter out which genes are relevant (have enough information)
# codon_from_fasta has multiple filtering parameters, the most relevant are -m (minimum percentage of known bases in a sequence for acceptance) -s (minimum number of accepted samples percentage)
python ./codon_from_fasta.py -d sorted

# Select all the gene filtered sequences
pattern="./sorted/*.filtered"
files=( $pattern )

# Write the PHYLIP header to the final.phy file
echo -n ${#files[@]} > $OUTPUT
echo -n " " >> $OUTPUT

wc -m < ${files[0]} >> $OUTPUT

# Write all the sequences to the final.phy file to generate a sequential PHYLIP file
for filename in sorted/*.filtered; do echo -n "$(basename $filename | cut -d '_' -f 1) "; cat $filename; echo ""; done >> $OUTPUT

# The final.phy file is ready to go into RAxML
