#!/bin/bash

# Test for input command line arguments
if [ $# -ne 5 ]; then
    echo "Please supply: <directory path> <target name> <quality> <TC cutoff> <irid cutoff>"
    exit
fi



for lig in $(cat "$1"/prep/"$3"/"$2"\_rocs_"$3"\_score_1.csv)
do 
    lig_id=$(echo "$lig" | cut -c 1-3)
    tc=$(echo "$lig" | cut -d , -f 4)
    pdb_id=$(grep "$lig_id" "$1"/"$2"\_"$3"\_score_stripped.csv | cut -c 1-4)
    irid_score=$(grep "$lig_id" "$1"/"$2"\_"$3"\_score_stripped.csv | cut -d , -f 7)
    if (( $(echo "$tc > $4" | bc -l) )) && (( $(echo "$irid_score < $5" | bc -l) ))
        then echo "$pdb_id" "$lig_id" "$irid_score" "$tc"
    fi
done
