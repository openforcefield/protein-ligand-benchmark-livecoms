#!/bin/bash


# Test for input command line arguments
if [ $# -ne 4 ]; then
    echo "Please supply: <directory path> <target name> <quality> <irid/dpi cutoff>"
    exit
fi

if [ ! -d "$1"/prep/"$3" ]
then
    mkdir "$1"/prep/"$3"
fi

if [ ! -d "$1"/prep/"$3"/ligands ]
then
    mkdir "$1"/prep/"$3"/ligands
fi

# figure out what directory the files are in
#layer_count=$(tr -dc '/' <<<"$1" | wc -c)
#(( layer_count=layer_count+2 ))
#(( directory_count=layer_count-1 ))
#echo "$layer_count" "$directory_count"


for file in $(cat "$1"/"$2"\_"$3"\_score_stripped.csv)
do 
    test=$(echo "$file" | cut -d , -f 3)
    pdb_id=$(echo "$file" | cut -c 1-4)
    dpi=$(echo "$file" | cut -d , -f 6)
    irid_score=$(echo "$file" | cut -d , -f 7)
    pack=$(echo "$file" | cut -d , -f 8)
    if [[ "$4" == ALL ]]
    then
        if (( $(echo "$irid_score < $4" | bc -l) )) && [[ "$pack" == *false* ]]
        then
            cp "$1"/prep/"$pdb_id"*__ligand__*.oeb "$1"/prep/"$3"/ligands/
        fi
    else
       cp "$1"/prep/"$pdb_id"*__ligand__*.oeb "$1"/prep/"$3"/ligands/
    fi
done
