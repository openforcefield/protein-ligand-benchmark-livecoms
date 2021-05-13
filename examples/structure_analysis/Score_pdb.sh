#!/bin/bash
# Prepare pdb files using OpenEye Spruce
# """
#
# Pipeline Step: Prepare raw PDB protein structures and score them.
#
# Description: This script pulls pdb files from a list and scores them using the Iridium criteria.
#
# Step Flowchart:
#     * Generate a text file with a pdb id per line for the structures to be scored.
#     * The inputs are the name of the pdb id text file and the Target name.
#     * Select which quality levels to structures to use:
#          * HT: Highly Trustworthy structures - a good model for the data and little or no missing density
#          * ALL: This include HT and MT (Moderarely Trustworthy) but no NT (Not Trustworthy) structures
#     * Output is a .csv file named target_HT/ALL_scoring_list.csv
# """

# This script requires installation of the
# OpenEye application Spruce, python OEChem toolkit and docopt python package

# Test for input command line arguments
if [ $# -ne 3 ]; then
    echo "Please supply: <pdb list file name> <target name> <Score Type (HT ALL)>"
    exit
fi
# Constants
zero=0
DPI_limit=0.75
LaD_limit=1.1
ASaD_limit=0.49
Iridium_score=0
Final_Iridium_score=1000000
HT_flag=0
MT_flag=0
target_name_count=$(tr -dc '_' <<<"$2" | wc -c)
(( final=target_name_count+1 ))

# Test for pdb data directory
# Create a pdb_data directory
# for downloading data into
if [ ! -d pdb_data ]
then
    mkdir pdb_data
    cd pdb_data || exit
else
    cd pdb_data || exit
    mkdir backup
    mv -- *.* backup/
fi

for file in $(cat ../"$1")
do
    getstructure "$file"
done
cd ../

# Test if there is a prep directory
# Create if not
if [ ! -d prep ]
then
    mkdir prep
    cd prep || return
else
    cd prep || return
    mkdir backup
    mv -- *.* backup/
fi

# Run spruce once to assess the quality of the structures
# Choosing to use .cif files versus pdb
# cif should have ligand bond order information while pdb does not
for file in ../pdb_data/*.cif
do
    tmp=$(echo "$file" | cut -d / -f 3)
    pdb_prefix=$(echo "$tmp" | cut -d . -f 1)
    map_prefix=$(echo "$file" | cut -d . -f 1-3)
    for list in "$map_prefix".m*
    do
        if [ -e "$list" ]
        then
            spruce -in "$file" -map "$list" -prefix "$2"-"$pdb_prefix"
        else
            spruce -in "$file" -prefix "$2"-"$pdb_prefix"
        fi
    done
done

# Determine the best structure based on
# ligand density, active site density, and DPI
# for this score LOWER is BETTER

selection="$3"

# print title for csv output file
if [[ "$selection" == *HT* ]]
then
    if [ -e ../"$2"\_HT_scoring_list.csv ]; then rm ../"$2"\_HT_scoring_list.csv; fi
    echo "PDB id,Ligand,Iridium Category,Ligand density,Active-site density,DPI,Iridium score,Packing residues" >> ../"$2"\_HT_scoring_list.csv
elif [[ "$selection" == *A[lL][lL]* ]]
then
    if [ -e ../"$2"\_ALL_scoring_list.csv ]; then rm ../"$2"\_ALL_scoring_list.csv; fi
    echo "PDB id,Ligand,Iridium Category,Ligand density,Active-site density,DPI,Iridium score,Packing residues" >> ../"$2"\_ALL_scoring_list.csv
else
    echo "Please select either HT or ALL (case sensitive)"
    exit
fi

HT_test=$(grep "Iridium" *output.log | grep -c "HT")
MT_test=$(grep "Iridium" *output.log | grep -c "MT")

if [ "$HT_test" -gt 0 ] || [ "$MT_test" -gt 0 ]
then
    for file in $(ls -1 *_output.log | cut -d \_ -f 1-"$final")
    do
        category=$(grep "Iridium" "$file"\_output.log | cut -d : -f 3 | cut -d , -f 1)
        LaD=$(grep "Iridium" "$file"\_output.log | cut -d : -f 4 | cut -d , -f 1)
        ASaD=$(grep "Iridium" "$file"\_output.log | cut -d : -f 5 | cut -d , -f 1)
        DPI=$(grep "Iridium" "$file"\_output.log | cut -d : -f 6 | cut -d , -f 1)
        file_name=$(grep "Iridium" "$file"\_output.log | cut -d : -f 2 | cut -d \> -f 1)
        ligand_name=$(grep "Iridium" "$file"\_output.log | cut -d : -f 2 | cut -d \> -f 2 | cut -d . -f 1 | cut -d , -f 1)
        packing_name=$(grep "Iridium" "$file"\_output.log | cut -d : -f 10 | cut -d , -f 1)
        h=0
        i=0
        j=0
        k=0
        l=0
        m=0
        n=0
        o=0
        C_tmp=""
        L_tmp="-1"
        A_tmp="-1"
        D_tmp="100"
        F_tmp=""
        Lig_tmp=""
        #P_tmp=""
        for job in ${LaD[*]}
        do
            if (( $(echo "$job > $L_tmp" |bc -l) )) && (( $(echo "$job <= 1.0" | bc -l) ))
            then
                L_tmp="$job"
                (( i=i+1 ))
            fi
        done
        for job in ${ASaD[*]}
        do
            if (( $(echo "$job > $A_tmp" |bc -l) )) && (( $(echo "$job <= 1.0" | bc -l) ))
            then
                A_tmp="$job"
                (( j=j+1 ))
            fi
        done
        for job in ${DPI[*]}
        do
            if (( $(echo "$job < $D_tmp" |bc -l) ))
            then
                D_tmp="$job"
                (( k=k+1 ))
            fi
        done
        for job in ${file_name[*]}
        do
            (( l=l+1 ))
            #echo $job $l $i $j $k
            if [ "$l" -eq "$i" ] && [ "$l" -eq  "$j" ] && [ "$l" -eq  "$k" ]
            then
                F_tmp="$job"
                n="$l"
            elif [ "$l" -eq  "$i" ] &&  [ "$l" -eq  "$k" ]
            then
                F_tmp="$job"
                n="$l"
            elif [ "$l" -eq  "$i" ] && [ "$l" -eq  "$j" ]
            then
                F_tmp="$job"
                n="$l"
            elif [ "$l" -eq  "$i" ]
            then
                F_tmp="$job"
                n="$l"
           fi
        done
        for job in ${category[*]}
        do
            (( h=h+1 ))
            #echo $i $j $k $job
            if [ "$h" -eq "$n" ]
            then
                C_tmp="$job"
            fi
        done
        for job in ${ligand_name[*]}
        do
            (( m=m+1 ))
            if [ "$m" -eq "$n" ]
            then
                Lig_tmp="$job"
            fi
        done
        for job in ${packing_name[*]}
        do
            (( o=o+1 ))
            if [ "$o" -eq "$n" ]
            then
                packing="$job"
            fi
        done
        LaD="$L_tmp"
        ASaD="$A_tmp"
        if (( $(echo "$D_tmp < $DPI_limit" | bc -l) ))
        then
            DPI="$D_tmp"
        else
            factor=$(echo 0.75 / "$D_tmp" | bc -l)
            factor=$(echo "$factor" / 101 | bc -l)
            DPI=$(echo "$DPI_limit" - "$factor" | bc -l)
            #echo $DPI
        fi
        LaD=$(echo "$LaD_limit" - "$LaD" | bc -l)
        ASaD=$(echo "$ASaD" - "$ASaD_limit" | bc -l)
        DPI=$(echo "$DPI_limit" - "$DPI" | bc -l)
        Numerator=$(echo "$LaD" / "$ASaD" | bc -l)
        Iridium_score=$(echo "$Numerator" / "$DPI" | bc -l)
        Iridium_cat="$C_tmp"

        if [[ "$selection" == HT && "$Iridium_cat" == *HT* ]]
        then
            echo "$F_tmp","$Lig_tmp","$Iridium_cat","$L_tmp","$A_tmp","$D_tmp","$Iridium_score","$packing" >> ../"$2"\_HT_scoring_list.csv
            if (( $(echo "$Iridium_score < $Final_Iridium_score" | bc -l) ))
            then
                Final_Iridium_score="$Iridium_score"
                Best_structure="$F_tmp"
                Best_ligand="$Lig_tmp"
            fi
            HT_flag=1
        fi

        if [[ "$selection" == A[lL][lL]  && "$Iridium_cat" == *[HM]T* ]]
        then
            echo "$F_tmp","$Lig_tmp","$Iridium_cat","$L_tmp","$A_tmp","$D_tmp","$Iridium_score","$packing" >> ../"$2"\_ALL_scoring_list.csv
            if (( $(echo "$Iridium_score < $Final_Iridium_score" | bc -l) ))
            then
                Final_Iridium_score="$Iridium_score"
                Best_structure="$F_tmp"
                Best_ligand="$Lig_tmp"
            fi
            MT_flag=1
        fi
    done
else
    echo "No good structures in this set"
    exit
fi

if [[ "$selection" == HT && "$HT_flag" -eq 1 ]]
then
    echo ""
    echo "The best structure is " "$Best_structure" "$Best_ligand" "$Final_Iridium_score" >> ../"$2"\_HT_scoring_list.csv
    uniq ../"$2"\_HT_scoring_list.csv tmp.csv
    mv tmp.csv ../"$2"\_HT_scoring_list.csv
elif [[ "$selection" == A[lL][lL] && "$MT_flag" -eq 1 ]]
then
    echo ""
    echo "The best structure is " "$Best_structure" "$Best_ligand" "$Final_Iridium_score" >> ../"$2"\_ALL_scoring_list.csv
    uniq ../"$2"\_ALL_scoring_list.csv tmp.csv
    mv tmp.csv ../"$2"\_ALL_scoring_list.csv
else
    if [[ "$HT_flag" -ne 1 && "$selection" == HT ]]
    then
        rm ../"$2"\_HT_scoring_list.csv
        echo "There are no HT structures in this data set"
    fi
    if [[ "$MT_flag" -ne 1 && "$selection" == A[lL][lL] ]]
    then
        rm ../"$2"\_ALL_scoring_list.csv
        echo "There are no HT/MT structures in this data set"
    fi
fi

cd ../




