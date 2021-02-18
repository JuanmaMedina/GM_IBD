#!/bin/bash

# Usage example:
# bash seq_extractor_mwpi_calc.sh accesories_top1000_UHGG accesories_top1000_results/PF00871_2-7-2-15_acc.out propk_acc_seqs.faa propk_acc_ipc.txt propk_acc_ipc.png
# bash seq_extractor_mwpi_calc.sh full_pangenomes_4644_UHGG total_acetolactate_decarboxilase.out acetlac_seqs.faa acetlac_acc_ipc.txt acetlac_ipc.png

# 6. Retrieve sequences from .fasta headers captured with HMMER
# faSomeRecords, from Jim Kent (https://www.biostars.org/p/329774/)
# Make it executable: chmod a+x faSomeRecords

# KENT VERSION does not seem to work properly on server, as it retrieves all sequences as matches
# Use python2 version (wget https://raw.githubusercontent.com/santiagosnchez/faSomeRecords/master/faSomeRecords.py)


# Input defined variables:
GENOMIC_FOLDER=$1  # Accesory or core folder: e.g. accesories_top1000_UHGG
HEADERLIST=$2      # HMMER output with headers: e.g. accesories_top1000_results/PF00871_2-7-2-15_acc.out
SEQS=$3            # .fa with hit sequences: e.g. propk_seqs.faa

# Change Conda environment
echo "Activating Python 2..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py2

for SPECIES in $GENOMIC_FOLDER/genomic_MGYG-HGUT*_final.faa; do
#for SPECIES in $GENOMIC_FOLDER/*.faa; do

	echo "Extracting hits from genome $SPECIES..."

	#./faSomeRecords $SPECIES $HEADERLIST "$(basename "$SPECIES" .faa)_indseqs.faa"
	python faSomeRecords.py --fasta $SPECIES --list $HEADERLIST --outfile "$(basename "$SPECIES" .faa)_indseqs.faa" # We could use --stdout option, but keep like this to keep script integrity
done

echo "Merging individual-species hits and cleaning..."
cat *indseqs.faa > $SEQS
rm *indseqs.faa

# Flags of "faSomeRecords"
# SPECIES: individual species proteome in fasta
# HEADERLIST: header list file (remove ">" if present)
# 3: hit sequences from every species

# 6A. Calculate pI/mW and perform 2D EF (ipc.py)
# From http://isoelectric.org/
# Runs with python2, requires numpy, scipy and matplotlib to elaborate the 2D plots
# Follow stand-alone installation instructions, pretty straightforward

IPC_TXT_OUT=$4   # .out with pi/mw: e.g. propk_acc_ipc.txt
IPC_PNG_OUT=$5   # .png with 2D EF: e.g. propk_acc_ipc.png

echo "Calculating mW/pI and drawing plot..."

# ipc $SEQS ALL $IPC_TXT_OUT $IPC_PNG_OUT
python2 ipc-1.0/ipc/ipc.py $SEQS ALL $IPC_TXT_OUT $IPC_PNG_OUT

rm $IPC_TXT_OUT  # Remove file with pI/mW data

# echo "Done"

echo "CONFIRMATION:"

echo "Number of lines in patterns .out file:"
wc -l $HEADERLIST

echo "Number of sequences in final .faa file:"
grep "^>" $SEQS | wc -l

# Flags:
# l_ldh_seqs.faa: protein sequences in fasta
# ALL: pKa models sets used to report pI
# l_ldh_ipc.txt: output with fasta header, pI and mW
# l_ldh.png: custom 2D EF plot

# Close Python2 environment
echo "Back to Python 3..."
conda deactivate

echo "All done"

