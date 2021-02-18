#!/bin/bash

# Usage example:

# Single sample:
# bash ~/HGM_MG_sample_analyzer.sh ~/loomba17_study/samples/SRR5275394_1.fastq.gz ~/GIF_DBs/acetolac_seqs_lenfil_fil.dmnd

# Multiple samples:
# for i in ~/ANALYSIS/IBD_cohort/final_samples/*.fastq.gz; do bash ~/HGM_MG_sample_analyzer.sh $i ~/GIF_DBs/ldh_seqs.dmnd; done

# Input defined variables:
FQ_SAMPLE=$1  # .fq left reads, .gz compressed: e.g. "SP-099.human.filtered_1.fq.gz" or "SRR5275394_1.fastq.gz"
DB=$2         # .dmnd file with refined-protein set DB: e.g. acetolac_seqs_lenfil_fil.dmnd
# Variables:
FA_SAMPLE="$(basename "$FQ_SAMPLE" .fastq.gz).faa"        # .fa left reads
DIAMOND_OUT="$(basename "$FQ_SAMPLE" _R1.fastq.gz)_matches.m8" # DIAMOND matches

##################################################################

echo "Analyzing sample $FQ_SAMPLE..."

#echo "$FQ_SAMPLE" >> lengths.txt
#echo "$FQ_SAMPLE" >> kkreads.txt

# 1. Illumina reads features:
#echo "1. Illumina reads features:"

# 1.A. Median read length (20 first reads):
#echo "1.A. Median read length (20 first reads):"
#gunzip -c $FQ_SAMPLE | awk 'NR%4==2{print length($0)}' | head -n20 | sort -n | awk '{a[i++]=$1;} END {print a[int(i/2)];}' >> lengths.txt

# 1.B. # reads:
#echo "1.B. # kk reads:"
#gunzip -c $FQ_SAMPLE | awk 'END {print NR/4000000}' >> kkreads.txt

####################################################################

# 2. Convert FQ to FA (SEQKIT)
echo "2. Convert sample $FQ_SAMPLE to .fasta format:"
seqkit fq2fa $FQ_SAMPLE -o $FA_SAMPLE

# 3. Translated alignment of reads to GIF-enzymatic reference database (DIAMOND)
echo "3. Translated alignment of sample $FQ_SAMPLE to GIF-enzymatic reference database:"
~/diamond/diamond blastx --evalue 0.001 --max-hsps 1 -k0 --id 80 \
--outfmt 6 qseqid sseqid evalue stitle pident length mismatch gapopen \
--db $DB --query $FA_SAMPLE --out $DIAMOND_OUT

rm $FA_SAMPLE

# --max-hsps 1: report a single HSP for a query/target pair
# But multiple alignments for a single query read with different targets are allowed
# -k0: report all alignments that were found
# --id 80: report only alignments above 80% of sequence identity
# --evalue 0.001: report only alignments below 0.001 expected value (default=0.001)

# Columns format details:
# 1. qseqid: read ID1
# 2. sseqid: target ID
# 3. evalue: expected value
# 4. stitle: target title
# 5. pident: % of identical matches
# 6. length: alignment length
# 7. mismatch: # missmatches
# 8. gapopen: # gap openings

echo ""
echo "Sample $FQ_SAMPLE analyzed!"
