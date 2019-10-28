### HYBRID GENOME ASSEMBLY
# the script merely focuses on point mutations 
# genome adjustment using arrays in bash #
# input: a filtered vcf file and a referebce fasta file
# alternative input: a file with positions and point mutations (SNPs), and a referebce fasta file

mkdir hybrid
# call a list of point mutations and their positions
cut -f2,5 Your_filtered.vcf | grep -P "\t[ACGT]{1}$" > Your_filtered_poin_mutations.csv
 
# create an array with one bacterial chromosome
declare Seq
eval Seq=( $( cat referece_genome.fa | sed '/^>/d' | sed '/^$/d' | tr -d '\n' | awk 'BEGIN{FS=""}; { for (i=1; i<length+1; i++) { print $i}} ' ))
# create another array with one bacterial chromosome
declare  SeqOrig
eval SeqOrig=( $( cat referece_genome.fa | sed '/^>/d' | sed '/^$/d' | tr -d '\n' | awk 'BEGIN{FS=""}; { for (i=1; i<length+1; i++) { print $i}} ' ))

# get the total number of SNPs
countSNP=$( wc -l Your_filtered.vcf | cut -d" " -f1) 
# for loop from the first to the last SNP
for ((j=1; j<$(($countSNP+1)); j++))
do
MLpos=$(sed -n ''$j'p' Your_filtered_poin_mutations.csv | cut -f1 )
# array starts from 0, positions in VCF files start from 1; check if this assumption is correct
ReindxPos=$(($MLpos-1))
Nt=$(sed -n ''$j'p' Your_filtered_poin_mutations.csv | cut -f2 )
Seq[$ReindxPos]=$Nt
# check if the script runs correctly
# echo -e $ReindxPos"\t"$Nt"\t"${Seq[$ReindxPos]}"\t"${SeqOrig[$ReindxPos]} >> check.log
# it prints a reindexed position, a nucleotide at the reindexed position, a nucleotited reaplaced in the array called Seq (hybrid), an original nucleotide (from the array SeqOrig)
done

# print a fasta header with a chromosome name
echo ">your_chromosome" > your_hybrid.fa
# print the whole sequence of the hybrid chromosome
echo ${Seq[@]} | tr -d ' ' >> your_hybrid.fa

# recall all the declared arrays
unset Seq
unset SeqOrig
# remove the help file with info about SNPs and positions
rm Your_filtered_poin_mutations.csv

exit 0
