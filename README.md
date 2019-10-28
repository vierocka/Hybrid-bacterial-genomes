The script is written in bash and uses arrays to adjust a reference fasta file for called SNPs (point mutations). The algorithm has been used to improve the mappability of RNA-seq reads for newly emerged bacterial hybrids.

**Problem**  
Two closely related bacterial strains (subspecies, species) recombined and created a new bacterial hybrid. The two bacterial strains are called a donor and a recipient. We have DNA-seq reads of the donor and the hybrid, as well as RNA-seq reads of the hybrid. We want to reconstruct the reference genome of the hybrid to increase the precision of RNA-seq reads mapping.

**Aim**  
to create a new reference of genome (chromosome) of the hybrid bacterial strain (that helps to increase mapping precision of RNA-seq reads)

**The prerequisites**  
1. the reference genome of the donor (a fasta file from a public DB)
2. a file with identified SNPs (a list of identified nucleotide differences between the reference and the hybrid)
 
**Pipeline to get the second prerequisity (the vcf file)**  
- bwa mem (-O 1 -B 1): a loose mapping of DNA-seq reads from the hybrid against the reference genome (donor)
- GATK: SNP calling and filtering
- custom-made script: filtering 

**Solution**  
one the bash scripts: adjustments of the reference genome (donor's fasta file)

The scripts created to prepare adjusted genome references of hybrid bacterial strains:  
- reference adjustment (only point mutation)  
See the script "Hybrid genome assembly"
- advance reference adjustment (point mutations and indels) 
- advance reference annotation (point mutations and indels)
