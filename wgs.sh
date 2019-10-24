# Script métagenomique TP 2

# Sample G

mkdir results/
mkdir results/sam
mkdir results/bam
mkdir results/camembert
mkdir results/genes
mkdir results/blast

soft/bowtie2-build databases/all_genome.fasta databases/all_genome.fasta

soft/bowtie2 -p 6 --fast --end-to-end -x databases/all_genome.fasta --1 fastq/EchG_R1.fastq.gz --2 fastq/EchG_R2.fastq.gz -S results/sam/EchG.sam

samtools view --threads 6 -b -1 results/sam/EchG.sam -o results/bam/EchG.bam 

samtools sort --threads 6 results/bam/EchG.bam -o results/bam/EchG.sorted.bam

samtools index results/bam/EchG.sorted.bam

samtools idxstats results/bam/EchG.sorted.bam > results/bam/EchG.idxstats.txt

grep ">" databases/all_genome.fasta|cut -f 2 -d ">" >results/camembert/association.txt

sed 's/ /\t/' results/camembert/association.txt > results/camembert/association.tsv 

soft/megahit --mem-flag 0 --k-list 21 -1 fastq/EchG_R1.fastq.gz -2 fastq/EchG_R2.fastq.gz -o results/assemblage

prodigal -i results/assemblage/final.contigs.fa -d results/genes/predicted_genes.fasta

sed "s:>:*\n>:g" results/genes/predicted_genes.fasta | sed -n "/partial=00/,/*/p"|grep -v "*" > results/genes/genes_full.fna

blastn -query results/genes/genes_full.fna -db databases/resfinder.fna -perc_identity 80 -qcov_hsp_perc 80 -evalue 0.001 -out results/blast/EchG_blastn.txt -outfmt '6 qseqid sseqid pident qcovs evalue' -best_hit_score_edge 0.001
