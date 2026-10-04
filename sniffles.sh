module load sniffles
sniffles --threads 4 --tmp-dir /lscratch/$SLURM_JOB_ID --allow-overwrite --phase \
  --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
  --input ABCA4.Arno_UCL.phased.bam --vcf ABCA4.Arno_UCL.sniffles.vcf
  
#Gavin's bam file used chr name as NC_000001.11 etc., what's the reference used?
#genome_assemblies_genome_fasta/ncbi-genomes-2020-09-07/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
#minimap2 -ax map-ont /Users/gavin/genome_assemblies_genome_fasta/ncbi-genomes-2020-09-07/GCF_000001405.39_GRCh38.p13_genomic.fna.gz ABCA4_KH.fastq
