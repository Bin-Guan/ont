#!/bin/bash
#SBATCH -c8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:50
#SBATCH --time=2:0:0

#Could submit with dependency then start automatically after the previous step completes. 
#sbatch --dependency=afterok:11254323 ~/git/NGS_genotype_calling/dev/bcm_RPGR.porechop.minimap2.alignment.sh OGL
#~20 min per sample, thus adjust the time accordingly. 
#fastq files have to be SampleID-LW and SampleID-MW. If missing either LW or MW, touch an empty file.

set -e
#library=$1 #QBplEx or OGL; will be part of @RG
metadata=$1 #excel file with sample name, target name, Fwd primer and reverse primer.
batchname=$(date '+%Y%m%d_%H%M') 
WORK_DIR=/lscratch/$SLURM_JOB_ID

##Have not started editing
module load R/4.3.0
Rscript ~/git/NGS_genotype_calling/ont/phasing_metadata.R $metadata $batchname.cDNA.metadata.tsv

module load porechop/0.2.4 minimap2/2.26 
module load samtools/1.21

mkdir -p bam porechop

while read -r line;
do
SAMPLE=$(echo $line | cut -d" " -f 1)
TARGET=$(echo $line | cut -d" " -f 2 | sed 's/_//' )
FWD_PRIMER=$(echo $line | cut -d" " -f 3 | tr 'atcg' 'ATCG')
RevCom=$(echo $FWD_PRIMER | rev | tr 'ATCG' 'TAGC')
export PORECHOP_ADAPTERS=$SAMPLE.$TARGET.porechop.adapters.py
sed -e "s/TARGETNAME/$TARGET/" -e "s/PRIMER_SEQ/$FWD_PRIMER/" -e "s/PRIMER_REVCOMPL/$RevCom/" ~/git/ont/amplicon/porechop/porechop.adapters.py > $SAMPLE.$TARGET.porechop.adapters.py
porechop --discard_unassigned --discard_middle --extra_end_trim 0 -i fastq/$SAMPLE.fastq.gz -b porechop/fwd/$SAMPLE
REV_PRIMER=$(echo $line | cut -d" " -f 4 | tr 'atcg' 'ATCG')
RevCom=$(echo $REV_PRIMER | rev | tr 'ATCG' 'TAGC')
sed -e "s/TARGETNAME/$TARGET/" -e "s/PRIMER_SEQ/$REV_PRIMER/" -e "s/PRIMER_REVCOMPL/$RevCom/" ~/git/ont/amplicon/porechop/porechop.adapters.py > $SAMPLE.$TARGET.porechop.adapters.py
porechop --discard_unassigned --discard_middle --extra_end_trim 0 -i porechop/fwd/$SAMPLE/BC$TARGET.fastq.gz -b porechop/trimmed/$SAMPLE
rm $SAMPLE.$TARGET.porechop.adapters.py
minimap2 -ax splice:hq -t 3 -R "@RG\tLB:OGL\tID:$SAMPLE\tSM:$SAMPLE\tPL:ONT" /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi porechop/trimmed/$SAMPLE/BC$TARGET.fastq.gz \
| samtools sort -@ 3 -o bam/$SAMPLE.bam 
samtools index -@ 3 bam/$SAMPLE.bam
done < $batchname.cDNA.metadata.tsv

rm -r porechop