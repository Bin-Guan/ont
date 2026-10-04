#!/bin/bash
#SBATCH -c32
#SBATCH --mem=32g
#SBATCH --partition=gpu
#SBATCH --gres=lscratch:100,gpu:1 --constraint='gpuv100|gpua100|gpuv100x'
#SBATCH --time=6:0:0


set -e
module load dorado/0.9.6

ref=$1
pod5Dir=$2
outputBam=$3

#followed protocol here https://www.protocols.io/view/bioinformatic-pipeline-for-analysing-variations-in-36wgqnoz3gk5/v1?step=1

dorado basecaller --reference $ref \
${DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
$pod5Dir --min-qscore 10 --recursive --modified-bases-models \
${DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v3 \
> $outputBam

dorado summary $outputBam > sequencing_summary_dorado.txt

#sbatch ~/git/ont/genome/dorado_pod5.sh resource/ref/mm39.TCR.mmi secondRun/R161H_RH6878/pod5/ analysis/R161H.dorado.bam