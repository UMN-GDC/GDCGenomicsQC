#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=24:00:00
#SBATCH -p msismall
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tran1327@umn.edu
#SBATCH -o genome_harmonizer.%j.out
#SBATCH -e genome_harmonizer.%j.err
#SBATCH --job-name genome_harmonizer
#SBATCH --array 0-23

module load java/openjdk-17.0.2
module load plink

REF=/panfs/jay/groups/9/gdc/shared/GDC_pipeline/Ref
WORK=/panfs/jay/groups/9/gdc/shared/GDC_pipeline/QC_example
NAME="SMILES_GDA"
OUTDIR="aligned"

mkdir -p $WORK/$OUTDIR

DATA=($(find "$WORK/lifted" -name "study.$NAME.lifted.chr*.bed" | sort -V))
LIST=${DATA[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename $LIST | sed "s/.bed//g")

java -Xmx1g -jar $REF/GenotypeHarmonizer/GenotypeHarmonizer.jar --input $WORK/lifted/$SAMPLE --inputType PLINK_BED --ref $REF/1000G_GRCh38/ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased --refType VCF --keep --output $WORK/$OUTDIR/${SAMPLE}.aligned
