#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=48:00:00
#SBATCH -p agsmall
#SBATCH -o rfmix_all.out
#SBATCH -e rfmix_all.err
#SBATCH --job-name rfmix_all

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4

module load R/4.4.0-openblas-rocky8

Rscript ${path_to_repo}/src/gai2.R ${WORK} ${NAME}

echo "Potential duplication below"
mkdir $WORK/PCA
cp $WORK/study.$NAME.unrelated.comm.popu $WORK/PCA/study.$NAME.unrelated.comm.popu
cd $WORK/PCA

# awk -F '\t' '{print $3}' $WORK/study.$NAME.unrelated.comm.popu | sort | uniq -c > subpop.txt
awk '{print $3}' study.$NAME.unrelated.comm.popu | sort | uniq -c > subpop.txt
awk '{print $1 "\t" $2 "\t" $3}' study.$NAME.unrelated.comm.popu > data.txt
Rscript ${path_to_repo}/src/subpop.R ${WORK} ${NAME}

mkdir $WORK/PCA
cp $WORK/study.$NAME.unrelated.comm.popu $WORK/PCA/study.$NAME.unrelated.comm.popu
cd $WORK/PCA

# awk -F '\t' '{print $3}' $WORK/study.$NAME.unrelated.comm.popu | sort | uniq -c > subpop.txt
awk '{print $3}' study.$NAME.unrelated.comm.popu | sort | uniq -c > subpop.txt
awk '{print $1 "\t" $2 "\t" $3}' study.$NAME.unrelated.comm.popu > data.txt
Rscript ${path_to_repo}/src/subpop.R ${WORK} ${NAME}
#rm *.dat
