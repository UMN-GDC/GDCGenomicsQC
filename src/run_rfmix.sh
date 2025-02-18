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

mkdir $WORK/phased
sbatch --time 4:00:00 --mem 16GB --array 1-22 --wait -N1 ${path_to_repo}/src/phase_individual.sh ${WORK} ${NAME}

mkdir $WORK/rfmix
sbatch --time 18:00:00 --mem 64GB --array 1-22 --wait -N1 ${path_to_repo}/src/rfmix_individual.sh ${WORK} ${NAME}

module load R/4.4.0-openblas-rocky8

Rscript ${path_to_repo}/src/gai.R ${WORK} ${NAME}

mkdir $WORK/PCA
cp $WORK/study.$NAME.unrelated.comm.popu $WORK/PCA/study.$NAME.unrelated.comm.popu
cd $WORK/PCA

# Added in below to make the 1000G.comm file
# Reference panel has no rsid nums...
# $REF/Fraposa/commvar.sh $REF/hgdp1kg.filtered.SNV_INDEL.38.phased.shapeit5/ALL.hgdp1kg.filtered.SNV_INDEL.38.phased.shapeit5 $WORK/relatedness/study.$NAME.unrelated hgdp1kg.comm study.$NAME.unrelated.comm 
# $REF/Fraposa/predstupopu.py hgdp1kg.comm study.$NAME.unrelated.comm
# $REF/Fraposa/plotpcs.py hgdp1kg.comm study.$NAME.unrelated.comm

$REF/Fraposa/commvar.sh $REF/PCA_ref/1000G.aligned $WORK/relatedness/study.$NAME.unrelated 1000G.comm study.$NAME.unrelated.comm
$REF/Fraposa/fraposa_runner.py --stu_filepref study.$NAME.unrelated.comm 1000G.comm #Main program for Fraposa 
$REF/Fraposa/predstupopu.py 1000G.comm study.$NAME.unrelated.comm # To get this to work will need to run all of the standard fraposa steps... 
$REF/Fraposa/plotpcs.py 1000G.comm study.$NAME.unrelated.comm
mv $WORK/PCA/study.$NAME.unrelated.comm.popu $WORK/PCA/fraposa_study.$NAME.unrelated.comm.popu
cp $WORK/study.$NAME.unrelated.comm.popu $WORK/PCA/study.$NAME.unrelated.comm.popu

# awk -F '\t' '{print $3}' $WORK/study.$NAME.unrelated.comm.popu | sort | uniq -c > subpop.txt
awk '{print $3}' study.$NAME.unrelated.comm.popu | sort | uniq -c > subpop.txt
awk '{print $1 "\t" $2 "\t" $3}' study.$NAME.unrelated.comm.popu > data.txt
Rscript ${path_to_repo}/src/subpop.R ${WORK} ${NAME}
#rm *.dat
