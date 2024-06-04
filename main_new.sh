!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH -p msismall
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=x500@umn.edu 
#SBATCH -o x.out
#SBATCH -e x.err
#SBATCH --job-name x

source /home/faird/shared/code/external/envs/miniconda3/load_miniconda3.sh
conda activate GDC_pipeline
module load plink
module load perl
module load R

# This pipeline assumes the input is in plink binary format.

#################################### Specifying paths #########################################

# Hard-code the path to the Reference folder (containing reference dataset, other bash scripts, and programs' executables like CrossMap, GenomeHarmonizer, PRIMUS, and fraposa)
REF=/home/gdc/shared/GDC_pipeline/Ref

# Hard-code an example dataset (SMILES). This should be obtained from user's input in the future
FILE=/home/gdc/shared/GDC_pipeline/Study
NAME="SMILES_GDA"

# Hard-code the working directory (where outputs should be located in)
WORK=/home/gdc/shared/GDC_pipeline/test
cd $WORK
#################################################################################################



############## Updating genome build and conducting strand alignment/allele flipping #############

echo "(Step 1) Matching data to NIH's GRCh38 genome build"
# LiftOver/CrossMap: results are so far identical for both programs. 
# Using LiftOver tool for now as it requires less steps given plink format.
# Since plink denote X chromosome's pseudo-autosomal region as a separate 'XY' chromosome, we want to merge to pass ontto LiftOver/CrossMap. 
# We also reformat the numeric chromsome {1-26} to {1-22, X, Y, MT} for LiftOver/CrossMap
plink --bfile $FILE/$NAME --merge-x --make-bed --out prep1
plink --bfile prep1 --recode --output-chr 'MT' --out prep2

rm prep.bed updated.snp updated.position updated.chr
awk '{print $1, $4-1, $4, $2}' prep2.map > prep.bed
## Stuck at this spot ... 
python $REF/CrossMap/CrossMap.py bed $REF/CrossMap/GRCh37_to_GRCh38.chain.gz prep.bed study.$NAME.lifted.bed3

awk '{print $4}' study.$NAME.lifted.bed3 > updated.snp
awk '{print $4, $3}' study.$NAME.lifted.bed3 > updated.position
awk '{print $4, $1}' study.$NAME.lifted.bed3 > updated.chr
plink --file prep2 --extract updated.snp --make-bed --out result1
plink --bfile result1 --update-map updated.position --make-bed --out result2
plink --bfile result2 --update-chr updated.chr --make-bed --out result3
plink --bfile result3 --recode --out study.$NAME.lifted
# Break the dataset by chromosomes for faster processing in the next step (genome harmonizer)
mkdir $WORK/lifted
for chr in {1..22} X Y; do plink --file study.$NAME.lifted --chr $chr --make-bed --out $WORK/lifted/study.$NAME.lifted.chr${chr};  done
rm prep1.* prep2.* result1.* result2.* result3.* prep.bed updated.snp updated.position updated.chr

# Using genome harmonizer, update strand orientation and flip alleles according to the reference dataset.
srun $REF/harmonizer.job
# Currently reference dataset does not have chrY for alignment, and ChrX has no match with study data
# Hence, we bring the unaligned ChrX and ChrY to the result folder, i.e. skipping alignment
cp $WORK/lifted/study.SMILES_GDA.lifted.chrX.bed $WORK/aligned/study.SMILES_GDA.lifted.chrX.aligned.bed
cp $WORK/lifted/study.SMILES_GDA.lifted.chrX.bim $WORK/aligned/study.SMILES_GDA.lifted.chrX.aligned.bim
cp $WORK/lifted/study.SMILES_GDA.lifted.chrX.fam $WORK/aligned/study.SMILES_GDA.lifted.chrX.aligned.fam
cp $WORK/lifted/study.SMILES_GDA.lifted.chrY.bed $WORK/aligned/study.SMILES_GDA.lifted.chrY.aligned.bed
cp $WORK/lifted/study.SMILES_GDA.lifted.chrY.bim $WORK/aligned/study.SMILES_GDA.lifted.chrY.aligned.bim
cp $WORK/lifted/study.SMILES_GDA.lifted.chrY.fam $WORK/aligned/study.SMILES_GDA.lifted.chrY.aligned.fam
#######################################################################################################


###################################### QC #############################################################

echo "(Step 2) Standard variants and samples filtering"
# Merge chromosomes for this step
cd $WORK/aligned
rm mergelist.txt
for chr in {2..22} X Y; do echo study.$NAME.lifted.chr${chr}.aligned >> mergelist.txt; done
plink --bfile study.$NAME.lifted.chr1.aligned --merge-list mergelist.txt --allow-no-sex --make-bed --out study.$NAME.lifted.aligned1
plink --bfile study.$NAME.lifted.aligned1 --split-x 'hg38' --make-bed --out study.$NAME.lifted.aligned
# Run standard_QC.job with the appropriate parameters (full path to dataset name + output folder name)
cd $WORK
DATATYPE=full
srun $REF/standard_QC.job $WORK/aligned/study.$NAME.lifted.aligned $DATATYPE
########################################################################################################


######################################## Pedigree ######################################################

echo "(Step 3) Relatedness check"
mkdir $WORK/relatedness
perl $REF/PRIMUS/bin/run_PRIMUS.pl --file $WORK/${DATATYPE}/$DATATYPE.QC8 --genome -t 0.2 -o $WORK/relatedness
OUT=$WORK/relatedness/$DATATYPE.QC8_cleaned.genome_unrelated_samples.txt
# Reformat the unrelated set text file in a suitable format for plink --keep
tail -n +2 "$OUT" > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
awk '{print "0", $1}' $OUT > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
# Keep only the unrelated set of individuals determined by PRIMUS
plink --bfile $WORK/$DATATYPE/$DATATYPE.QC8 --keep $WORK/relatedness/$DATATYPE.QC8_cleaned.genome_unrelated_samples.txt --make-bed --out $WORK/relatedness/study.$NAME.unrelated
#########################################################################################################


######################################## Ethnicity ######################################################
echo "(Step 4) PCA"
mkdir $WORK/PCA
cd $WORK/PCA

# Load the required packages for fraposa
pip install numpy pandas scikit-learn pyplink matplotlib rpy2
Rscript $REF/fraposaRpackage.R

# fraposa operations
$REF/Fraposa/commvar.sh $REF/PCA_ref/1000G.aligned $WORK/relatedness/study.$NAME.unrelated 1000G.comm study.$NAME.unrelated.comm
$REF/Fraposa/fraposa_runner.py --stu_filepref study.$NAME.unrelated.comm 1000G.comm
$REF/Fraposa/predstupopu.py 1000G.comm study.$NAME.unrelated.comm
$REF/Fraposa/plotpcs.py 1000G.comm study.$NAME.unrelated.comm

awk -F '\t' '{print $3}' study.$NAME.unrelated.comm.popu | sort | uniq -c > subpop.txt
awk '{print $1 "\t" $2 "\t" $3}' study.$NAME.unrelated.comm.popu > data.txt
Rscript $REF/subpop.R
rm *.dat
#########################################################################################################


################### Subset data based on Ethnicity and Rerun QC (Step 2) on the subsets #################
DATATYPE="EUR"
cd $WORK
plink --bfile $WORK/aligned/study.$NAME.lifted.aligned --keep $WORK/PCA/$DATATYPE --make-bed --out $WORK/aligned/study.$NAME.$DATATYPE.lifted.aligned
srun $REF/standard_QC.job $WORK/aligned/study.$NAME.$DATATYPE.lifted.aligned $DATATYPE

DATATYPE="AFR"
cd $WORK
plink --bfile $WORK/aligned/study.$NAME.lifted.aligned --keep $WORK/PCA/$DATATYPE --make-bed --out $WORK/aligned/study.$NAME.$DATATYPE.lifted.aligned
srun $REF/standard_QC.job $WORK/aligned/study.$NAME.$DATATYPE.lifted.aligned $DATATYPE

DATATYPE="AMR"
cd $WORK
plink --bfile $WORK/aligned/study.$NAME.lifted.aligned --keep $WORK/PCA/$DATATYPE --make-bed --out $WORK/aligned/study.$NAME.$DATATYPE.lifted.aligned
srun $REF/standard_QC.job $WORK/aligned/study.$NAME.$DATATYPE.lifted.aligned $DATATYPE

###########################################################################################################

#file_to_read=/home/miran045/shared/projects/SCAN_PD/experiments/make_dconns/subj_ids_trio.txt
#
#num_lines=$(cat ${file_to_read} | wc -l)
#
#for (( i=1; i<=$num_lines; i++)); do
#    #Subject
#    line_content=$(sed "${i}q;d" "$file_to_read")
#    participant_name=$(echo "$line_content" | tr -d '[:space:]')  # Remove spaces if needed
#    participant_name_1=sub-${participant_name}
#    
#    #Session
#    line2_content=$(sed "${i}q;d" "$sessions_to_follow")
#    session_using=$(echo "$line2_content" | tr -d '[:space:]')
#    session_using_1=ses-${session_using}
#done
