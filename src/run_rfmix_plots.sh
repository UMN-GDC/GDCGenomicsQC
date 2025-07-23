#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=1:00:00
#SBATCH -p agsmall
#SBATCH -o rfmix_plot.out
#SBATCH -e rfmix_plot.err
#SBATCH --job-name rfmix_plot

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4

mkdir -p $WORK/visualization
cd $WORK/visualization

n_rfmix_rows=$(wc -l < "$WORK/rfmix/ancestry_chr1.rfmix.Q")

# Extract individual IDs from ancestry_chr1.rfmix.Q, skipping first 2 lines
mapfile -t sample_ids < <(tail -n +3 "$WORK/rfmix/ancestry_chr1.rfmix.Q" | cut -f1)

# Write header for the mapping file
mapping_file="$WORK/visualization/ancestry_index_map.tsv"
echo -e "individual_index\tindividual_id" > "$mapping_file"

# Process ancestry_chr*.rfmix.Q files
for chr in {1..22}; do
    input_file="$WORK/rfmix/ancestry_chr${chr}.rfmix.Q"

    for ind in $(seq 3 "$n_rfmix_rows"); do
        individual_index=$((ind - 2))
        output_file="$WORK/visualization/ancestry${individual_index}_chr${chr}.rfmix.Q"
        sed -n -e "1p" -e "2p" -e "${ind}p" "$input_file" > "$output_file"

        if [ "$chr" -eq 1 ]; then
            sample_index=$((individual_index - 1))
            echo -e "${individual_index}\t${sample_ids[$sample_index]}" >> "$mapping_file"
        fi
    done
done

# GAP Pipeline
python ${REF}/RFMIX2-Pipeline-to-plot/GAP/Scripts/RFMIX2ToBed4GAP.py --prefix $WORK/visualization/ancestry --output $WORK/visualization
python ${REF}/RFMIX2-Pipeline-to-plot/GAP/Scripts/BedToGap.py --input ancestry.bed --out ancestry_GAP.bed

input_bed="ancestry_GAP.bed"
gap_script="${REF}/RFMIX2-Pipeline-to-plot/GAP/Scripts/GAP.py"
output_dir="./GAP_individual_plots"

mkdir -p "$output_dir"
head -n 2 "$input_bed" > header.tmp
tail -n +3 "$input_bed" | cut -f1 | sort -u > individual_ids.txt

batch_num=1
ids_batch=()

while read -r id; do
    ids_batch+=("$id")
    if [ "${#ids_batch[@]}" -eq 10 ]; then
        batch_file="$output_dir/batch_${batch_num}.bed"
        batch_pdf="$output_dir/batch_${batch_num}.pdf"
        cat header.tmp > "$batch_file"
        for id_in_batch in "${ids_batch[@]}"; do
            grep -w "$id_in_batch" "$input_bed" >> "$batch_file"
        done
        echo "Generating $batch_pdf..."
        python "$gap_script" --input "$batch_file" --output "$batch_pdf"
        rm "$batch_file"
        ids_batch=()
        ((batch_num++))
    fi
done < individual_ids.txt

if [ "${#ids_batch[@]}" -gt 0 ]; then
    batch_file="$output_dir/batch_${batch_num}.bed"
    batch_pdf="$output_dir/batch_${batch_num}.pdf"
    cat header.tmp > "$batch_file"
    for id_in_batch in "${ids_batch[@]}"; do
        grep -w "$id_in_batch" "$input_bed" >> "$batch_file"
    done
    echo "Generating $batch_pdf..."
    python "$gap_script" --input "$batch_file" --output "$batch_pdf"
    rm "$batch_file"
fi

mv ./GAP_individual_plots $WORK/GAP_plots

# Replace numeric IDs in ancestry.bed using ancestry_index_map.tsv
module load R/4.4.0-openblas-rocky8
Rscript ${path_to_repo}/src/combine_ids_script.R $WORK
cp ancestry.txt $WORK/ancestry_${NAME}.txt

# Copy final result
cp ancestry_labeled.bed $WORK/GAP_plots/ancestry_GAP_posterior_labeled.bed
cp ancestry.bed $WORK/GAP_plots/ancestry_GAP_posterior.bed


# Process MSPs for LAP
n_subs=${#sample_ids[@]}

for chr in {1..22}; do
    input_file="$WORK/rfmix/ancestry_chr${chr}.msp.tsv"
    for ind in $(seq 1 "$n_subs"); do
        ind1=$(( (2 * ind - 1) + 6 ))
        ind2=$(( (2 * ind) + 6 ))
        output_file="$WORK/visualization/ancestry${ind}_chr${chr}.msp.tsv"
        echo "Processing: $output_file (Columns: 1-6, $ind1, $ind2)"
        cut -f1-6,$ind1,$ind2 "$input_file" > "$output_file"
    done
done

# LAP Pipeline
python ${REF}/RFMIX2-Pipeline-to-plot/LAP/Scripts/RFMIX2ToBed.py --prefix $WORK/visualization/ancestry --output $WORK/visualization
mkdir -p $WORK/LAP_plots

for ind in $(seq 1 "$n_subs"); do
    input_file_1="$WORK/visualization/ancestry${ind}_hap1.bed"
    input_file_2="$WORK/visualization/ancestry${ind}_hap2.bed"
    output_file="$WORK/LAP_plots/ancestry${ind}.bed"
    output_LAP="$WORK/LAP_plots/ancestry${ind}.pdf"
    echo "Processing: $output_file"
    python ${REF}/RFMIX2-Pipeline-to-plot/LAP/Scripts/BedToLAP.py --bed1 "$input_file_1" --bed2 "$input_file_2" --out "$output_file"
    python ${REF}/RFMIX2-Pipeline-to-plot/LAP/Scripts/LAP.py -I "$output_file" -O "$output_LAP"
done
