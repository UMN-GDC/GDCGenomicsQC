#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=48:00:00
#SBATCH -p agsmall
#SBATCH -o rfmix_plot.out
#SBATCH -e rfmix_plot.err
#SBATCH --job-name rfmix_plot

WORK=$1
REF=$2
NAME=$3
path_to_repo=$4

mkdir $WORK/visualization
cd $WORK/visualization

n_subs=$(wc -l < $WORK/relatedness/study.$NAME.unrelated.fam)
n_rfmix_rows=$(wc -l < $WORK/rfmix/ancestry_chr1.rfmix.Q)

for chr in {1..22}; do
    input_file="$WORK/rfmix/ancestry_chr${chr}.rfmix.Q"

    for ind in $(seq 3 "$n_rfmix_rows"); do
        individual_index=$((ind - 2))
        output_file="$WORK/visualization/ancestry${individual_index}_chr${chr}.rfmix.Q"
        sed -n -e "1p" -e "2p" -e "${ind}p" "$input_file" > "$output_file"
    done
done

python /home/gdc/shared/RFMIX2-Pipeline-to-plot/GAP/Scripts/RFMIX2ToBed4GAP.py --prefix $WORK/visualization/ancestry --output $WORK/visualization
python /home/gdc/shared/RFMIX2-Pipeline-to-plot/GAP/Scripts/BedToGap.py --input ancestry.bed --out ancestry_GAP.bed
# python /home/gdc/shared/RFMIX2-Pipeline-to-plot/GAP/Scripts/GAP.py --input ancestry_GAP.bed --output ancestry_GAP.pdf

input_bed="ancestry_GAP.bed"
gap_script="/home/gdc/shared/RFMIX2-Pipeline-to-plot/GAP/Scripts/GAP.py"
output_dir="./GAP_individual_plots"

mkdir -p "$output_dir"

# Extract header
head -n 2 "$input_bed" > header.tmp

# Get unique individual IDs (assumes 1st column is the ID)
tail -n +3 "$input_bed" | cut -f1 | sort -u > individual_ids.txt

batch_num=1
ids_batch=()

# Process IDs in batches of 10
while read -r id; do
    ids_batch+=("$id")

    if [ "${#ids_batch[@]}" -eq 10 ]; then
        batch_file="$output_dir/batch_${batch_num}.bed"
        batch_pdf="$output_dir/batch_${batch_num}.pdf"

        # Extract data for all individuals in this batch
        cat header.tmp > "$batch_file"
        for id_in_batch in "${ids_batch[@]}"; do
            grep -w "$id_in_batch" "$input_bed" >> "$batch_file"
        done

        # Run GAP
        echo "Generating $batch_pdf..."
        python "$gap_script" --input "$batch_file" --output "$batch_pdf"

        # Cleanup
        rm "$batch_file"
        ids_batch=()
        ((batch_num++))
    fi
done < individual_ids.txt

# Handle final batch (if fewer than 10 left)
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

# Final cleanup
rm header.tmp individual_ids.txt

mv ./GAP_individual_plots $WORK/GAP_plots

mv ancestry.bed ancestry_GAP_posterior.bed
cp ancestry_GAP_posterior.bed $WORK/GAP_plots/ancestry_GAP_posterior.bed


for chr in {1..22}; do
    input_file="/scratch.global/and02709/gdc/rfmix/ancestry_chr${chr}.msp.tsv"

    for ind in $(seq 1 "$n_subs"); do
        ind1=$(( (2 * ind - 1) + 6 ))
        ind2=$(( (2 * ind) + 6 ))
        output_file="/scratch.global/and02709/gdc/visualization/ancestry${ind}_chr${chr}.msp.tsv"
        
        # Debugging
        echo "Processing: $output_file (Columns: 1-6, $ind1, $ind2)"
        
        cut -f1-6,$ind1,$ind2 "$input_file" > "$output_file"
    done
done


python /home/gdc/shared/RFMIX2-Pipeline-to-plot/LAP/Scripts/RFMIX2ToBed.py --prefix $WORK/visualization/ancestry --output $WORK/visualization

mkdir $WORK/LAP_plots

for ind in $(seq 1 "$n_subs"); do
    input_file_1="/scratch.global/and02709/gdc/visualization/ancestry${ind}_hap1.bed"
    input_file_2="/scratch.global/and02709/gdc/visualization/ancestry${ind}_hap2.bed"
    output_file="/scratch.global/and02709/gdc/LAP_plots/ancestry${ind}.bed"
    output_LAP="/scratch.global/and02709/gdc/LAP_plots/ancestry${ind}.pdf"
        
    # Debugging
    echo "Processing: $output_file"
        
    python /home/gdc/shared/RFMIX2-Pipeline-to-plot/LAP/Scripts/BedToLAP.py --bed1 "$input_file_1" --bed2 "$input_file_2" --out "$output_file"
    python /home/gdc/shared/RFMIX2-Pipeline-to-plot/LAP/Scripts/LAP.py -I "$output_file" -O "$output_file"
done
# python /home/gdc/shared/RFMIX2-Pipeline-to-plot/LAP/Scripts/BedToLAP.py --input ancestry.bed --out ancestry_GAP.bed
# python /home/gdc/shared/RFMIX2-Pipeline-to-plot/LAP/Scripts/GAP.py --input ancestry_GAP.bed --output ancestry_GAP.pdf


