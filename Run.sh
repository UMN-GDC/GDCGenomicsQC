#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH -p msismall
#SBATCH --mail-type=FAIL  
#SBATCH --mail-user=and02709@umn.edu 
#SBATCH -o GDCgenomics.out
#SBATCH -e GDCgenomics.err
#SBATCH --job-name GDCgenomics

show_help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "--set_working_directory		Provide a path to where you'd like the outputs to be stored"
  echo "--input_directory		Provide the path to where the bim/bed/fam data is stored"
  echo "--input_file_name 		Provide the common name that ties the bim/bed/fam files together"
  echo "--path_to_github_repo 		Provide the path to the GDCGenomicsQC pipeline repository"
  echo "--user_x500			Provide your x500 samp213@umn.edu so you may receive email updates regarding sbatch submissions"
  echo "--use_crossmap			Enter '1' for if you would like to update your reference genome build from GRCh37 to GRCh38"
  echo "--use_genome_harmonizer 	Enter '1' if you would like to update strand allignment by using genome harmonizer"
  echo "--use_king			Enter '1' if you would like to use king to estimate relatedness"
  echo "--use_rfmix			Enter '1' if you would like to use rfmix to estimate ancestry"
  echo "--make_report			Enter '1' if you would like an automated report to be generated of the qc steps and what was changed"
  echo "--custom_qc			Enter '1' if you would like to use your own settings for the qc steps such as marker and sample filtering"
  echo "--custom_ancestry		Enter '1' if you would like to use your own ancestry assignment algorithm"
  echo "					When providing this flag you will need to answer all of the questions prompted by the terminal"
  echo "Default settings: 		The pipeline by default if flags are not provided will use crossmap, genome harmonizer, fraposa and will generate the automated reports"
  echo "--help				Display this help message."
}


# Check for command line arguments
if [ $# -eq 0 ]; then
  echo "No arguments provided."
  show_help
  exit 1
fi

NYS=1
set_working_directory=$(pwd)
input_directory=$(pwd)
input_file_name=99
path_to_github_repo=$(pwd)
user_x500=99
use_crossmap=1
use_genome_harmonizer=1
use_king=1
use_rfmix=1
make_report=1
custom_qc=0
custom_ancestry=0
flag=0


# *** Make sure you have a new enough getopt to handle long options (see the man page)
getopt -T &>/dev/null
if [[ $? -ne 4 ]]; then echo "Getopt is too old!" >&2 ; exit 1 ; fi

declare {set_working_directory,input_directory,input_file_name,path_to_github_repo,user_x500,use_crossmap,use_genome_harmonizer,use_rfmix,make_report,custom_qc,custom_ancestry,help}
OPTS=$(getopt -u -o '' -a --longoptions 'set_working_directory:,input_directory:,input_file_name:,path_to_github_repo:,user_x500:,use_crossmap:,use_genome_harmonizer:,use_rfmix:,make_report:,custom_qc:,custom_ancestry:,help' -n "$0" -- "$@")
    # *** Added -o '' ; surrounted the longoptions by ''
if [[ $? -ne 0 ]] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi
    # *** This has to be right after the OPTS= assignment or $? will be overwritten

set -- $OPTS
    # *** As suggested by chepner

while [[ $# -gt 0 ]]; do
  key=$1
  case $key in
	--set_working_directory )
		set_working_directory=$2
		shift 2
		;;
	--input_directory )
		input_directory=$2
		shift 2
		;;
	--input_file_name )
		input_file_name=$2
		shift 2
		;;
	--path_to_github_repo )
		path_to_github_repo=$2
		shift 2
		;;
	--user_x500 )
		user_x500=$2
		shift 2
		;;
	--use_crossmap )
        	use_crossmap=$2
        	shift 2
        	;;
	--use_genome_harmonizer )
        	use_genome_harmonizer=$2
        	shift 2
        	;;
	--use_king )
        	use_king=$2
        	shift 2
        	;;
	--use_rfmix )
        	use_rfmix=$2
        	shift 2
        	;;
	--make_report )
        	make_report=$2
        	shift 2
        	;;
	--custom_qc )
        	custom_qc=$2
        	shift 2
        	;;
	 --custom_ancestry )
        	custom_ancestry=$2
        	shift 2
        	;;
    --help )
			show_help
			shift 2
			exit 1
			;;
    --)
			shift
			break
			;;
    	*)
  esac
done

if (( $input_file_name == 99)); then  
	echo "Please provide plink binary file root name" 
	flag=99
fi

if (( $user_x500 == 99)); then
	echo "Please provide user x500" 
	flag=99
fi

if  (( $flag == 99)); then 
	echo "Please see --help" 
	exit 1 
fi

if [ ! -f "${path_to_github_repo}/README.md" ]; then
            echo "Incorrect path to github repo provided: ${path_to_github_repo}"
            exit 1
fi


echo "working directory: $set_working_directory"
echo "input directory: $input_directory"
echo "input file name: $input_file_name"
echo "github repository path: $path_to_github_repo"
echo "user: $user_x500"
echo "crossmap: $use_crossmap"
echo "genome harmonizer: $use_genome_harmonizer"
echo "king: $use_king"
echo "genome rfmix: $use_rfmix"
echo "make report: $make_report"
echo "custom qc: $custom_qc"

if [ ${custom_qc} -eq 1 ]; then
          echo "You have chosen to customize the standard qc steps, please answer all of the following questions"
          source ${path_github_repo}/src/gather_custom_qc_inputs.sh ${path_github_repo} ${desired_working_directory}
fi

#source /home/gdc/and02709/QCmja/create_main_wrapper.sh \
source ${path_to_github_repo}/src/create_main_wrapper.sh \
${input_directory} \
${input_file_name} \
${path_to_github_repo} \
${user_x500} \
${set_working_directory} \
${use_crossmap} \
${use_genome_harmonizer} \
${use_king} \
${use_rfmix} \
${make_report} \
${custom_qc}

sleep 0.5

#cp /home/gdc/and02709/QCmja/temp.sh /home/gdc/and02709/QCmja/SMILES_GDA_folder/temp.sh
source ${set_working_directory}/${input_file_name}_wrapper.sh

exit 0
