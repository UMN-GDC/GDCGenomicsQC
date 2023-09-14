#Load the necessary modules
module load python3/3.9.3_anaconda2021.11_mamba
module load R/4.2.2-openblas

#Download fraposa and reference data
git clone https://github.com/daviddaiweizhang/fraposa.git
cd fraposa
wget https://upenn.app.box.com/v/fraposa-demo/file/1170205195226

#Be sure to have the study data output 9 in the fraposa folder

#Extract common variants between the reference data and the study data
./commvar.sh thousandGenomes output9 thousandG_comm output10
#Run fraposa
./fraposa_runner.py --stu_filepref output10 thousandG_comm
#Predict ancestry for the study data
./predstupopu.py thousandG_comm output10
#Principle component plot that superimpose the study data on the reference data
./plotpcs.py thousandG_comm output10
