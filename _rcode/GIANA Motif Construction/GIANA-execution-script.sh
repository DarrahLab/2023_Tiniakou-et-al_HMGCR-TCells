#$ -cwd
#$ -m e
#$ -M agirgis3@jhmi.edu
#$ -l mem_free=16G,h_vmem=16G,h_fsize=150G
#$ -R y
#$ -pe local 4
#$ -o /users/agirgis/job-logs/
#$ -e /users/agirgis/job-logs/

# Run GIANA4.1 on a combined file input 1 which contains NON-unique CDR3 sequences from (a) CD4+ CD154+ anti-HMGCR TCRs (b) Healthy Control PBMCs 
# And also on a combined file input 2 which conatins  NON-unique CDR3 sequences from (a) CD4+ CD154+ anti-CMV TCRs (b) Healthy Control PBMCs  
# Combined input file generated in script 2023_08_28_GIANA-ComboInputFormat.r
# 28 August 2023 AAG 



GIANA_DIR='/users/agirgis/GIANA'
INPUT_DIR='/dcs04/fertig/data/agirgis/ET_myositis'
OUTPUT_DIR='/dcs04/fertig/data/agirgis/GIANA_outputs'

conda activate faiss_1.7.3

for fname in HMGCR CMV; do
	python3 $GIANA_DIR/GIANA4.1.py -f $INPUT_DIR/2023-08-28_GIANA-ComboInput-${fname}.txt -v True -N 4 -M True -o $OUTPUT_DIR
done
