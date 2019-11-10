# 1
# find the locations of genes in the reference
# ***need to select the top alignments by hand***
echo "run_blast"
bash src/bash/run_blast.sh

# # 2
# # process all the data into intermediate forms
# echo "run_R_data_processing"
# bash src/bash/run_R_data_processing.sh

# # 3
# # use the intermediate data to create tables and figures
# echo "run_R_results"
# bash src/bash/run_R_results.sh