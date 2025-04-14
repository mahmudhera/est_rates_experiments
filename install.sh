g++ src/mutate_genome.cpp -std=c++17 -O3 -o mutate_genome
export PATH=`pwd`:$PATH
conda create -n est_rates --file requirements.txt -c conda-forge -c bioconda -y
conda activate est_rates