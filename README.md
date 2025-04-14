# Summary
This repository contains the experiments verifying our estimator for the mutation rates works well. The verifications are through simulations.

# Installation
Run the following in the repo's home. Installation requires g++ and conda.
```
g++ src/mutate_genome.cpp -std=c++17 -O3 -o mutate_genome
export PATH=`pwd`:$PATH
conda create -n est_rates --file requirements.txt -c conda-forge -c bioconda -y
conda activate est_rates
```


# Running everything
Run the following in the repo's home. Installation is required.
```
bash run_everything.sh
```