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

# Some commands
```
python src/run_for_varying_fAs.py --L 1000000 --fA_min 0.15 --fA_max 0.35 --fA_step 0.01 --num_simulations 10 --output_file results/estimated_rates_by_varying_fAs.csv --working_dir data_for_varying_fA --rate 0.01

python src/plot_results_varying_fAs.py --estimated_rates_file results/estimated_rates_by_varying_fAs_rates_set_to_0.05.csv --output_filename plots/random_estimated_rates_varying_fA_true_rates_set_to
_0.05.pdf

python src/plot_results_varying_fAs.py --estimated_rates_file results/estimated_rates_by_varying_fAs_rates_set_to_0.01.csv --output_filename plots/random_estimated_rates_varying_fA_true_rates_set_to
_0.01.pdf
```