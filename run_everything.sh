# experiment by varying k (Fig 1)

# expt
python src/run_for_varying_k.py --num_simulations 20 --output_file results/estimated_rates_by_varying_ks_rates_set_to_0.05.csv --working_dir data_for_varying_k --rate 0.05

# plot
python src/plot_results_varying_k.py --estimated_rates_file results/estimated_rates_by_varying_ks_rates_set_to_0.05.csv --output_filename plots/random_estimated_rates_varying_k_true_rates_set_to_0.05.pdf --true_rate 0.05

# experiment by varying L (Fig 2)

# expt
python src/run_for_varying_L.py --num_simulations 20 --output_file results/estimated_rates_by_varying_Ls_rates_set_to_0.05.csv --working_dir data_for_varying_L --rate 0.05
python src/run_for_varying_L.py --num_simulations 20 --output_file results/estimated_rates_by_varying_Ls_rates_set_to_0.01.csv --working_dir data_for_varying_L --rate 0.01

# plot
python src/plot_results_varying_L.py --estimated_rates_file results/estimated_rates_by_varying_Ls_rates_set_to_0.01.csv --output_filename plots/random_estimated_rates_varying_L_true_rates_set_to_0.01.pdf --true_rate 0.01
python src/plot_results_varying_L.py --estimated_rates_file results/estimated_rates_by_varying_Ls_rates_set_to_0.05.csv --output_filename plots/random_estimated_rates_varying_L_true_rates_set_to_0.05.pdf --true_rate 0.05


# experiment by varying fAs (Fig 3)

# expt
python src/run_for_varying_fAs.py --L 1000000 --fA_min 0.15 --fA_max 0.35 --fA_step 0.01 --num_simulations 10 --output_file results/estimated_rates_by_varying_fAs_rates_set_to_0.01.csv --working_dir data_for_varying_fA --rate 0.01
python src/run_for_varying_fAs.py --L 1000000 --fA_min 0.15 --fA_max 0.35 --fA_step 0.01 --num_simulations 10 --output_file results/estimated_rates_by_varying_fAs_rates_set_to_0.05.csv --working_dir data_for_varying_fA --rate 0.05

# plot
python src/plot_results_varying_fAs.py --estimated_rates_file results/estimated_rates_by_varying_fAs_rates_set_to_0.01.csv --output_filename plots/random_estimated_rates_varying_fA_true_rates_set_to_0.01.pdf --true_rate 0.01
python src/plot_results_varying_fAs.py --estimated_rates_file results/estimated_rates_by_varying_fAs_rates_set_to_0.05.csv --output_filename plots/random_estimated_rates_varying_fA_true_rates_set_to_0.05.pdf --true_rate 0.05


# run simulations for the other plots
python src/simulate_mutations.py data/staphylococcus.fasta data/
python src/simulate_mutations.py data/random.fasta data/

# compute rates (Fig 4, 5, 6)
python estimate_rates.py data/staphylococcus.fasta data/ 128 results/est_rates_staph.csv
python estimate_rates.py data/random.fasta data/ 128 results/est_rates_random_using_known.csv --use_true_values

# plot (Fig 4)
python src/plot_estimated_rates.py --estimated_rates_file results/est_rates_random_using_known.csv --output_filename plots/random_estimated_rates_using_known.pdf --method linear

# plot (Fig 5)
python src/plot_estimated_rates.py --estimated_rates_file results/est_rates_staph.csv --output_filename plots/staph_estimated_rates.pdf --method linear

# plot (Fig 6)
python src/plot_against_smm.py --estimated_rates_file results/est_rates_staph.csv --output_filename plots/staph_against_smm.pdf --method linear