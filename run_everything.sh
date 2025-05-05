# experiment by varying L


# experiment by varying fAs
python src/run_for_varying_fAs.py --L 1000000 --fA_min 0.15 --fA_max 0.35 --fA_step 0.01 --num_simulations 10 --output_file results/estimated_rates_by_varying_fAs_rates_set_to_0.01.csv --working_dir data_for_varying_fA --rate 0.01
python src/run_for_varying_fAs.py --L 1000000 --fA_min 0.15 --fA_max 0.35 --fA_step 0.01 --num_simulations 10 --output_file results/estimated_rates_by_varying_fAs_rates_set_to_0.05.csv --working_dir data_for_varying_fA --rate 0.05
python src/plot_results_varying_fAs.py --estimated_rates_file results/estimated_rates_by_varying_fAs_rates_set_to_0.01.csv --output_filename plots/random_estimated_rates_varying_fA_true_rates_set_to_0.01.pdf
python src/plot_results_varying_fAs.py --estimated_rates_file results/estimated_rates_by_varying_fAs_rates_set_to_0.05.csv --output_filename plots/random_estimated_rates_varying_fA_true_rates_set_to_0.05.pdf


# run simulations
python src/simulate_mutations.py data/ndl.fasta data/
python src/simulate_mutations.py data/staphylococcus.fasta data/

# compute rates
python estimate_rates.py data/ndl.fasta data/ 128 results/est_rates_ndl.csv
python estimate_rates.py data/staphylococcus.fasta data/ 128 results/est_rates_staph.csv
python estimate_rates.py data/ndl.fasta data/ 128 results/est_rates_ndl_using_known.csv --use_true_values
python estimate_rates.py data/staphylococcus.fasta data/ 128 results/est_rates_staph_using_known.csv --use_true_values

# plot
python src/plot_estimated_rates.py --estimated_rates_file results/est_rates_ndl.csv --output_filename plots/ndl_estimated_rates.pdf --method linear
python src/plot_estimated_rates.py --estimated_rates_file results/est_rates_staph.csv --output_filename plots/staph_estimated_rates.pdf --method linear
python src/plot_estimated_rates.py --estimated_rates_file results/est_rates_ndl_using_known.csv --output_filename plots/ndl_estimated_rates_using_known.pdf --method linear
python src/plot_estimated_rates.py --estimated_rates_file results/est_rates_staph_using_known.csv --output_filename plots/staph_estimated_rates_using_known.pdf --method linear