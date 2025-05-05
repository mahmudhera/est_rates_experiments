# experiment by varying fAs


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