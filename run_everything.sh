# run simulations
python simulate_mutations.py data/ndl.fasta data/
python simulate_mutations.py data/staphylococcus.fasta data/

# compute rates
python estimate_rates.py data/ndl.fasta data/ 128 est_rates_ndl.csv
python estimate_rates.py data/staphylococcus.fasta data/ 128 est_rates_staph.csv

# plot