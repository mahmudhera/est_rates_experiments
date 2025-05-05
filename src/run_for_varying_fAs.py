# workflow
# set L = 1M
# fA = 0.15 to 0.35, varies by 0.01
# num_simulation for each setting is 10
# create random genome
# for each random genome, create a mutated genome using ps, pd, d = 0.01
# estimate rates
# record estimates rates against fA


import os
import subprocess
import argparse
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Run simulations for varying fA.")
    parser.add_argument("--L", type=int, default=1000000, help="Length of the genome.")
    parser.add_argument("--fA_min", type=float, default=0.15, help="Minimum frequency of A.")
    parser.add_argument("--fA_max", type=float, default=0.35, help="Maximum frequency of A.")
    parser.add_argument("--fA_step", type=float, default=0.01, help="Step size for frequency of A.")
    parser.add_argument("--num_simulations", type=int, default=10, help="Number of simulations for each setting.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    parser.add_argument("--output_file", type=str, default="results_by_varying_fA", help="Where to save the results.")
    parser.add_argument("--working_dir", type=str, default="data/", help="Working directory.")
    
    return parser.parse_args()


def create_random_genome(fA, L, output_filename, seed):
    """
    Create a random genome with the given frequency of A, length L, and seed
    """
    # Create the random genome
    subprocess.run(["python", "src/create_random_genome.py", str(fA), str(L), output_filename, "--seed", str(seed)], check=True)
    
    
def create_mutated_genome(input_filename, output_filename, ps, pd, d, ksize, seed):
    """
    Create a mutated genome using the input genome and specified mutation parameters.
    """
    # command: f"mutate_genome {seed} {input_file} {ps} {pd} {d} {output_filename} {ksize}"
    # Create the mutated genome
    cmd = f"mutate_genome {seed} {input_filename} {ps} {pd} {d} {output_filename} {ksize}"
    subprocess.run(cmd, shell=True, check=True)
    
    
def run_simulations_est_scores_and_record(args):
    ksizes = [21, 31, 41]
    ps, pd, d = 0.01, 0.01, 0.01
    
    for fA in [round(fA, 2) for fA in list(np.arange(args.fA_min, args.fA_max + args.fA_step, args.fA_step))]:
        for i in range(args.num_simulations):
            # Create random genome
            random_genome_filename = os.path.join(args.working_dir, f"random_genome_fA_{fA}_sim_{i}.txt")
            create_random_genome(fA, args.L, random_genome_filename, args.seed)
            
            # Create mutated genome
            mutated_genome_filename = os.path.join(args.working_dir, f"mutated_genome_fA_{fA}_sim_{i}.txt")
            create_mutated_genome(random_genome_filename, mutated_genome_filename, ps, pd, d, ksizes[0], args.seed)
            
            # Run simulations and estimate scores
            # Assuming there's a function or script to run simulations and estimate scores
            # For example:
            # estimate_scores(mutated_genome_filename)
            
            # Record results
            # Assuming we have a function to record results
            # record_results(fA, mutated_genome_filename)
    