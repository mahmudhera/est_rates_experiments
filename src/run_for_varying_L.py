# workflow
# fA = 0.3, fC, fG, fT = 0.2, 0.2, 0.2
# vary L: 10K, 3.3K, 100K, 333K, 1M
# take rates as input
# create random genome
# create mutated genome using ps, pd, d = the entered rate
# estimate rates
# record estimates rates against L


import os
import subprocess
import argparse
import numpy as np
from Bio import SeqIO
from numpy.linalg import solve


def read_genome(genome_file):
    """
    Reads a genome file and returns the genome as a string
    """
    genome = ""
    for record in SeqIO.parse(genome_file, "fasta"):
        genome += str(record.seq)
    return clean_genome_string(genome)

def clean_genome_string(genome_string):
    """
    Removes all non-alphabet characters from a genome string
    """
    alphabet = set('ACGT')
    return ''.join(filter(alphabet.__contains__, genome_string))

def get_genome_length(genome_filename):
    """
    Returns the length of the genome in the genome file
    """
    genome = read_genome(genome_filename)
    return len(clean_genome_string(genome))

def get_kmers(genome_string, k):
    """
    Returns a list of all k-mers in a genome string
    """
    kmers = []
    for i in range(len(genome_string)-k+1):
        kmers.append(genome_string[i:i+k])
    return kmers


def estimate_rates_linear(L, L2, N, D, S, fA, fA_mut, k):
    
    a1 = 1.0 * (L - fA) / 3.0 - fA
    b1 = - fA
    c1 = L/4.0
    d1 = fA_mut - fA
    
    a2 = 0
    b2 = -1
    c2 = 1
    d2 = 1.0*L2/L - 1

    a3 = 1
    b3 = 1.0 * N * k / D + 1
    c3 = 0
    d3 = 1.0

    A = np.array([[a1, b1, c1], [a2, b2, c2], [a3, b3, c3]])
    b = np.array([d1, d2, d3])
    x = solve(A, b)

    subst_rate, del_rate, ins_rate = x
    
    #p_s = 3 * ( N*k*(L2-4*fA_mut-L+4*fA) + D*(L2-4*fA_mut) ) / ( (L - 4*fA) * (3 - 4*N*k - 4*D))
    #print(p_s, subst_rate)

    return subst_rate, del_rate, ins_rate


def compute_mutation_rates_by_true_values(genome_filename1, genome_filename2, k):
    orig_string = read_genome(genome_filename1)
    mutated_string = read_genome(genome_filename2)
    
    L = len(orig_string)
    L2 = len(mutated_string)
    
    fA = orig_string.count('A')
    fA_mut = mutated_string.count('A')
    
    # read the mutated file's first line
    S, I, D, N = None, None, None, None
    with open(genome_filename2) as f:
        first_line = f.readline().strip()
        # the line looks like: "> mutated_13110_11715_12488_60632"
        # parse the numbers
        numbers = first_line.split('_')[1:]
        S, D, I, N = [int(x) for x in numbers]
        
    # if S, I, D, N are None, return None
    if S is None or I is None or D is None or N is None:
        print(f"Error: S, I, D, N are None")
        return None, None, None, None, None, None
    
    # compute the rates
    subst_rate_lin, del_rate_lin, ins_rate_lin = estimate_rates_linear(L, L2, N, D, S, fA, fA_mut, k)
    
    return subst_rate_lin, del_rate_lin, ins_rate_lin


def parse_arguments():
    parser = argparse.ArgumentParser(description="Run simulations for varying fA.")
    parser.add_argument("--num_simulations", type=int, default=10, help="Number of simulations for each setting.")
    parser.add_argument("--output_file", type=str, default="results_by_varying_L", help="Where to save the results.")
    parser.add_argument("--working_dir", type=str, default="data/", help="Working directory.")
    parser.add_argument("--rate", type=float, default=0.01, help="Mutation rate, ps, pd, d will all be set to this.")
    
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
    ps, pd, d = args.rate, args.rate, args.rate
    L_values = [10000, 33333, 100000, 333333, 1000000]
    fA = 0.3
    
    for L in L_values:
        for i in range(args.num_simulations):
            for ksize in ksizes:
                
                # Create random genome
                random_genome_filename = os.path.join(args.working_dir, f"random_genome_L_{L}_sim_{i}.fasta")
                seed = i
                create_random_genome(fA, L, random_genome_filename, seed)
                
                # Create mutated genome
                mutated_genome_filename = os.path.join(args.working_dir, f"random_genome_L_{L}_sim_{i}.fasta")
                create_mutated_genome(random_genome_filename, mutated_genome_filename, ps, pd, d, ksize, seed)
                
                # Compute mutation rates
                subst_rate, del_rate, ins_rate = compute_mutation_rates_by_true_values(random_genome_filename, mutated_genome_filename, ksize)
                if subst_rate is None or del_rate is None or ins_rate is None:
                    print(f"Error: Could not compute mutation rates for fA = {fA}, simulation {i}")
                    continue
                
                # Save results
                with open(args.output_file, "a") as f:
                    f.write(f"{L}\t{i}\t{ksize}\t{subst_rate}\t{del_rate}\t{ins_rate}\n")
                        
    print("All simulations completed.")
    return


def main():
    args = parse_arguments()
    # Create the working directory if it doesn't exist
    if not os.path.exists(args.working_dir):
        os.makedirs(args.working_dir)
        
    # Create the output file
    if os.path.exists(args.output_file):
        os.remove(args.output_file)
    with open(args.output_file, "w") as f:
        f.write("L\tsimulation\tksize\tsubst_rate\tdel_rate\tins_rate\n")
    print(f"Output file created: {args.output_file}")
    
    # Run the simulations
    run_simulations_est_scores_and_record(args)
    return

if __name__ == "__main__":
    main()