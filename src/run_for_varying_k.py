import os
import subprocess
import argparse
import numpy as np
from Bio import SeqIO
from numpy.linalg import solve
from tqdm import tqdm
from multiprocessing import Pool
import edlib



def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in reversed(seq)])


def run_cuttlefish(genome_filename, k, num_threads, outoput_prefix):
    #rm random_mutated.fasta_unitigs*
    #cuttlefish build -s random_mutated.fasta -k 21 -t 128 -o random_mutated.fasta_unitigs -w . --ref
    
    cmd = f"rm {outoput_prefix}*"
    # run the command
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    
    # run cuttlefish
    cmd = f"cuttlefish build -s {genome_filename} -k {k} -t {num_threads} -o {outoput_prefix} -w . --ref"
    print(cmd)
    
    with open(os.devnull, 'wb') as devnull:
        subprocess.check_call(cmd.split(' '), stdout=devnull, stderr=subprocess.STDOUT)


def split_unitigs(unitigs, k):
    return_list = []
    for u in unitigs:
        if len(u) < 10000:
            return_list.append(u)
        else:
            num_splits = len(u) / 8000
            # round up
            num_splits = int(num_splits) + 1
            for i in range(num_splits):
                start = i * 8000
                end = min((i+1) * 8000, len(u))
                if start >= end:
                    break
                return_list.append(u[start:end])
                start = end - (k-1)
                end = min(end + k, len(u))
                if start >= end:
                    break
                return_list.append(u[start:end])
    return return_list




def read_unitigs(unitigs_file):
    unitigs = set()
    with open(unitigs_file) as f:
        for line in f:
            if line[0] == '>':
                continue
            else:
                unitigs.add(line.strip())
    return list(unitigs)


def compute_S_D_I_N(u1, unitig_set_mutd, k):
    num_kmers_single_subst, num_kmers_single_delt, num_kmers_no_mutation = 0, 0, 0
    num_kmers_single_insertion = 0

    for u2 in unitig_set_mutd:
        alignment, distance, st1, st2 = None, 9999999999, None, None
        
        r1 = edlib.align(u1, u2, mode = "HW", task = "path")
        r2 = edlib.align(u2, u1, mode = "HW", task = "path")
        
        u3 = reverse_complement(u1)
        r3 = edlib.align(u3, u2, mode = "HW", task = "path")
        r4 = edlib.align(u2, u3, mode = "HW", task = "path")
        
        for i, r in enumerate([r1, r2, r3, r4]):
            if r['editDistance'] < distance:
                alignment, distance = r, r['editDistance']
                if i == 0:
                    st1, st2 = u1, u2
                    flip = False
                elif i == 1:
                    st1, st2 = u2, u1
                    flip = True
                elif i == 2:
                    st1, st2 = u3, u2
                    flip = False
                else:
                    st1, st2 = u2, u3
                    flip = True
        
        nice = edlib.getNiceAlignment(alignment, st1, st2)
        seqA, seqB = nice['query_aligned'], nice['target_aligned']
        assert len(seqA) == len(seqB)
        
        if flip:
            seqB, seqA = seqA, seqB
            
        alphabet = set('ACGT')
        num_chars = len(seqA)
        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqA[i] != seqB[i]:
                if seqA[i] in alphabet and seqB[i] in alphabet:
                    in_numbers[i] = 1
                else:
                    in_numbers[i] = 2

        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                num_kmers_single_subst += 1

        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqB[i] == '-' and seqA[i] in alphabet:
                in_numbers[i] = 1
            elif seqA[i] != seqB[i]:
                in_numbers[i] = 2

        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                num_kmers_single_delt += 1
            if sum(in_numbers[i:i+k]) == 0:
                num_kmers_no_mutation += 1
        
        in_numbers = [0 for i in range(num_chars)]
        for i in range(num_chars):
            if seqB[i] in alphabet and seqA[i] == '-':
                in_numbers[i] = 1
            elif seqA[i] != seqB[i]:
                in_numbers[i] = 2

        for i in range(num_chars-k+1):
            if sum(in_numbers[i:i+k]) == 1:
                num_kmers_single_insertion += 1
                    
                    
    return num_kmers_single_subst, num_kmers_single_delt, num_kmers_single_insertion, num_kmers_no_mutation

def wrapper(args):
    return compute_S_D_I_N(*args)

def compute_S_D_I_N_all(unitig_set_orig, unitig_set_mutd, k, num_threads=64):
    arg_list = [(u1, unitig_set_mutd, k) for u1 in unitig_set_orig]
    
    S, D, I, N = 0, 0, 0, 0

    with Pool(num_threads) as pool:
        for result in tqdm(pool.imap(wrapper, arg_list), total=len(arg_list), desc="Processing"):
            S_, D_, I_, N_ = result
            S += S_
            D += D_
            I += I_
            N += N_

    return S, D, I, N


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


def compute_mutation_rates(genome_filename1, genome_filename2, k, num_threads = 255):
    orig_string = read_genome(genome_filename1)
    mutated_string = read_genome(genome_filename2)
    
    L = len(orig_string)
    L2 = len(mutated_string)
    
    fA = orig_string.count('A')
    fA_mut = mutated_string.count('A')
    
    genome1_cuttlefish_prefix = genome_filename1+f"_{k}_"+"_unitigs"
    genome1_unitigs_filename = genome1_cuttlefish_prefix + ".fa"
    genome2_cuttlefish_prefix = genome_filename2+"_unitigs"
    genome2_unitigs_filename = genome2_cuttlefish_prefix + ".fa"
    
    run_cuttlefish(genome_filename1, k, num_threads, genome1_cuttlefish_prefix)
    run_cuttlefish(genome_filename2, k, num_threads, genome2_cuttlefish_prefix)
    
    assert os.path.exists(genome1_unitigs_filename), f"Mutated unitigs file {genome1_unitigs_filename} not found"
    assert os.path.exists(genome2_unitigs_filename), f"Original unitigs file {genome2_unitigs_filename} not found"
    
    # read two sets of unitigs
    unitig_set_orig = read_unitigs(genome1_unitigs_filename)
    unitig_set_mutd = read_unitigs(genome2_unitigs_filename)
    
    # split unitigs into smaller unitigs
    unitig_set_orig = split_unitigs(unitig_set_orig, k)
    unitig_set_mutd = split_unitigs(unitig_set_mutd, k)
    
    # compute S, D, I, N
    S, D, I, N = compute_S_D_I_N_all(unitig_set_orig, unitig_set_mutd, k, num_threads)
    
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
    ksizes = [15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51]
    ps, pd, d = args.rate, args.rate, args.rate
    L = 1000000
    fA = 0.3
    
    for i in range(args.num_simulations):
        for ksize in ksizes:
            
            # Create random genome
            random_genome_filename = os.path.join(args.working_dir, f"random_genome_k_{ksize}_sim_{i}.fasta")
            seed = i
            create_random_genome(fA, L, random_genome_filename, seed)
            
            # Create mutated genome
            mutated_genome_filename = os.path.join(args.working_dir, f"random_mutated_k_{ksize}_sim_{i}.fasta")
            create_mutated_genome(random_genome_filename, mutated_genome_filename, ps, pd, d, ksize, seed)
            
            # Compute mutation rates
            subst_rate, del_rate, ins_rate = compute_mutation_rates(random_genome_filename, mutated_genome_filename, ksize, num_threads=128)
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
        f.write("simulation\tksize\tsubst_rate\tdel_rate\tins_rate\n")
    print(f"Output file created: {args.output_file}")
    
    # Run the simulations
    run_simulations_est_scores_and_record(args)
    return

if __name__ == "__main__":
    main()