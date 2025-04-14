from Bio import SeqIO
import argparse
import os
import sys
import time
from tqdm import tqdm
import numpy as np
import pandas as pd
from multiprocessing import Pool
from numpy.linalg import solve
import edlib
import subprocess
import random
import math


def run_cuttlefish(genome_filename, k, num_threads, outoput_prefix):
    #rm random_mutated.fasta_unitigs*
    #cuttlefish build -s random_mutated.fasta -k 21 -t 128 -o random_mutated.fasta_unitigs -w . --ref
    
    cmd = f"rm {outoput_prefix}*"
    os.system(cmd)

    cmd = f"cuttlefish build -s {genome_filename} -k {k} -t {num_threads} -o {outoput_prefix} -w . --ref"
    
    with open(os.devnull, 'wb') as devnull:
        subprocess.check_call(cmd.split(' '), stdout=devnull, stderr=subprocess.STDOUT)


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in reversed(seq)])

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

    best_seqA, best_seqB = None, None
    smallest_distance = 9999999999

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


def estimate_rates_polynomial(L, L2, S, D, I, N, k):
    K1 = L - k + 1
    K2 = L2 - k + 1
    
    S_norm = 1.0 * S / (K1 * k)
    D_norm = 1.0 * D / (K1 * k)
    I_norm = 1.0 * I / (K1 * k - K1)
    
    coeffs = [0 for i in range(k+1)]
    coeffs[0] = (S_norm + D_norm + I_norm) * D_norm**k
    coeffs[-1] = 1
    coeffs[-2] = -D_norm
    
    roots = np.polynomial.polynomial.polyroots(coeffs)
    
    p_d_ests = (D_norm - roots)/(S_norm + D_norm + I_norm)
    p_d_ests = [np.real(p_d_est) for p_d_est in p_d_ests if not np.iscomplex(p_d_est)]
    
    if len(p_d_ests) == 0:
        return None, None, None
    
    p_d_ests.sort()
    d_ests = [ (D_norm - (S_norm + D_norm) * p_d_est)/(D_norm - (S_norm + D_norm + I_norm) * p_d_est) - 1.0 for p_d_est in p_d_ests ]
    p_s_ests = [ (S_norm * p_d_est)/(D_norm) for p_d_est in p_d_ests ]
    
    all_solutions = list( zip(p_s_ests, p_d_ests, d_ests) )
    solution_N_ratio = 0.0
    solution = (None, None, None)
    
    for p_s_est, p_d_est, d_est in all_solutions:
        if p_s_est < 0 or p_d_est < 0 or d_est < 0:
            continue
        
        N_using_these = (L - k + 1) * (1 - p_s_est - p_d_est)**k / ( (d_est+1)**(k-1) )
        ratio_this = N_using_these / N
        if ratio_this > 1.0:
            ratio_this = 1.0 / ratio_this
        if ratio_this > solution_N_ratio:
            solution_N_ratio = ratio_this
            solution = (p_s_est, p_d_est, d_est)    
    
    return solution


def estimate_rates_linear(L, L2, N, D, S, fA, fA_mut, k):
    
    if 0.21 <= fA/L <= 0.29:
        a1 = 1
        b1 = -1.0*S/D
        c1 = 0
        d1 = 0
        
    else:
        # use the equations to estimate the rates
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


def split_unitigs(unitigs, k):
    return_list = []
    for u in unitigs:
        if len(u) < 6000:
            return_list.append(u)
        else:
            num_splits = len(u) / 5000
            # round up
            num_splits = int(num_splits) + 1
            for i in range(num_splits):
                start = i * 5000
                end = min((i+1) * 5000, len(u))
                if start >= end:
                    break
                return_list.append(u[start:end])
                start = end - (k-1)
                end = min(end + k, len(u))
                if start >= end:
                    break
                return_list.append(u[start:end])
    return return_list


def compute_mutation_rates(genome_filename1, genome_filename2, k, num_threads = 255):
    orig_string = read_genome(genome_filename1)
    mutated_string = read_genome(genome_filename2)
    
    L = len(orig_string)
    L2 = len(mutated_string)
    
    fA = orig_string.count('A')
    fC = orig_string.count('C')
    fG = orig_string.count('G')
    fT = orig_string.count('T')
    fA_mut = mutated_string.count('A')
    fC_mut = mutated_string.count('C')
    fG_mut = mutated_string.count('G')
    fT_mut = mutated_string.count('T')
    
    genome1_cuttlefish_prefix = genome_filename1+"_unitigs"
    genome1_unitigs_filename = genome1_cuttlefish_prefix + ".fa"
    genome2_cuttlefish_prefix = genome_filename2+"_unitigs"
    genome2_unitigs_filename = genome2_cuttlefish_prefix + ".fa"
    
    run_cuttlefish(genome_filename1, k, 64, genome1_cuttlefish_prefix)
    run_cuttlefish(genome_filename2, k, 64, genome2_cuttlefish_prefix)
    
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
    
    # DEBUG: print L, L2, S, D, I, N, fA, fA_mut, k
    #print(f"DBG: L: {L}, L2: {L2}, S: {S}, D: {D}, I: {I}, N: {N}, fA: {fA}, fA_mut: {fA_mut}, k: {k}")
    # DEBUG: show fA, fC, fG, fT
    #print(f"DBG: fA: {fA}, fC: {fC}, fG: {fG}, fT: {fT}")
    # DEBUG: show fA_mut, fC_mut, fG_mut, fT_mut
    #print(f"DBG: fA_mut: {fA_mut}, fC_mut: {fC_mut}, fG_mut: {fG_mut}, fT_mut: {fT_mut}")
    
    # compute the rates
    subst_rate_lin, del_rate_lin, ins_rate_lin = estimate_rates_linear(L, L2, N, D, S, fA, fA_mut, k)
    subst_rate_poly, del_rate_poly, ins_rate_poly = estimate_rates_polynomial(L, L2, S, D, I, N, k)
    
    return subst_rate_lin, del_rate_lin, ins_rate_lin, subst_rate_poly, del_rate_poly, ins_rate_poly


def compute_subst_rate_smm(genome_filename1, genome_filename2, k):
    # read genome files, clean strings, get kmers
    orig_string = read_genome(genome_filename1)
    mutated_string = read_genome(genome_filename2)
    orig_string = clean_genome_string(orig_string)
    mutated_string = clean_genome_string(mutated_string)
    
    kmers_orig = get_kmers(orig_string, k)
    kmers_mutated = get_kmers(mutated_string, k)
    kmers_orig_set = set(kmers_orig)
    kmers_mutated_set = set(kmers_mutated)
    
    num_intersection = len(kmers_orig_set.intersection(kmers_mutated_set))
    containment1 = num_intersection / len(kmers_orig_set)
    containment2 = num_intersection / len(kmers_mutated_set)
    max_containment = max(containment1, containment2)
    ani = max_containment ** (1.0/k)
    
    return 1.0-ani


def main():
    # arguments: original genome filename, data_directory where all mutated genomes are
    # output: a file. for each combination of ps, pd, d, seed, k: compute rates and write to file
    parser = argparse.ArgumentParser(description="Estimate mutation rates from genome files.")
    parser.add_argument("genome_filename1", type=str, help="Path to the original genome file.")
    parser.add_argument("dir_name", type=str, help="Directory where all mutated genomes are.")
    parser.add_argument("num_threads", type=int, help="Number of threads to use.", default=32)
    parser.add_argument("output_filename", type=str, help="Output filename.")
    args = parser.parse_args()
    
    genome_filename1 = args.genome_filename1
    dir_name = args.dir_name
    num_threads = args.num_threads
    output_filename = args.output_filename
    
    # check if genome_filename1 exists
    if not os.path.exists(genome_filename1):
        print(f"Error: {genome_filename1} does not exist.")
        sys.exit(1)
        
    # check if dir_name exists
    if not os.path.exists(dir_name):
        print(f"Error: {dir_name} does not exist.")
        sys.exit(1)
        
    # check if dir_name is a directory
    if not os.path.isdir(dir_name):
        print(f"Error: {dir_name} is not a directory.")
        sys.exit(1)
        
    # check if genome_filename1 is a fasta file
    if not genome_filename1.endswith(".fasta") and not genome_filename1.endswith(".fa"):
        print(f"Error: {genome_filename1} is not a fasta file.")
        sys.exit(1)
        
    # check if dir_name is empty
    if len(os.listdir(dir_name)) == 0:
        print(f"Error: {dir_name} is empty.")
        sys.exit(1)
        
    rates = [0.01, 0.02, 0.03, 0.04, 0.05]
    ps_rates = list(rates)
    pd_rates = list(rates)
    d_rates = list(rates)
    seeds = list(range(10))
    ksizes = [21, 31, 41, 51]
    
    num_completed = 0
    total = len(ps_rates) * len(pd_rates) * len(d_rates) * len(seeds) * len(ksizes)
    print(f"Total combinations: {total}")
    
    for ps in ps_rates:
        for pd in pd_rates:
            for d in d_rates:
                for seed in seeds:
                    for ksize in ksizes:
                        # show progress
                        num_completed += 1
                        progress = (num_completed / total) * 100
                        print(f"Progress: {progress:.2f}%")
                        
                        # create output filename
                        output_filename = f"{os.path.splitext(os.path.basename(genome_filename1))[0]}_mut_{ps}_{pd}_{d}_{seed}_{ksize}.fasta"
                        output_path = os.path.join(dir_name, output_filename)
                        
                        # check if the file exists
                        if not os.path.exists(output_path):
                            print(f"Error: {output_path} does not exist.")
                            continue
                            
                        # compute rates
                        subst_rate_lin, del_rate_lin, ins_rate_lin, subst_rate_poly, del_rate_poly, ins_rate_poly = compute_mutation_rates(genome_filename1, output_path, ksize, num_threads)
                        
                        # compute subst rate using SMM
                        subst_rate_smm = compute_subst_rate_smm(genome_filename1, output_path, ksize)
                        
                        # write to file
                        # header
                        with open(output_filename, "w") as f:
                            f.write("ps\tpd\td\tseed\tksize\tsubst_rate_lin\tdel_rate_lin\tins_rate_lin\tsubst_rate_poly\tdel_rate_poly\tins_rate_poly\tsubst_rate_smm\n")
                            
                        with open(output_filename, "a") as f:
                            f.write(f"{ps}\t{pd}\t{d}\t{seed}\t{ksize}\t{subst_rate_lin}\t{del_rate_lin}\t{ins_rate_lin}\t{subst_rate_poly}\t{del_rate_poly}\t{ins_rate_poly}\t{subst_rate_smm}\n")
                            
    print(f"Successfully completed!")
    
    
if __name__ == "__main__":
    main()