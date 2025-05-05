import argparse
import random
import string

# arguments: fA, L, output_filename
def parse_arguments():
    parser = argparse.ArgumentParser(description="Create a random genome.")
    parser.add_argument("fA", type=float, help="Frequency of A in the genome.")
    parser.add_argument("L", type=int, help="Length of the genome.")
    parser.add_argument("output_filename", type=str, help="Output filename for the random genome.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    
    return parser.parse_args()

def create_random_genome(fA, L, output_filename, seed=None):
    if seed is not None:
        random.seed(seed)  # Set the random seed for reproducibility
    if not (0 <= fA <= 1):
        raise ValueError("Frequency of A (fA) must be between 0 and 1.")
    if L <= 0:
        raise ValueError("Length of the genome (L) must be a positive integer.")
    if not output_filename.endswith('.fasta'):
        raise ValueError("Output filename must end with '.fasta'.")
    if not output_filename:
        raise ValueError("Output filename cannot be empty.")
    if not isinstance(output_filename, str):
        raise ValueError("Output filename must be a string.")
    if not isinstance(L, int):
        raise ValueError("Length of the genome (L) must be an integer.")
    if not isinstance(fA, float):
        raise ValueError("Frequency of A (fA) must be a float.")
    if fA < 0 or fA > 1:
        raise ValueError("Frequency of A (fA) must be between 0 and 1.")
    if L <= 0:
        raise ValueError("Length of the genome (L) must be a positive integer.")
    if not isinstance(seed, (int, type(None))):
        raise ValueError("Seed must be an integer or None.")
    if seed is not None and (seed < 0 or seed > 2**32 - 1):
        raise ValueError("Seed must be a non-negative integer less than 2^32.")
    if not isinstance(output_filename, str):
        raise ValueError("Output filename must be a string.")
    # Create a random genome with the given frequency of A and length L
    
    # Generate a random genome
    genome = ''.join(random.choices(['A', 'C', 'G', 'T'], weights=[fA, (1-fA)/3, (1-fA)/3, (1-fA)/3], k=L))
    
    # Write the genome to the output file
    with open(output_filename, 'w') as f:
        f.write(f">{output_filename}\n")
        f.write(genome + "\n")
    
    print(f"Random genome created and saved to {output_filename}.")
    return



if __name__ == "__main__":
    args = parse_arguments()
    
    # Create the random genome
    create_random_genome(args.fA, args.L, args.output_filename, args.seed)
    print(f"Random genome created and saved to {args.output_filename}.")
    # Example usage: python create_random_genome.py 0.2 1000 random_genome.fasta --seed 0