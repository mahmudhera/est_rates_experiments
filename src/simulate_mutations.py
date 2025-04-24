import argparse 
import os

def parse_arguments():
    # input: one filename: fasta. example: data/name.fasta
    # outputs: many files: data/name_mut_ps_pd_d_seed.fasta, for varying ps, pd, d, seed
    parser = argparse.ArgumentParser(description="Simulate mutations in a FASTA file.")
    parser.add_argument("input_file", type=str, help="Path to the input FASTA file.")
    parser.add_argument("output_dir", type=str, help="Directory to save the output files.")
    
    return parser.parse_args()

def simulate_mutations(input_file, output_dir):
    ps_rates = [0.01, 0.02, 0.03, 0.04, 0.05]
    pd_rates = [0.01, 0.02, 0.03, 0.04, 0.05]
    d_rates = [0.01, 0.02, 0.03, 0.04, 0.05]
    seeds = list(range(10))
    ksizes = [21, 31, 41, 51]
    
    # iterate over all combinations of ps, pd, d, seed, and ksize
    num_completed = 0
    total = len(ps_rates) * len(pd_rates) * len(d_rates) * len(seeds) * len(ksizes)
    for ps in ps_rates:
        for pd in pd_rates:
            for d in d_rates:
                for seed in seeds:
                    for ksize in ksizes:
                        # create output filename
                        output_filename = f"{os.path.splitext(os.path.basename(input_file))[0]}_mut_{ps}_{pd}_{d}_{seed}_{ksize}.fasta"
                        output_path = os.path.join(output_dir, output_filename)
                        
                        # simulate mutations
                        # mutate_genome seed orig_genome_filename p_s p_d d output_filename kmer_size
                        cmd = f"mutate_genome {seed} {input_file} {ps} {pd} {d} {output_path} {ksize}"
                        print(f"Running command: {cmd}")
                        
                        # check if call is successful
                        ret = os.system(cmd)
                        if ret != 0:
                            print(f"Error: Command failed with return code {ret}")
                        else:
                            print(f"Successfully created {output_filename}")
                        
                        # show progress
                        num_completed += 1
                        progress = (num_completed / total) * 100
                        print(f"Progress: {progress:.2f}%\r", end="")
                        
    print("Completed all simulations.")
    return

if __name__ == "__main__":
    args = parse_arguments()
    
    # create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    simulate_mutations(args.input_file, args.output_dir)
    print(f"Simulations completed. Output files are saved in {args.output_dir}.")
    
                        