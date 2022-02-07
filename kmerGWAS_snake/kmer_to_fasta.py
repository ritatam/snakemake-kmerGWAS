#!/g/data/xf3/ht5438/miniconda3/envs/kmerGWAS/bin/python
# convert kmers passing 5% threshold to FASTA format

import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="This function appends and converts the fetched kmers into a FASTA file.")
parser.add_argument("-n", "--Yr_name", metavar="", required=True, help="Name of the Yr strain.")
parser.add_argument("-d", "--date", metavar="", required=True, help="Name of the datetime directory containing kmerGWAS output.")
parser.add_argument("-f", "--kmer_filename", metavar="", required=True, help="Absolute path to pass_threshold_5per.")
parser.add_argument("-k", "--kmer_length", metavar="", default="31", help="kmer length. default = 31.")
parser.add_argument("-o", "--outdir", metavar="", required=True, help="Path to desired output directory.")
args = parser.parse_args()


def kmers_to_fasta(Yr_name, date, kmer_fn, k="31"):
    with open(kmer_fn, "r") as kmers:
        kmers.readline()
        kmerlist = [line[2:int(k)+2] for line in kmers.readlines()]
    kmers_fasta_dir = os.path.join(args.outdir, f"{date}_analysis", Yr_name, "kmers_fasta")
    os.makedirs(kmers_fasta_dir, exist_ok=True)
    with open(os.path.join(kmers_fasta_dir, f"{Yr_name}_kmers.fasta"), "w") as fasta_file:
        for n in range(1,len(kmerlist)+1):
            print(f">{date}_{Yr_name}_kmer{n}\n{kmerlist[n-1]}", file=fasta_file)


if __name__ == "__main__":
    kmers_to_fasta(args.Yr_name, args.date, args.kmer_filename, args.kmer_length)
    
    