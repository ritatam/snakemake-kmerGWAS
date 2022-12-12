import os
import pandas as pd
from glob import glob
import argparse


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    
    prepfiles_parser = subparsers.add_parser("PrepFiles", help="Writes phenotype_value.tsv and kmers_list_paths.txt required by downstream kmerGWAS.")
    prepfiles_parser.add_argument("--xlsx", required=True, help="Path to phenotype book spreadsheet (.xlsx).")
    prepfiles_parser.add_argument("--phenotype", required=True, help="Phenotype to analyse.")
    prepfiles_parser.add_argument("--sample_csv", required=True, help="Path to a csv file listing samples to use.")
    prepfiles_parser.add_argument("--strand_dir", required=True, help="Absolute path to directory storing sample subdirectories.")
    prepfiles_parser.add_argument("-o", "--output", required=True, help="Output path to write out files.")
    prepfiles_parser.add_argument("-rf", "--refhetfilt", required=False, action="store_true", help="[Optional] Include this flag to use refHetFreqSpec sheet (if available) in phenotype book spreadsheet to filter samples. Default: False.")
    prepfiles_parser.set_defaults(func=prepare_files)
    
    
    kmers2fasta_parser = subparsers.add_parser("kmers2fasta", help="Converts kmerGWAS output pass_threshold_5per to fasta file for a phenotype.")
    kmers2fasta_parser.add_argument("--kmers_in", required=True, help="Input file pass_threshold_5per")
    kmers2fasta_parser.add_argument("--phenotype", required=True, help="Phenotype that the k-mers significantly associated with.")
    kmers2fasta_parser.add_argument("-o", "--output", required=True, help="Output path to write out kmers FASTA.")
    kmers2fasta_parser.set_defaults(func=kmers2fasta)
    
    args = parser.parse_args()
    
    if args.command == "PrepFiles":
        prepare_files(args.xlsx, args.phenotype, args.sample_csv, args.strand_dir, args.output, args.refhetfilt)
    if args.command == "kmers2fasta":
        kmers2fasta(args.kmers_in, args.phenotype, args.output)



def refHetInclusionList(pb_fn, sample_list):
    """
    Use refHetFreqSpec sheet to filter samples from a input sample list.
    Returns a list of samples that passed the filter.
    """
    refHetDf = pd.read_excel(pb_fn, sheet_name="refHetFreqSpec", engine="openpyxl", dtype=str)
    passed = [x for x in sample_list if x in refHetDf[refHetDf["Passed"] == "True"]["Identifier"].tolist()]
    failed = [x for x in sample_list if x in refHetDf[refHetDf["Passed"] == "False"]["Identifier"].tolist()]
    ambiguous = [x for x in sample_list if x in refHetDf[refHetDf["Passed"] == "AMBIGUOUS"]["Identifier"].tolist()]
    na = [x for x in sample_list if x not in refHetDf["Identifier"].tolist()]
    
    print("\nrefHetFreqSpec filter enabled:")
    print(f"Failed: {len(failed)}")
    print(f"Ambiguous: {len(ambiguous)}")
    print(f"NA: {len(na)}")
    print(f"Passed: {len(passed)}\n")
    
    return passed



def WritePhenotypeFile(pb_fn, phenotype, sample_list, output_path, refHetFreqSpec_filter=False):
    """
    Reads in a sample list, filters based on phenotype availability, 
    and writes out a phenotype value file required by downstream kmerGWAS. 
    
    refHetFreqSpec_filter: If set to be True, it further filters the sample list using the refHetFreqSpec sheet.
    
    Terminates if encounter the following:
    - Fewer than 30 samples after filtering (requirement by kmerGWAS)
    - Only one type of phenotype found across all samples (inappropriate for any GWA)
    """
    
    pbdf = pd.read_excel(pb_fn, sheet_name="Pathotype_book", engine="openpyxl",dtype=str).fillna(int(-1)).set_index("Identifier")
    phenodf = pbdf[phenotype].replace({"FALSE":int(0), "False":int(0), "TRUE":int(1), "True":int(1)}).to_frame(name="Phenotype")
    
    print(f"\nTotal no. of input samples: {len(sample_list)}")
    
    # exclude samples without phenotype data
    phenodf = phenodf[phenodf["Phenotype"] != -1]
    fil_list = [sample for sample in sample_list if sample in phenodf.index]
    print(f"Samples with phenotype data for {phenotype}: {len(fil_list)}")
    
    # filter sample list if refHetFreqSpec filter is enabled
    assert type(refHetFreqSpec_filter) == bool, "Please use either True or False for refHetFreqSpec_filter."
    if refHetFreqSpec_filter == True:
        fil_list = refHetInclusionList(pb_fn, fil_list)
    
    # terminate if filtered sample list contains <30 samples
    assert len(fil_list) > 30, f">>> WARNING: Must have at least 30 passed samples to run kmerGWAS. Terminating. <<<"
    
    phenofile_rows = []
    for sample in fil_list:
        new_row = []
        new_row.append(sample)
        new_row.append(int(phenodf.loc[sample,]))
        phenofile_rows.append(new_row)

    output_fn = os.path.join(output_path, f"{phenotype}_phenotype_values.tsv")
    out_df = pd.DataFrame(phenofile_rows, columns=["accession_id", "phenotype_value"])
    # terminate if only one phenotype across all samples
    assert len(set(out_df["phenotype_value"].tolist())) > 1, f">>> WARNING: Inappropriate to run kmerGWAS for {phenotype} - Only one phenotype across all samples. Terminating. <<<"
    out_df.to_csv(output_fn, sep="\t", index=False)
    print(f"\nPhenotype value file for {len(fil_list)} samples written to: {output_fn}.")
    return out_df["accession_id"]



def kmers_list_paths(sample_list, strand_dir, output_path):
    """
    Writes out kmers_list_paths.txt that lists paths to .kmers_with_strand for all samples in the input sample list.
    """
    kmers_wstrand_fns = glob(os.path.join(strand_dir, "*/*.kmers_with_strand"))
    kmers_wstrand_fns = [x for x in kmers_wstrand_fns \
                        if any([os.path.basename(x).startswith(sample) for sample in sample_list])]
    
    samp = [os.path.basename(x).split(".kmers_with_strand")[0] for x in kmers_wstrand_fns]
    join_path_samp = [x+"\t"+y for x,y in zip(kmers_wstrand_fns, samp)]
    
    kmers_list_paths_fn = os.path.join(output_path, "kmers_list_paths.txt")
    with open(kmers_list_paths_fn, "w") as f:
        for x in join_path_samp:
            print(x, file=f)
    print(f"kmer strand list paths file written to: {kmers_list_paths_fn}\n.")
    return kmers_list_paths_fn


def prepare_files(pb_fn, phenotype, sample_csv, strand_dir, output_path, refHetFreqSpec_filter=False):
    sample_list = pd.read_csv(sample_csv, header=None)[0].tolist()
    fil_list = WritePhenotypeFile(pb_fn, phenotype, sample_list, output_path, refHetFreqSpec_filter)
    kmers_list_paths(fil_list, strand_dir, output_path)



def kmers2fasta(kmer_fn, phenotype, output_path):
    df = pd.read_csv(kmer_fn, sep="\t")
    kmers = [x.split("_")[0] for x in df["rs"].tolist()]
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    with open(os.path.join(output_path, f"{phenotype}_kmers.fasta"), "w") as fasta:
        for i in range(1, len(kmers)+1):
            print(f">{phenotype}_kmer_{i}\n{kmers[i-1]}", file=fasta)

    



if __name__ == "__main__":
    main()
