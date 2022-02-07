#!/g/data/xf3/ht5438/miniconda3/envs/kmerGWAS/bin/python

import pandas as pd
import os
import subprocess
import glob
import argparse


parser = argparse.ArgumentParser(description="Filter the Pst accessions for a specified Yr strain and run the kmerGWAS pipeline accordingly. \
                                                Extract the kmers passing the 5% threshold which are automatically converted into a FASTA file \
                                                for downstream analysis.")
parser.add_argument("-d", "--date", metavar="", required=True, help="Date as prefix. YYYYMMDD.")
parser.add_argument("-n", "--host_strain_name", metavar="", required=True, help="Annotations of host strains, e.g. 'Yr8'.")
parser.add_argument("-p", "--pathotype_book", metavar="", required=True, help="Absolute path to pathotype book.")
parser.add_argument("-L", "--library", metavar="", required=True, help="Absolute path to kmerGWAS library.")
parser.add_argument("-O", "--outdir", metavar="", required=True, help="Absolute path to the directory where kmerGWAS output is generated.")
parser.add_argument("-S", "--strand_dir", metavar="", required=True, help="Absolute path to sample strands.")
parser.add_argument("-r", "--refHetfilt", metavar="", default="True", help="(Optional) Filter for PST accessions based on refHetFreqSpec. Can be either True/False. Default = True.")
parser.add_argument("--mac", metavar="", type=int, default=2, help="(Optional) Minor allele count (minimum allowed appearance of a k-mer). Default = 2.")
parser.add_argument("-m", "--min_per", metavar="", type=float, default=0.1, help="(Optional) Minimum percent of appearance in each strand form. Default = 0.1.")
parser.add_argument("-k", "--kmer_length", metavar="", type=int, default=31, help="(Optional) kmer length in bp. Default = 31.")
parser.add_argument("-t", "--threads", metavar="", help="(Optional) Number of threads. Default = 8.")
args = parser.parse_args()




# ========= (I) Quality Control ========= #

def refHetInclusionList(pb_fn, filt='True'):
    """
    Reads in the pathotype book in xlsx format and returns a list of accessions to include
    based on the refHetFreqSpec and the selected filter, which can be either 'True', 'False' or 'NOT SURE'.
    """
    refHetDf = pd.read_excel(pb_fn, sheet_name='refHetFreqSpec', engine='openpyxl')
    refHetDf['usable?'] = refHetDf['usable?'].astype(str)
    return refHetDf[refHetDf['usable?'] == filt]['Identifier'].to_list()


def rm_SRA_not_in_pathobook(pb_fn, refHetInclusionList):
    """
    Pass in the inclusion list of accessions from refHetInclusionList.
    Removes accessions that are not included in the pathotype book from the inclusion list.
    Returns a list of filtered accessions that are present in the pathotype book.
    """
    pathobookDf = pd.read_excel(pb_fn, sheet_name='Pathotype_book', engine='openpyxl')
    return [SRA for SRA in refHetInclusionList if SRA in pathobookDf['Identifier'].tolist()]


def get_phenotype_value(pb_fn, Yr, outpath, refHet_filt='True'):
    """
    Reads in the filtered accession inclusion list based on the selected filter for refHetFreqSpec and pathotype book.
    Takes in a specified Yr strain (as string, e.g. '17' for Yr17) and writes out the corresponding phenotype file in tsv
    format in the output path.
    """
    pathobookDf = pd.read_excel(pb_fn, sheet_name='Pathotype_book', engine='openpyxl').fillna(-1).set_index('Identifier')
    pathobookDf.columns = pathobookDf.columns.astype(str)
    inclusion_list = rm_SRA_not_in_pathobook(pb_fn, refHetInclusionList(pb_fn, refHet_filt))
    pheno_val_list = []
    for SRA in inclusion_list:
        new_row = []
        if pathobookDf.loc[SRA,Yr] != -1.0:
            new_row.append(SRA)
            new_row.append(int(pathobookDf.loc[SRA,Yr]))
            pheno_val_list.append(new_row)
        output_fn = os.path.join(outpath, f"Yr{Yr}_phenotype_values.tsv")
        pd.DataFrame(pheno_val_list, columns=['accession_id','phenotype_value']).to_csv(output_fn,sep='\t',index=False)
        
        
        
# ========= (II) kmerGWAS ========= #       

def create_date_Yr_fn(Yr_name, date, pb_fn, refHet_filt='True'):
    """
    This functionn provides the basic frame of directories for running kmerGWAS.
    Yr_name should be a string that specfies the name of a strain, e.g. 'Yr17'.
    pb_fn is the absolute path to the pathotype book.
    filt refers to filter applied to refHetFreqSpec.
    -
    Creates a directory named as the currrent date.
    Inside the directory, it creates a folder named as the specified Yr name.
    In this folder it creates "{Yr_name}_phenoTable" directory where 
    {Yr_name}_phenotype_values.tsv of the specified Yr gene is stored.
    """
    TODAY_OUTDIR = os.path.join(OUTDIR, date, Yr_name)
    print(TODAY_OUTDIR)
    if not os.path.exists(TODAY_OUTDIR):
        os.makedirs(TODAY_OUTDIR)
    
    os.mkdir(f"{TODAY_OUTDIR}/{Yr_name}_phenoTable")
    get_phenotype_value(pb_fn, Yr_name[2:], f"{TODAY_OUTDIR}/{Yr_name}_phenoTable", refHet_filt)
    return TODAY_OUTDIR


def kmers_list_paths(inclu_list=[]):
    """
    Creates a txt file listing the absolute paths of all sample strands based on the accession inclusion list
    at the datetime directory.
    Each path is followed by the corresponding accession number of the PST strain, separated by tab.
    """
    kmers_wstrand_fns = glob.glob(f"{STRAND_DIR}/*/*.kmers_with_strand")
    sample_list_fn = os.path.join(TODAY_OUTDIR, 'kmers_list_paths.txt')
    kmers_wstrand_fns = [x for x in kmers_wstrand_fns \
                        if any([os.path.basename(x).startswith(y) for y in inclu_list])]
    sample_names = [os.path.basename(l)[:-18] for l in kmers_wstrand_fns]
    join_path_names = [x+'\t'+y for x, y in zip(kmers_wstrand_fns,sample_names)]
    with open(sample_list_fn,'w') as fh:
        for x in join_path_names:
            print(x,file=fh)
    return sample_list_fn


def run_kmer_GWAS(kmers_list_paths, pheno_val_path, LIB_DIR, OUTDIR,\
                  X=2, P=0.1, K=31, THREADS=8):
    """
    Runs the kmerGWAS pipeline. Output generated in directory 'gwas_results' under the datetime directory.
    Generates 4 strings of commands to be run in bash (kmers_to_use, kmers_table, kmers_kinship, and kmerGWAS).
    kmers_to_use, kmers_table and kmers_kinship are generated with Mac{X}P{P}as prefix.
    
    kmers_list_paths: path to kmers_list_paths.txt
    pheno_val_path: path to phenotype file in tsv
    X: minor allele count (minimum allowed appearance of a k-mer).
    P: minimum percent of appearance in each strand form.
    K: kmer length
    THREADS: number of threads
    """
    # paths for output of kmers_to_use and kmers_table
    kmers_to_use = os.path.join(TODAY_OUTDIR, f"Mac{X}P{P}.kmers_to_use")
    kmers_table = os.path.join(TODAY_OUTDIR, f"Mac{X}P{P}.kmers_table")
    
    # --- generate strings of bash commands ---
    
    cmd_kmermultiple = f"{os.path.join(LIB_DIR, 'bin', 'list_kmers_found_in_multiple_samples')}\
    -l {kmers_list_paths} -k {K} --mac {X} -p {P} -o {kmers_to_use}"
    
    cmd_kmer_table = f"{os.path.join(LIB_DIR, 'bin', 'build_kmers_table')}\
    -l {kmers_list_paths} -k {K} -a {kmers_to_use} -o {kmers_table}"

    cmd_kmer_kinship = f"{os.path.join(LIB_DIR, 'bin', 'emma_kinship_kmers')}\
    -t {kmers_table} -k {K} --maf 0.05 > {os.path.join(TODAY_OUTDIR, f'Mac{X}P{P}.kmers_table.kinship')}"
    
    cmd_kmerGWAS = f"python2.7 {os.path.join(LIB_DIR, 'kmers_gwas.py')} --pheno {pheno_val_path}\
    --kmers_table {kmers_table} -l {K} -p {THREADS} --outdir {os.path.join(TODAY_OUTDIR, 'gwas_results')}\
    >> {os.path.join(TODAY_OUTDIR,'gwas_run.log.txt')} 2>&1"
    
    # log all output messages
    kmers_log = os.path.join(TODAY_OUTDIR, f"Mac{X}P{P}.log")
    
    # run all commands in series
    for cmd in [cmd_kmermultiple, cmd_kmer_table, cmd_kmer_kinship, cmd_kmerGWAS]: 
        try:
            print(f"######Running######\n{cmd}\n")
            output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
            print(f"####Done with output#####\n{output}")
            with open(kmers_log, 'a') as fh:
                print(output, file=fh, end='')
            print(f'command logged')
        except:
            print(F"Issues with the following command: \n{cmd}")
            
            
def kmerGWAS_for_Yr(Yr_name, date, pb_fn, LIB_DIR, OUTBASE_DIR, STRAND, refHet_filt='True', X=2, P=0.1, K=31, THREADS=8):
    global STRAND_DIR
    STRAND_DIR = STRAND
    global OUTDIR
    OUTDIR = OUTBASE_DIR
    global TODAY_OUTDIR
    TODAY_OUTDIR = create_date_Yr_fn(Yr_name, date, pb_fn, refHet_filt)
    inclu_list = pd.read_csv(f"{TODAY_OUTDIR}/{Yr_name}_phenoTable/{Yr_name}_phenotype_values.tsv",sep='\t').iloc[:,0].tolist()
    pheno_val_path = f"{TODAY_OUTDIR}/{Yr_name}_phenoTable/{Yr_name}_phenotype_values.tsv"
    run_kmer_GWAS(kmers_list_paths(inclu_list), pheno_val_path, LIB_DIR, OUTDIR, X, P, K, THREADS)
    


if __name__ == "__main__":
    kmerGWAS_for_Yr(args.host_strain_name, args.date, args.pathotype_book, args.library, args.outdir, args.strand_dir, args.refHetfilt, args.mac, args.min_per, args.kmer_length, args.threads)