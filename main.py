#! /bin/python3

from argparse import ArgumentParser
import sys
import os
import pandas as pd
from Bio import SeqIO
import time
import glob
import re




"""
    SelectionPy is a python application to produce a selection pressure analysis (dn/ds) 
    from multiple genomes.
    Author: Metz Sebastian
    email: seba.metz91 [at] gmail.com
    Day: 7/03/2023
    Reference: 
    GitHub: github/sebametz/SelectionPy
"""

def arg_parser():
    parser = ArgumentParser(
        prog = "SelectionPy",
        # usage = "%(prog)s [options]",
        description = "SelectionPy is a python application to produce a selection pressure analysis (dn/ds)",
        epilog = "if any problem please contact me: seba.metz91 [at] gmail.com",
    )

    parser.add_argument(
        "--input_dir", "-i",
        help = "Working directory with the genomes in FASTA format with extension .fasta with their respectively gff files",
        action = "store",
        dest = "input"
    )

    parser.add_argument(
        "--reference", "-r",
        help = "Input name of the reference organisms",
        action = "store",
        dest = "reference"
    )

    parser.add_argument(
        "--evalue", "-e",
        help = "Blastp min E-value",
        action = "store",
        dest = "evalue",
        default = 1e-6,
        type = float,
    )

    parser.add_argument(
        "--threads", "-N",
        help = "number of threads used",
        action = "store",
        dest = "n_cores",
        default = 1,
        type = int,
    )

    parser.add_argument(
        "--version", "-v", 
        action='version', 
        version="%(prog)s version 1.0.0",)

    
    args = parser.parse_args()
    return args


# GLOBAL VARIABLES

INPUT = arg_parser().input
CORES = arg_parser().n_cores
REFERENCE = arg_parser().reference

   
# run a cmd and exit if status is different to 0
def run_cmd(cmd, menssage = None):
    res = os.system(cmd)
    time.sleep(1)
    if menssage != None:
        print(f"[{menssage}] finished with code({res})")
    if res != 0:
        return sys.exit(0)

# set working directory
def set_working_directory():
    cmd = f'mkdir -p {os.path.join(INPUT, "proteomes")} {os.path.join(INPUT, "genes")} {os.path.join(INPUT, "tmp")} \
        {os.path.join(INPUT, "tmp/blastdb")} {os.path.join(INPUT, "tmp/BRH")} {os.path.join(INPUT, "tmp/clusters")} {os.path.join(INPUT, "tmp/alignments")} \
        {os.path.join(INPUT, "tmp/codeml")} {os.path.join(INPUT, "tmp/kaks")} {os.path.join(INPUT, "reports")}'
    run_cmd(cmd, "set working directory")

# Parse species from genomes files, test gff gtf files number is iqual to fasta
def parse_files():
    if os.path.exists(INPUT):
        genomes = {}
        gffs = {} 
        for i in os.listdir(INPUT):
            genome = i.split(".fasta")
            gff = i.split(".gff")
            if len(genome) > 1:
                genomes[genome[0]] = os.path.join(INPUT, i)
            elif len(gff) > 1:
                gffs[gff[0]] = os.path.join(INPUT, i)
        if len(genomes) == len(gffs):
            for k,v in gffs.items():
                # get transcripts
                tmp = os.path.join(INPUT, f"genes/{k}.fna")
                cmd = f"gffread -w {tmp} -g {genomes[k]} {v}"
                run_cmd(cmd, f"{k} gff to transcript")
                # get proteins
                tmp = os.path.join(INPUT, f"proteomes/{k}.faa")
                cmd = f"gffread -y {tmp} -g {genomes[k]} {v}"
                run_cmd(cmd, f"{k} gff to proteins")
            print("[parse files] finished with code(0)")
            return gffs.keys()
        else:
            print(f"{INPUT} FASTA and GTF length doesn't match: Exit!")
            return sys.exit(0)    
    else:
        print(f"{INPUT} direcotry not found: Exit!")
        return sys.exit(0)

# make blast dbs in temporal folder tmp/blastdb/specie
def makeblastdbs(species):
    for k in species.keys():
        tmp = os.path.join(INPUT, f"tmp/blastdb/{k}")
        prot = os.path.join(INPUT, f"proteomes/{k}.faa")
        cmd =  f"makeblastdb -in {prot} -dbtype prot -parse_seqids -out {tmp}"
        run_cmd(cmd, f"make blast db {k}")


# run blast
def run_blast(query, subject, evalue = 1e-6, num_hits = 1000, ncores = CORES):

    db = os.path.join(INPUT, f"tmp/blastdb/{subject}")
    query_prot = os.path.join(INPUT, f"proteomes/{query}.faa")
    tmp = os.path.join(INPUT, f"tmp/BRH/{query}_vs_{subject}.blastp")

    cmd = f'blastp -db {db} -query {query_prot} -out {tmp} -outfmt "6 std" \
        -evalue {evalue} -num_threads {ncores} -max_target_seqs {num_hits}'
    
    run_cmd(cmd, f"blastp {query} vs {subject}")

# Get best hit from the blastp result
def get_besthits(query, subject):
    
    file = os.path.join(INPUT, f"tmp/BRH/{query}_vs_{subject}.blastp")

    # get best hit from results of reference vs specie
    results = pd.read_table(file, header=None, names = ["query", "subject", "pident", "len", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    sorted = results.sort_values(by=["evalue", "bitscore"], ascending=True)
    filter = sorted.groupby("query").first().reset_index()
    
    # save results in a tmp output file

    out = os.path.join(INPUT, f"tmp/BRH/{query}_vs_{subject}.fblastp")
    filter.to_csv(f"{out}", sep="\t")
    
    print(f"[{out} created] exit with code(0)")

# # Execution of reciprocal best hits between reference and specie x
def get_BHR(specie1, specie2):
    
    # read files
    tmp = os.path.join(INPUT, f"tmp/BRH/{specie1}_vs_{specie2}.fblastp")
    df1 = pd.read_table(tmp)

    tmp = os.path.join(INPUT, f"tmp/BRH/{specie2}_vs_{specie1}.fblastp")
    df2 = pd.read_table(tmp)
    
    # filter only reciprocal hits
    newdf2 = df2[["subject","query"]].rename(columns = {"query":"subject", "subject":"query"})
    newdf1 = df1[["query","subject"]]
    final = pd.merge(newdf1, newdf2, how="inner")

    final = final.assign(reference = specie1, specie = specie2)
    print(f"[BHR between {specie1} and {specie2}] finished with code(0)")
    return final

# get clusters from BHR file VERY SLOW! NEED OPTIMIZATION
def get_clusters():

    # read BHR file
    tmp = os.path.join(INPUT, f"tmp/BRH/{REFERENCE}.bhr")
    file = open(tmp, "r")
    file.readline() # ignore first line

    # for each line extract variables
    for line in file:

        splite_line = line.split("\t")
        cluster = splite_line[0] # protein name
        sequences = splite_line[1] # id subject
        organisms = splite_line[3].split("\n")[0] # name subject
    
        
        clust_prot = os.path.join(INPUT, f"tmp/clusters/{cluster}.faa")
        clust_nuc = os.path.join(INPUT, f"tmp/clusters/{cluster}.fna")
        
        # if file not exist add cluster id sequences
        if(not os.path.exists(clust_prot)):

            # # get sequences from reference
            # protein sequence
            prot = os.path.join(INPUT, f"proteomes/{REFERENCE}.faa")
            seq_dict = SeqIO.to_dict(SeqIO.parse(prot, "fasta"))
            with open(clust_prot, "a") as output_handle:
                SeqIO.write(seq_dict[cluster], output_handle, "fasta")

            # nuc sequence
            nuc = os.path.join(INPUT, f"genes/{REFERENCE}.fna")
            seq_dict = SeqIO.to_dict(SeqIO.parse(nuc, "fasta"))
            with open(clust_nuc, "a") as output_handle:
                SeqIO.write(seq_dict[cluster], output_handle, "fasta")


        # # get sequences from organisms
        # protein sequence
        prot = os.path.join(INPUT, f"proteomes/{organisms}.faa")
        seq_dict = SeqIO.to_dict(SeqIO.parse(prot, "fasta"))
        with open(clust_prot, "a") as output_handle:
            SeqIO.write(seq_dict[sequences], output_handle, "fasta")
        
        # nuc sequence
        nuc = os.path.join(INPUT, f"genes/{organisms}.fna")
        seq_dict = SeqIO.to_dict(SeqIO.parse(nuc, "fasta"))
        with open(clust_nuc, "a") as output_handle:
            SeqIO.write(seq_dict[sequences], output_handle, "fasta")
        
        print(f"[{cluster}] finished with code(0)")

    print("[get sequences] finished with code(0)")
    file.close()

# function to generate a report
def generate_report():
    path_codeml = os.path.join(INPUT, f"tmp/codeml/*.txt")
    
    
    
    path_report = os.path.join(INPUT, f"reports/{REFERENCE}.report.txt")

    with open(path_report, "w") as out:
        head = "Cluster\tT\tS\tN\tdN/dS\tdN\tdS\n"
        out.write(head)
        for file in glob.glob(path_codeml):
            cluster = file.split(".txt")[0].split("/")[-1]
            with open(file, 'r') as f:
                last_line = f.readlines()[-1]
            values = re.split("\s+",last_line)
            out.write(f"{cluster}\t{values[1]}\t{values[3]}\t{values[5]}\t{values[7]}\t{values[10]}\t{values[13]}\n")

    print(f"[Report generation {cluster}] finished with code (0)")
            

if __name__ == '__main__':
    
    # # set working directories
    set_working_directory() 

    # get dictionary of species and save proteins and genes in their respectively folders
    species = {key: None for key in parse_files()} # these are dictionaries with "specie: path_to_file"
    
    # # # make blast dbs with proteins
    makeblastdbs(species)

    # # # run blasts between Reference and Subject and viceversa
    table_bhr = pd.DataFrame()
    for k in species.keys():
        if k != REFERENCE:
            run_blast(REFERENCE, k, arg_parser().evalue)
            get_besthits(REFERENCE,k)
            run_blast(k, REFERENCE, arg_parser().evalue)
            get_besthits(k,REFERENCE)
            table_bhr = pd.concat([get_BHR(REFERENCE, k), table_bhr])
    
    #sort by query=clusters and save
    table_bhr = table_bhr.sort_values(by="query")
    tmp = os.path.join(INPUT, f"tmp/BRH/{REFERENCE}.bhr")
    table_bhr.to_csv(tmp, sep="\t", index=False)
    
    # read BHR table get sequences and create a list to extract from transcriptome
    get_clusters()

    # for each clusters
    path_clusters = os.path.join(INPUT, "tmp/clusters/*.faa")
    list_faa = glob.glob(path_clusters)
    
    for path_faa in list_faa:
        
        cluster = path_faa.split(".faa")[0].split("/")[-1]
        
        # alignment
        path_aln = os.path.join(INPUT, f"tmp/alignments/{cluster}.aln.faa")
        cmd = f"muscle -align {path_faa} -output {path_aln} -threads {CORES}"
        run_cmd(cmd, f"muscle {cluster} done!")

        # pal2nl
        path_nt = os.path.join(INPUT, f"tmp/clusters/{cluster}.fna")
        path_pl2nl = os.path.join(INPUT, f"tmp/codeml/{cluster}.paml")
        cmd = f"pal2nal.pl {path_aln} {path_nt} -output paml -nogap > {path_pl2nl}"
        run_cmd(cmd, f"pal2nal {cluster} done!")

        # codeml
        # create codml.ctl
        path_codeml = os.path.join(INPUT, f"tmp/codeml/{cluster}.txt")

        if(len(species) > 2):
            # iqtree
            path_tree = os.path.join(INPUT, f"tmp/alignments/{cluster}")
            cmd = f"iqtree -s {path_aln} -m GTR+G -nt AUTO -ntmax {CORES} -alrt 1000 -o {path_tree}"
            run_cmd(cmd, f"Tree {cluster}")
            var = f"""seqfile = {path_pl2nl} \noutfile = {path_codeml}\nmodel = 0\nNSsites = 0\ntreefile = {path_aln}.treefile\n"""
        else:
            var = f"""seqfile = {path_pl2nl} \noutfile = {path_codeml}\nrunmode = -2\nseqtype = 1\nmodel = 0\nNSsites = 0\n"""

        path_codeml_ctl = os.path.join(INPUT, f"tmp/codeml/codeml.ctl")
        ctl = open(path_codeml_ctl, "w")
        ctl.write(var)
        ctl.close()

        # run codeml
        mv = os.path.join(INPUT, f"tmp/codeml/")
        os.chdir(mv)
        cmd = "codeml"
        run_cmd(cmd, f"{cluster} codeml")
        os.chdir(INPUT)

    generate_report()
        
sys.exit(0)




