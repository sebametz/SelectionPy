import sys
import os
import time
import glob
import re
import pandas as pd
from Bio import SeqIO


# ------------------- Essential Functions -------------------------------

# run a cmd and exit if the status is different to 0


def run_cmd(cmd, message=None):
    """Run a system command line and return the message of the exit code."""
    res = os.system(cmd)
    time.sleep(1)
    if message is not None:
        print(f"[{message}] finished with code({res})")
        return 0
    if res != 0:
        return sys.exit(0)

# set working directory


def set_working_directory(workdir=os.getcwd()):
    """Set the working directory."""
    cmd = f'mkdir -p {os.path.join(workdir, "proteomes")} \
        {os.path.join(workdir, "genes")} \
        {os.path.join(workdir, "tmp")} \
        {os.path.join(workdir, "tmp/blastdb")} \
        {os.path.join(workdir, "tmp/RBH")} \
        {os.path.join(workdir, "tmp/clusters")} {os.path.join(workdir, "tmp/alignments")} \
        {os.path.join(workdir, "tmp/codeml")} {os.path.join(workdir, "reports")}'
    run_cmd(cmd, "set working directory")


# ------------- MODULE 1: Extract proteins and genes and BLAST sequences -------------------

# Parse species from genomes files, test GFF files number is equal to fasta
def parse_files(workdir=os.getcwd()):
    """Parse genomes and GFF to extract organisms names, proteins and transcripts"""
    if os.path.exists(workdir):
        genomes = {}
        gffs = {}
        for i in os.listdir(workdir):
            genome = i.split(".fasta")
            gff = i.split(".gff")
            if len(genome) > 1:
                genomes[genome[0]] = os.path.join(workdir, i)
            elif len(gff) > 1:
                gffs[gff[0]] = os.path.join(workdir, i)
        if len(genomes) == len(gffs):
            for k, v in gffs.items():
                # get transcripts
                tmp = os.path.join(workdir, f"genes/{k}.fna")
                cmd = f"gffread -w {tmp} -g {genomes[k]} {v}"
                run_cmd(cmd, f"{k} gff to transcript")
                # get proteins
                tmp = os.path.join(workdir, f"proteomes/{k}.faa")
                cmd = f"gffread -y {tmp} -g {genomes[k]} {v}"
                run_cmd(cmd, f"{k} gff to proteins")
            print("[parse files] finished with code(0)")
            return gffs.keys()
        else:
            print(f"{workdir} FASTA and GFF length doesn't match: Exit!")
            return sys.exit(0)
    else:
        print(f"{workdir} direcotry not found: Exit!")
        return sys.exit(0)

# make blast dbs in temporal folder tmp/blastdb/specie


def makeblastdbs(species, workdir=os.getcwd()):
    """Make a blast database for each organism"""
    for k in species.keys():
        tmp = os.path.join(workdir, f"tmp/blastdb/{k}")
        prot = os.path.join(workdir, f"proteomes/{k}.faa")
        cmd = f"makeblastdb -in {prot} -dbtype prot -parse_seqids -out {tmp}"
        run_cmd(cmd, f"make blast db {k}")


# run blast
def run_blast(query, subject, evalue=1e-6, num_hits=1000, workdir=os.getcwd(), ncores=1):
    """Run blast"""
    db = os.path.join(workdir, f"tmp/blastdb/{subject}")
    query_prot = os.path.join(workdir, f"proteomes/{query}.faa")
    tmp = os.path.join(workdir, f"tmp/RBH/{query}_vs_{subject}.blastp")

    cmd = f'blastp -db {db} -query {query_prot} -out {tmp} -outfmt "6 std" \
        -evalue {evalue} -num_threads {ncores} -max_target_seqs {num_hits}'

    run_cmd(cmd, f"blastp {query} vs {subject}")

# Get best hit from the blast result


def get_besthits(query, subject, workdir=os.getcwd()):
    """get best hit from the blast results"""
    file = os.path.join(workdir, f"tmp/RBH/{query}_vs_{subject}.blastp")

    # get best hit from results of reference vs specie
    results = pd.read_table(file, header=None, names=[
                            "query", "subject", "pident", "len", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    sorted = results.sort_values(by=["evalue", "bitscore"], ascending=True)
    filter = sorted.groupby("query").first().reset_index()

    # save results in a tmp output file

    out = os.path.join(workdir, f"tmp/RBH/{query}_vs_{subject}.fblastp")
    filter.to_csv(f"{out}", sep="\t")

    print(f"[{out} created] exit with code(0)")

# ------------- MODULE 2: Extract Reciprocal Best Hits -------------------

# Execution of reciprocal best hits between reference and specie x


def get_rbh(specie1=str(), specie2=str(), workdir=os.getcwd()):
    """Get Reciprocal Best Hits from the blast results"""
    # read files
    tmp = os.path.join(workdir, f"tmp/RBH/{specie1}_vs_{specie2}.fblastp")
    df1 = pd.read_table(tmp)

    tmp = os.path.join(workdir, f"tmp/RBH/{specie2}_vs_{specie1}.fblastp")
    df2 = pd.read_table(tmp)

    # filter only reciprocal hits
    newdf2 = df2[["subject", "query"]].rename(
        columns={"query": "subject", "subject": "query"})
    newdf1 = df1[["query", "subject"]]
    final = pd.merge(newdf1, newdf2, how="inner")

    final = final.assign(reference=specie1, specie=specie2)
    print(f"[RBH between {specie1} and {specie2}] finished with code(0)")
    return final

# get clusters from RBH file (Need Optimization)


def get_clusters(workdir=os.getcwd(), reference=str()):
    """Extract the different clusters"""
    # read RBH file
    tmp = os.path.join(workdir, f"tmp/RBH/{reference}.rbh")
    file = open(tmp, "r")
    file.readline()  # ignore first line

    # for each line extract variables
    for line in file:

        splite_line = line.split("\t")
        cluster = splite_line[0]  # protein name
        sequences = splite_line[1]  # id subject
        organisms = splite_line[3].split("\n")[0]  # name subject

        clust_prot = os.path.join(workdir, f"tmp/clusters/{cluster}.faa")
        clust_nuc = os.path.join(workdir, f"tmp/clusters/{cluster}.fna")

        # if file not exist add cluster id sequences
        if (not os.path.exists(clust_prot)):

            # get sequences from reference
            # protein sequence
            prot = os.path.join(workdir, f"proteomes/{reference}.faa")
            seq_dict = SeqIO.to_dict(SeqIO.parse(prot, "fasta"))
            with open(clust_prot, "a") as output_handle:
                SeqIO.write(seq_dict[cluster], output_handle, "fasta")

            # nuc sequence
            nuc = os.path.join(workdir, f"genes/{reference}.fna")
            seq_dict = SeqIO.to_dict(SeqIO.parse(nuc, "fasta"))
            with open(clust_nuc, "a") as output_handle:
                SeqIO.write(seq_dict[cluster], output_handle, "fasta")

        # get sequences from organisms
        # protein sequence
        prot = os.path.join(workdir, f"proteomes/{organisms}.faa")
        seq_dict = SeqIO.to_dict(SeqIO.parse(prot, "fasta"))
        with open(clust_prot, "a") as output_handle:
            SeqIO.write(seq_dict[sequences], output_handle, "fasta")

        # nuc sequence
        nuc = os.path.join(workdir, f"genes/{organisms}.fna")
        seq_dict = SeqIO.to_dict(SeqIO.parse(nuc, "fasta"))
        with open(clust_nuc, "a") as output_handle:
            SeqIO.write(seq_dict[sequences], output_handle, "fasta")

        print(f"[{cluster}] finished with code(0)")

    print("[get sequences] finished with code(0)")
    file.close()

# ------------- MODULE 3: Genearate reports -------------------
# function to generate a report


def generate_report(workdir=os.getcwd(), reference=str()):
    """Generate final report from the codeml results"""
    path_codeml = os.path.join(workdir, f"tmp/codeml/*.txt")

    path_report = os.path.join(workdir, f"reports/{reference}.report.txt")

    with open(path_report, "w") as out:
        head = "Cluster\tT\tS\tN\tdN/dS\tdN\tdS\n"
        out.write(head)
        for file in glob.glob(path_codeml):
            cluster = file.split(".txt")[0].split("/")[-1]
            with open(file, 'r') as f:
                last_line = f.readlines()[-1]
            values = re.split("\s+", last_line)
            out.write(
                f"{cluster}\t{values[1]}\t{values[3]}\t{values[5]}\t{values[7]}\t{values[10]}\t{values[13]}\n")

    print(f"[Report generation {cluster}] finished with code (0)")
