"""
########################################################
#     SelectionPy is a python application to estimating#
#     the non-synonymous to synonymous rate ratio of   #
#     all protein-coding genes in a genome.            #
########################################################
#     Author: Metz Sebastian
#     email: seba.metz91 [at] gmail.com
#     Created: 7/03/2023
#     Reference:
#     GitHub: github/sebametz/SelectionPy
---------------------------------------------------------
"""

from argparse import ArgumentParser
import sys
import os
import glob
import pandas as pd
from modules import *

def arg_parser():
    """Reading arguments from command line"""
    parser = ArgumentParser(
        prog="SelectionPy",
        # usage = "%(prog)s [options]",
        description="SelectionPy is a python application to estimating \
            the non-synonymous to synonymous rate ratio of all protein-coding \
            genes in a genome.",
        epilog="if any problem please contact me: seba.metz91 [at] gmail.com",
    )

    parser.add_argument(
        "--input_dir", "-i",
        help="Working directory with the genomes in FASTA format with extension \
            .fasta with their respective gff files",
        action="store",
        dest="input"
    )

    parser.add_argument(
        "--reference", "-r",
        help="Input name of the reference genome",
        action="store",
        dest="reference"
    )

    parser.add_argument(
        "--evalue", "-e",
        help="Blastp min E-value",
        action="store",
        dest="evalue",
        default=1e-6,
        type=float,
    )

    parser.add_argument(
        "--threads", "-N",
        help="number of threads used",
        action="store",
        dest="n_cores",
        default=1,
        type=int,
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
EVAL = arg_parser().evalue

if __name__ == '__main__':

    # set working directories
    set_working_directory(workdir=INPUT)

    # get dictionary of species/genomes names
    species = {key: None for key in parse_files(workdir=INPUT)}

    # make blast dbs with proteins
    makeblastdbs(species, workdir=INPUT)

    # run blasts between Reference and Subject and viceversa
    table_rbh = pd.DataFrame()
    for k in species.keys():
        if k != REFERENCE:
            # Reference vs Others
            run_blast(query=REFERENCE, subject=k, evalue=EVAL,
                      workdir=INPUT, ncores=CORES)
            get_besthits(REFERENCE, k, workdir=INPUT)

            # Others vs Reference
            run_blast(query=k, subject=REFERENCE, evalue=EVAL,
                      workdir=INPUT, ncores=CORES)
            get_besthits(query=k, subject=REFERENCE, workdir=INPUT)
            table_rbh = pd.concat(
                [get_rbh(REFERENCE, k, workdir=INPUT), table_rbh])

    # Sort by query=clusters and save
    table_rbh = table_rbh.sort_values(by="query")
    tmp = os.path.join(INPUT, f"tmp/RBH/{REFERENCE}.rbh")
    table_rbh.to_csv(tmp, sep="\t", index=False)

    # Read RBH table get sequences and create a list to extract from transcriptome
    get_clusters(reference=REFERENCE, workdir=INPUT)

    # for each clusters calculate dN/dS
    path_clusters = os.path.join(INPUT, "tmp/clusters/*.faa")
    list_faa = glob.glob(path_clusters)

    for path_faa in list_faa:

        cluster = path_faa.split(".faa")[0].split("/")[-1]

        # alignment
        path_aln = os.path.join(INPUT, f"tmp/alignments/{cluster}.aln.faa")
        command = f"muscle -align {path_faa} -output {path_aln} -threads {CORES}"
        run_cmd(command, f"muscle {cluster} done!")

        # pal2nl
        path_nt = os.path.join(INPUT, f"tmp/clusters/{cluster}.fna")
        path_pl2nl = os.path.join(INPUT, f"tmp/codeml/{cluster}.paml")
        command = f"pal2nal.pl {path_aln} {path_nt} -output paml -nogap > {path_pl2nl}"
        run_cmd(command, f"pal2nal {cluster} done!")

        # codeml
        # create codml.ctl
        path_codeml = os.path.join(INPUT, f"tmp/codeml/{cluster}.txt")

        if (len(species) > 2):
            # iqtree
            path_tree = os.path.join(INPUT, f"tmp/alignments/{cluster}")
            command = f"iqtree -s {path_aln} -m GTR+G -nt AUTO -ntmax {CORES} -alrt 1000 -o {path_tree}"
            run_cmd(command, f"Tree {cluster}")
            var = f"""seqfile = {path_pl2nl} \noutfile = {path_codeml}\nmodel = 0\nNSsites = 0\ntreefile = {path_aln}.treefile\n"""
        else:
            var = f"""seqfile = {path_pl2nl} \noutfile = {path_codeml}\nrunmode = -2\nseqtype = 1\nmodel = 0\nNSsites = 0\n"""

        path_codeml_ctl = os.path.join(INPUT, "tmp/codeml/codeml.ctl")
        ctl = open(path_codeml_ctl, "w")
        ctl.write(var)
        ctl.close()

        # run codeml
        mv = os.path.join(INPUT, "tmp/codeml/")
        os.chdir(mv)
        command = "codeml"
        run_cmd(command, f"{cluster} codeml")
        os.chdir(INPUT)

    # Generate report
    generate_report(workdir=INPUT, reference=REFERENCE)

sys.exit(0)
