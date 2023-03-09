<img src = "logo/SelectionPy2.png" alt = "SelectionPy logo" width = 250px>

----------------------------------

### About
SelectionPy is a python application (script for now) to estimate the non-synonymous to synonymous rate ratio (dN/dS) of all protein-coding genes in a genome. SelectionPy aims to allow users to analyse dN/dS quickly and easily. 

### Dependencies
- Python version > 3
- Biopython
- Python Pandas librery
- BLAST 2.13.0+
- Muscle version > 5
- IQ-TREE version 1.6.12
- pal2nal
- codeml

For now, please install the dependencies manually, but I will work to solve it with a __Docker__ container.

### To test it:

    git clone https://github.com/sebametz/SelectionPy.git

    cd SelectionPy

    python selectionPy.py -i test/data -r ProtA.fasta --threads 4

The results of the analysis are in 
    
    /test/reports

### To do list:
I wrote this script in the last couple of days, so there is much to do.
- Test with multiple genomes; I only tested it with two genomes.
- Add analysis with KaKs_Calculator 2.0 to validate CODEML results.
- Optimise configuration files for MUSCLE, IQ-TREE and CODEML.
- Create a user interfase with Tk.
- Create a Docker container to solve dependencies.

### Contact information:

Sebastian Metz: seba.metz91@gmail.com