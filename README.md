# Dawdlib
Create a cost sensitive degenerate codon sequence

## Requirments
* python >= 3.6

## Install

1. clone the repo: `git clone https://github.com/Fleishman-Lab/dawdlib.git`
2. Make sure the rust compiler is available with `rustc --version`
- In case it's not you try `module load rust` or install following: [How to install Rust](https://www.rust-lang.org/tools/install) 
3. Install using pip: `pip install dawdlib/`

## Usage
Please follow the example notebook in example/gg_oligo_design.ipynb
If you have any questions you know where to find me :)

*General instructions*
1. set the variable GG_temp in box 1.2 to 25 in order to follow the GGA reaction conditions in NEB's paper (ADD_REF) 
2. set the following variables in box 1.2:
a. W_PATH to be the working dir
b. resfile_path to your input refile. Make sure that your file does not contain positions with a single option (they are treated as variable positions)
c. dna_path to the sequence file. It should contain DNA seq of a parent design of your repertoire in fasta format. 
Do not forget to add lack of density before translating the protein back to DNA.
