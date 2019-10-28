#!/usr/bin/env python3
"""
a notebook for finding the best gates and degenerate codons for a given library
"""
#%%
import sys
#
import os
from flab.rosetta.rosetta_output.resfile import ResFile
from flab.sequences.dnaseq import DNASeq, parse_dnaseqs

#%% [markdown]
# setup all input

#%%
W_PATH = ""
dna_seq = parse_dnaseqs()
resfile = ResFile()

#%% [markdown]
# this is where the degenerate codon stuff comes in

#%%
