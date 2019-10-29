"""
a notebook for finding the best gates and degenerate codons for a given library
"""
#%%
import numpy as np
#%%
import sys
#
import os
from gate_data import GGData
from typing import List

#%% [markdown]
# setup all input

#%%
W_PATH = ""
dna_seq = parse_dnaseqs()
resfile = ResFile()

#%% [markdown]
# this is where the degenerate codon stuff comes in

#%%
def get_all_self_binding_gates() -> List[str]:
    gg_data = GGData()
    all_gates = gg_data.lig_df.columns
    return [g for g in all_gates if gg_data.score_gate1_rc_gate2(g, g) > 2000]