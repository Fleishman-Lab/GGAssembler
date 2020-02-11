OLIGO_PREFIX = "CGTGCGGTCTCG"
OLIGO_SUFFIX = "CGAGACCGCGCCGGGC"

# Variable oligo addition cost
VAR_ADD_COST = len(OLIGO_PREFIX) + len(OLIGO_SUFFIX)

# Constant region require primers to amplify
# therefore the incur a constant cost of a pair of primers
CONST_COST = 40

__all__ = ['OLIGO_PREFIX', 'OLIGO_SUFFIX', 'CONST_COST', 'VAR_ADD_COST']