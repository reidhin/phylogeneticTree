from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.Alphabet import IUPAC
from itertools import product
import pandas as pd

df = pd.DataFrame()

# Get all permutations of length 3
df['sequence'] = [
    ''.join(perm) for perm in product(['A', 'C', 'G', 'T'], repeat=3)
]

# Get translation
df['amino acid (1 letter)'] = df['sequence'].apply(lambda x: Seq(x, IUPAC.unambiguous_dna).translate())
df['amino acid (3 letter)'] = df['amino acid (1 letter)'].apply(lambda x: seq3(x))

df.groupby(['amino acid (1 letter)', 'amino acid (3 letter)'])['sequence'].apply(', '.join).reset_index()
