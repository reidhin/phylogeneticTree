# libraries
import os, re
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# read the sequence accession numbers
viruses = pd.read_csv(os.path.join('..', 'data', 'viruses.csv'), index_col='Accession number')

# get the alignment
align = AlignIO.read(os.path.join("..", "data", "alignment.fasta"), "fasta")

# calculate distance - https://biopython.org/wiki/Phylo
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(align)

# construct a tree
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)

# remove the names for the non-terminals for better visual appeal
for non_terminal in tree.get_nonterminals():
    non_terminal.name = ''

# change accession numbers into human more understandable names
for terminal in tree.get_terminals():
    terminal.name = viruses.loc[re.match("(^\S*)(?=\.)", terminal.name)[0], 'Virus']

# permutate the leaves of the tree into better order
tree.ladderize()

# plot the tree
fig, ax = plt.subplots(1, 1)
# draw the resulting tree
Phylo.draw(tree, show_confidence=False, axes=ax, do_show=False)
ax.set_xlim(right=1.0)
plt.show()
fig.savefig(os.path.join('..', 'figures', 'tree.png'))

print('finished')
