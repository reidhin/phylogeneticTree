# libraries
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from Bio import Phylo
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyScorer, \
    ParsimonyTreeConstructor, NNITreeSearcher
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO

# read the sequence accession numbers
viruses = pd.read_csv(os.path.join('..', 'data', 'viruses.csv'), index_col='Accession number')

# download the sequences
Entrez.email = 'hans@orikami.nl'  # Always tell NCBI who you are
search = " ".join(viruses.index.to_list())
result = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
seq_records = []
for id in result['IdList']:
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
    seq_records.append(SeqIO.read(handle, "fasta"))

# plot some basic information on the sequences
fig, axs = plt.subplots(1, 2)
axs[0].bar(range(len(seq_records)), [len(seq_record) for seq_record in seq_records])
axs[0].set_ylabel('Length')
axs[1].bar(range(len(seq_records)), [GC(seq_record.seq) for seq_record in seq_records])
axs[1].set_ylabel('GC content')
for ax in axs:
    ax.set_xticks(range(len(seq_records)))
    ax.set_xticklabels(
        [viruses.loc[re.match("(^\S*)(?=\.)", rec.id)[0], 'Virus'] for rec in seq_records],
        rotation=45,
        ha='right'
    )
plt.tight_layout()

# write them in file for later upload
SeqIO.write(seq_records, os.path.join("..", "data", "downloads.fasta"), "fasta")

# run  muscle to align all sequences

# specify where the muscle.exe is located
muscle_exe = os.path.join("..", "muscle.exe")

# define the command line for muscle
muscle_cline = MuscleCommandline(muscle_exe, input=os.path.join("..", "data", "downloads.fasta"))

# use 2 iterations; when sequences are far apart, the attempt to reach a more finer alignment leads to an error
muscle_cline.maxiters = 2

# execute the command
stdout, stderr = muscle_cline()

# get the alignment
align = AlignIO.read(StringIO(stdout), "fasta")
print(align)

# make sequence and id lists from the aln object
mat = [
    [0 if nucleotide == '-' else 1 for nucleotide in rec.seq] for rec in align
]
cmap = ListedColormap(['w', 'r'])
fig, ax = plt.subplots(1, 1)
ax.matshow(mat, cmap=cmap)
ax.set_aspect(ax.get_xlim()[1] / ax.get_ylim()[0] / 3)
ax.set_yticks(range(len(align)))
ax.set_yticklabels(
    [viruses.loc[re.match("(^\S*)(?=\.)", rec.id)[0], 'Virus'] for rec in seq_records],
    fontsize=8
)
ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
plt.tight_layout()

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

print(Phylo.draw_ascii(tree))

# plot the tree
fig, ax = plt.subplots(1, 1)
# draw the resulting tree
Phylo.draw(tree, show_confidence=False, axes=ax, do_show=False)
ax.set_xlim(right=0.8)

print('finished')
