# libraries
import os, re
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline

# read the sequence accession numbers
viruses = pd.read_csv(os.path.join('..', 'data', 'viruses.csv'), index_col='Accession number')

# specify where the muscle executable is located, and the exact name of the executable
muscle_exe = os.path.join("..", "muscle3.8.31_i86linux64")

# define the command line for muscle
muscle_cline = MuscleCommandline(muscle_exe, input=os.path.join("..", "data", "downloads.fasta"))

# use 2 iterations; when sequences are far apart, the attempt to reach a more finer alignment leads to an error
muscle_cline.maxiters = 2

# run muscle to align all sequences
stdout, stderr = muscle_cline()

# get the alignment
align = AlignIO.read(StringIO(stdout), "fasta")
print(align)

# make sequence lists from the alignment object
mat = [
    [0 if nucleotide == '-' else 1 for nucleotide in rec.seq] for rec in align
]
cmap = ListedColormap(['w', 'r'])
fig, ax = plt.subplots(1, 1)
ax.matshow(mat, cmap=cmap)
ax.set_aspect(ax.get_xlim()[1] / ax.get_ylim()[0] / 3)
ax.set_yticks(range(len(align)))
ax.set_yticklabels(
    [viruses.loc[re.match("(^\S*)(?=\.)", rec.id)[0], 'Virus'] for rec in align],
    fontsize=8
)
ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
plt.tight_layout()

# save for later faster processing or testing
with open(os.path.join('..', 'data', "virus_alignment.fasta"), "w") as alignment_file:
    alignment_file.write(stdout)

print('finished')
