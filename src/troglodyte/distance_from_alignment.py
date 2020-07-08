import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator


# get the alignment
alignments = AlignIO.read(os.path.join("..", "..", "data", "alignment.fasta"), "fasta")

# calculate distance - https://biopython.org/wiki/Phylo
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignments)

# create a translation dictionary for human understandable labels
trans_dict = dict(
    (alignment.id, " ".join(alignment.description.split()[1:3])) for alignment in alignments
)

# create dataframe from distance matrix for easier plotting
df = pd.DataFrame(
    dm.matrix,
    index=[trans_dict[name] for name in dm.names],
    columns=[trans_dict[name] for name in dm.names]
)

# create a heatmap of the distance dataframe
plt.figure()
sns.heatmap(
    df * 100,
    fmt='3.2f',
    annot=True,
    linewidths=0.5,
    cmap=sns.light_palette("navy"),
    cbar=False,
    square=True
)
plt.title("Percent difference")
plt.tight_layout()
plt.show()