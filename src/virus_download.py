# libraries
import os, re
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from Bio.SeqUtils import GC

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

print('finished!')