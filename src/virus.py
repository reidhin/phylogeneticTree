# libraries
import os
import re
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio import Entrez
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO


# https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&SeqType_s=Nucleotide&Flags_ss=refseq&VirusLineage_ss=Middle%20East%20respiratory%20syndrome-related%20coronavirus%20(MERS-CoV),%20taxid:1335626
# https://www.cdc.gov/flu/about/viruses/types.htm

# please cite us: https://www.ncbi.nlm.nih.gov/pubmed/27899678


def download_data(genomeAccessionNumbers: list):
    Entrez.email = 'hans@orikami.nl'  # Always tell NCBI who you are
    search = " ".join(genomeAccessionNumbers)
    result = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
    seq_records = []
    for id in result['IdList']:
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
        seq_records.append(SeqIO.read(handle, "fasta"))
    return seq_records

'''
h = Entrez.einfo()
record = Entrez.read(h)
h = Entrez.esearch(
    db='nucleotide',
    term='ebola AND refseq[type]'
#    term="NC_002549.1"
)
result = Entrez.read(h)

with Entrez.efetch(
    db="nucleotide", rettype="fasta", retmode="text", id=result['IdList'][0]
) as handle:
    seq_record = SeqIO.read(handle, "fasta")


h.close()
'''

if __name__ == '__main__':
    print(os.getcwd())

    # find sequence accession numbers on https://www.ncbi.nlm.nih.gov/labs/virus

    accession_numbers = {
        'NC_026436': 'H1N1 - swine flu',
        "NC_007373": "H3N2 - Hong Kong flu",
        "NC_007381": "H2N2 - Asian flu",
        'NC_007360': 'H5N1 - bird flu',
    #    'NC_019843': 'MERS',
    #    'NC_045512': 'Corona',
    #    "MK062183": "SARS",
        "NC_006432": "Ebola",
        "NC_024781": "Marburg"
    }
    '''
    accession_numbers = {
        'NC_045512': 'Corona',
        'NC_044932': 'Norwalk',
        'NC_043585': 'Ilesha'
    }
    '''
    # download the sequences
    seqs = download_data(list(accession_numbers.keys()))

    # write them in file for later upload
    SeqIO.write(seqs, os.path.join('..', 'data', "downloads.fasta"), "fasta")

    # run   muscle to align all sequences
    # can also be ran online:
    # https://www.ebi.ac.uk/Tools/services/web/toolresult.ebi?jobId=muscle-I20200329-210908-0063-87869209-p2m
    muscle_exe = os.path.join('..', 'muscle.exe')
    muscle_cline = MuscleCommandline(muscle_exe, input=os.path.join('..', 'data', "downloads.fasta"))
    muscle_cline.maxiters = 2
    print(muscle_cline)
    stdout, stderr = muscle_cline()

    # save for faster processing or testing
    with open(os.path.join('..', 'data', "alignment.fasta"), "w") as alignment_file:
        alignment_file.write(stdout)

    # read the aligned sequences
    align = AlignIO.read(StringIO(stdout), "fasta")
    print(align)

    # calculate distance - https://biopython.org/wiki/Phylo
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(align)
    dm.names = [accession_numbers[re.match("(^\S*)(?=\.)", x)[0]] for x in dm.names]
    print(dm)

    # construct a tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    print(tree)
    print(Phylo.draw_ascii(tree))

    fig, ax = plt.subplots(1, 1)
    for non_terminal in tree.get_nonterminals():
        non_terminal.name = ''
    Phylo.draw(tree, show_confidence=False, axes=ax)
    ax.set_xlim(right=0.8)


    print('finished')

