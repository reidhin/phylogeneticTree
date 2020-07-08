# libraries
import os
import re
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


# https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner
# https://plotly.com/~johnchase/22/visualizing-bioinformatics-data-with-plo/
# https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&SeqType_s=Nucleotide&Flags_ss=refseq&VirusLineage_ss=Middle%20East%20respiratory%20syndrome-related%20coronavirus%20(MERS-CoV),%20taxid:1335626
# https://www.cdc.gov/flu/about/viruses/types.htm
# https://piazza.com/class_profile/get_resource/ij0z6c2otfh5go/ilj276fchj86lg

# please cite us: https://www.ncbi.nlm.nih.gov/pubmed/27899678


def download_data(genome_accession_numbers: list) -> list:
    """
    Downloads the data from NCBI
    :param genome_accession_numbers: list of Entrez accession numbers
    :return: list of downloaded fasta sequences
    """
    Entrez.email = 'hans@orikami.nl'  # Always tell NCBI who you are
    search = " ".join(genome_accession_numbers)
    result = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
    seq_records = []
    for id in result['IdList']:
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
        seq_records.append(SeqIO.read(handle, "fasta"))
    return seq_records


def align_sequences(input_file: str, output_file: str = "virus_alignment.fasta") -> MultipleSeqAlignment:
    """
    Aligns the sequences using the muscle algorithm
    :param input_file: fasta-file with the input sequences
    :param output_file: save as aligned fasta-file
    :return: MultipleSeqAlignment with the alignment result
    """
    # run   muscle to align all sequences
    # can also be ran online:
    # https://www.ebi.ac.uk/Tools/services/web/toolresult.ebi?jobId=muscle-I20200329-210908-0063-87869209-p2m
    '''
    /nfs/public/ro/es/appbin/linux-x86_64/muscle-3.8.31/muscle -in muscle-I20200329-210908-0063-87869209-p2m.upfile -verbose -log muscle-I20200329-210908-0063-87869209-p2m.output -quiet -fasta -out muscle-I20200329-210908-0063-87869209-p2m.fasta -tree2 muscle-I20200329-210908-0063-87869209-p2m.dnd
    '''

    # specify where the muscle.exe is located
    muscle_exe = os.path.join('..', 'muscle.exe')

    # define the command line for muscle
    muscle_cline = MuscleCommandline(muscle_exe, input=os.path.join('..', 'data', input_file))

    # use 2 iterations; when sequences are far apart, the attempt to reach a more finer alignment leads to an error
    muscle_cline.maxiters = 2

    # report the final command line
    print(muscle_cline)

    # execute the command
    stdout, stderr = muscle_cline()

    # save for later faster processing or testing
    with open(os.path.join('..', 'data', output_file), "w") as alignment_file:
        alignment_file.write(stdout)

    # return the aligned sequences
    return AlignIO.read(StringIO(stdout), "fasta")


def plot_basics(recs, accession_numbers):
    fig, axs = plt.subplots(1, 2)
    axs[0].bar(range(len(recs)), [len(rec) for rec in recs])
    axs[0].set_ylabel('Length')
    axs[1].bar(range(len(recs)), [GC(rec.seq) for rec in recs])
    axs[1].set_ylabel('GC content')
    for ax in axs:
        ax.set_xticks(range(len(recs)))
        ax.set_xticklabels(
            [accession_numbers[re.match("(^\S*)(?=\.)", rec.id)[0]] for rec in recs],
            rotation=45,
            ha='right'
        )
    plt.tight_layout()
    return fig


def view_alignment(aln, accession_numbers):
    """matplotlib sequence alignment view"""

    # make sequence and id lists from the aln object
    mat = [
        [0 if nucleotide == '-' else 1 for nucleotide in rec.seq] for rec in aln
    ]
    cmap = ListedColormap(['w', 'r'])
    fig, ax = plt.subplots(1, 1)
    ax.matshow(mat, cmap=cmap)
    ax.set_aspect(ax.get_xlim()[1] / ax.get_ylim()[0] / 3)
    ax.set_yticks(range(len(aln)))
    ax.set_yticklabels(
        [accession_numbers[re.match("(^\S*)(?=\.)", rec.id)[0]] for rec in aln],
        fontsize=8
    )
    ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
    plt.tight_layout()
    return fig


def plot_phylo_tree(align: MultipleSeqAlignment, accession_numbers: dict):
    """
    Plots a phylogenetic tree
    :param align: MultipleSeqAlignment with the alignment result to be plotted
    :param accession_numbers: dict of accession numbers and their translation to human-understandable names
    :return: figure-handle of the plotted phylogenetic tree
    """
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
        terminal.name = accession_numbers[re.match("(^\S*)(?=\.)", terminal.name)[0]]

    print(Phylo.draw_ascii(tree))

    # plot the tree
    fig, ax = plt.subplots(1, 1)
    # draw the resulting tree
    Phylo.draw(tree, show_confidence=False, axes=ax, do_show=False)
    ax.set_xlim(right=0.8)
    return fig


def plot_phylo_tree_pars(align: MultipleSeqAlignment, accession_numbers: dict):
    """
    Plots a phylogenetic tree
    :param align: MultipleSeqAlignment with the alignment result to be plotted
    :param accession_numbers: dict of accession numbers and their translation to human-understandable names
    :return: figure-handle of the plotted phylogenetic tree
    """
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor = ParsimonyTreeConstructor(searcher)
    tree = constructor.build_tree(align)
    print(Phylo.draw_ascii(tree))

    # plot the tree
    fig, ax = plt.subplots(1, 1)
    # remove the names for the non-terminals for better visual appeal
    for non_terminal in tree.get_nonterminals():
        non_terminal.name = ''
    # draw the resulting tree
    Phylo.draw(tree, show_confidence=False, axes=ax, do_show=False)
    ax.set_xlim(right=0.8)
    return fig


if __name__ == '__main__':
    # find sequence accession numbers on https://www.ncbi.nlm.nih.gov/labs/virus
    accession_numbers = {
        "NC_026436": "H1N1 - swine flu",
        "NC_007373": "H3N2 - Hong Kong flu",
        "NC_007381": "H2N2 - Asian flu",
        "NC_007360": "H5N1 - bird flu",
        "NC_019843": "MERS",
        "NC_045512": "Corona",
        "MK062183": "SARS",
        "NC_006432": "Ebola",
        "NC_024781": "Marburg",
        "NC_001802": "HIV 1",
        "NC_001722": "HIV 2"
    }

    # download the sequences
    recs = download_data(list(accession_numbers.keys()))

    fig = plot_basics(recs)
    fig.savefig(os.path.join('..', 'figures', 'basic.png'))

    '''
    import seaborn as sns
    p1 = sns.scatterplot(
        [len(seq_record) for seq_record in seq_records],  # Horizontal axis
        [GC(seq_record.seq) for seq_record in seq_records],  # Vertical axis
        legend=False
    )

    for line in range(len(seq_records)):
        p1.text(
            len(seq_records[line]) + 0.01,
            GC(seq_records[line].seq),
            viruses.loc[re.match("(^\S*)(?=\.)", seq_records[line].id)[0], 'Virus'],
            horizontalalignment='left'
        )
    '''

    # write them in file for later upload
    #SeqIO.write(recs, os.path.join('..', 'data', "downloads.fasta"), "fasta")

    # align the sequences
    #align = align_sequences("downloads.fasta")
    align = AlignIO.read(os.path.join("..", "data", "virus_alignment.fasta"), "fasta")
    print(align)
    fig = view_alignment(align, accession_numbers)
    fig.savefig(os.path.join('..', 'figures', 'alignment.png'), bbox_inches='tight', dpi=300)

    # plot the resulting tree
    # plot_phylo_tree_pars(align, accession_numbers)
    fig = plot_phylo_tree(align, accession_numbers)
    fig.savefig(os.path.join('..', 'figures', 'tree.png'))

    plt.show()

    print('finished')
