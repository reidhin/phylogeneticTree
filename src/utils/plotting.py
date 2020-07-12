"""
Functions for plotting
"""
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import pandas as pd
from Bio.SeqUtils import GC
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment


def plot_basics(recs, labels):
    fig, axs = plt.subplots(1, 2)
    axs[0].bar(range(len(recs)), [len(rec) for rec in recs])
    axs[0].set_ylabel('Length')
    axs[1].bar(range(len(recs)), [GC(rec.seq) for rec in recs])
    axs[1].set_ylabel('GC content')
    for ax in axs:
        ax.set_xticks(range(len(recs)))
        ax.set_xticklabels(
            labels,
            rotation=45,
            ha='right'
        )
    plt.tight_layout()
    return fig


def view_alignment(aln, labels):
    #Todo: letters that are not the same should have another color
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
    ax.set_yticklabels(labels, fontsize=8)
    ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
    plt.tight_layout()
    return fig


def plot_phylo_tree(align: MultipleSeqAlignment, accession_numbers: dict, title=""):
    """
    Plots a phylogenetic tree
    :param align: MultipleSeqAlignment with the alignment result to be plotted
    :param accession_numbers: dict of accession numbers and their translation to human-understandable names
    :param title: string of title
    :return: figure-handle of the plotted phylogenetic tree
    """
    # calculate distance - https://biopython.org/wiki/Phylo
    calculator = DistanceCalculator('trans')
    dm = calculator.get_distance(align)

    # construct a tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    # order terminals
    tree.ladderize()

    # remove the names for the non-terminals for better visual appeal
    for non_terminal in tree.get_nonterminals():
        non_terminal.name = ''

    # change accession numbers into human more understandable names
    for terminal in tree.get_terminals():
        terminal.name = accession_numbers[terminal.name]

    print(Phylo.draw_ascii(tree))

    # plot the tree
    fig, ax = plt.subplots(1, 1)
    # draw the resulting tree
    Phylo.draw(tree, show_confidence=False, axes=ax, do_show=False)
    ax.set_xlim(right=1.3*ax.get_xlim()[1])
    ax.set_title(title)
    return fig


def print_gene_record(gene_record, line_length=int(100)):
    from src.utils.record_utils import print_record
    # print one gene record
    fig = plt.figure(figsize=(8, 8.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.text(0.05, 0.95, print_record(gene_record, line_length),
            horizontalalignment='left',
            verticalalignment='top',
            family='monospace',
            fontsize=8,
            transform=ax.transAxes)
    ax.set_axis_off()
    return fig


def plot_alignment_heatmap(alignments, trans_dict=None, title="Percent difference"):
    # calculate distance - https://biopython.org/wiki/Phylo
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignments)

    if trans_dict is None:
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
    plt.title(title)
    plt.tight_layout()
    return plt.gcf()

