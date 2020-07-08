# libraries
import matplotlib.pyplot as plt
import os
import seaborn as sns

from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

from src.utils.plotting import plot_phylo_tree, print_gene_record, plot_alignment_heatmap
from src.utils.record_utils import create_records
from src.utils.alignment_utils import get_distance_dataframe, align_sequences

# check dit
# https://www.ncbi.nlm.nih.gov/gene/?term=(%22Ape%22%5BOrganism%5D)+AND+HBB%5BGene%5D

# https://www.ensembl.org/Homo_sapiens/Gene/Compara_Tree?db=core;g=ENSG00000244734;r=11:5225464-5229395


if __name__ == '__main__':
    organism_list = [
        "Homo sapiens",
        "Pan troglodytes",
        "Gorilla gorilla",
        "Pongo abelii",
        "Pan paniscus"
#        "Hylobates moloch"
    ]

    gene_name = 'HBB'

    '''
    # get gene ids
    gene_ids = [get_gene_id(organism, gene_name) for organism in organism_list]

    # get gene records
    gene_records = [get_gene_record(gene_id) for gene_id in gene_ids]

    # write them in file for later upload
    SeqIO.write(gene_records, os.path.join('..', 'data', "trog_downloads.gb"), "genbank")
    '''

    gene_records = list(SeqIO.parse(os.path.join('..', 'data', "trog_downloads.gb"), "genbank"))

    # translate dict
    trans_dict = dict(
        zip(
            [gr.id for gr in gene_records], organism_list
        )
    )

    # print one gene record
    fig = print_gene_record(gene_records[0])
    fig.savefig(os.path.join("..", "figures", "human_gene_DNA.png"))

    # get amino-acid sequence
    amino_records = create_records(gene_records, 'amino')
    #fig = plot_basics(amino_records, organism_list)

    # save as fasta-file for alignment
    SeqIO.write(amino_records, os.path.join('..', 'data', "trog_before_alignment.fasta"), "fasta")

    # align the sequences
    align_amino = align_sequences("trog_before_alignment.fasta", output_file="trog_after_alignment.fasta")
    print(align_amino)

    aligner = Align.PairwiseAligner()
    alignments = aligner.align(amino_records[0].seq, amino_records[1].seq)
    for alignment in sorted(alignments):
        print("Score = %.1f:" % alignment.score)
        print(alignment)

    # get cds records
    cds_records = create_records(gene_records, 'cds')
    #fig = plot_basics(cds_records, organism_list)

    # save as fasta-file for alignment
    SeqIO.write(cds_records, os.path.join('..', 'data', "trog_before_alignment.fasta"), "fasta")

    # align the sequences
    align_cds = align_sequences("trog_before_alignment.fasta", output_file="trog_after_alignment.fasta")
    print(align_cds)

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignments = aligner.align(cds_records[0].seq, cds_records[1].seq)
    for alignment in sorted(alignments):
        print("Score = %.1f:" % alignment.score)
        print(alignment)

    # get introns sequence
    full_records = create_records(gene_records, 'cds_and_introns')
    #fig = plot_basics(full_records, organism_list)

    # save as fasta-file for alignment
    SeqIO.write(full_records, os.path.join('..', 'data', "trog_before_alignment.fasta"), "fasta")

    # align the sequences
    align_full = align_sequences("trog_before_alignment.fasta", output_file="trog_after_alignment.fasta")
    print(align_full)

    for al, name in zip([align_amino, align_cds, align_full], ["amino-acids", "exons", "exons and introns"]):
        fig = plot_alignment_heatmap(al, f'Percent difference based on {name}')
        fig.savefig(os.path.join('..', 'figures', f'distance_{name}_{gene_name}.png'))

    # plot the resulting tree
    fig = plot_phylo_tree(align_cds, trans_dict, title="Phylogenetic tree based on the CDS of the HBB gene")
    fig.savefig(os.path.join('..', 'figures', f'tree_{gene_name}.png'))

    df_cds = get_distance_dataframe(align_cds, trans_dict)
    df_full = get_distance_dataframe(align_full, trans_dict)
    df_cds = df_cds.fillna(0) + df_cds.transpose().fillna(0)
    df_full = df_full.fillna(0) + df_full.transpose().fillna(0)
    df_intron = (df_full * align_full.get_alignment_length() - df_cds * align_cds.get_alignment_length()) / \
                (align_full.get_alignment_length() - align_cds.get_alignment_length())
    sns.heatmap(
        df_intron*100,
        fmt='3.2f',
        annot=True,
        linewidths=0.5,
        cmap=sns.light_palette("navy"),
        cbar=False,
        square=True
    )
    plt.title('Percent difference')
    plt.tight_layout()
    plt.show()

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # save as fasta-file for alignment
    out_records = list()
    for i in range(2):
        cds = cds_records[i]
        cds.id = f"cds_{i}"
        full = full_records[i]
        full.id = f"full_{i}"
        out_records = out_records + [cds, full]

    SeqIO.write(
        out_records,
        os.path.join('..', 'data', "trog_before_alignment.fasta"),
        "fasta"
    )

    # align the sequences
    align_print = align_sequences("trog_before_alignment.fasta", output_file="trog_after_alignment.fasta")
    print(align_print)

    # todo: check if nucleotide is part of the cds by position
    d = dict(zip([a.id for a in align_print], [i for i in range(len(align_print))]))
    comb_0 = "".join([f.lower() if c == '-' else f for f, c in zip(align_print[d['full_0']], align_print[d['cds_0']])])
    comb_1 = "".join([f.lower() if c == '-' else f for f, c in zip(align_print[d['full_1']], align_print[d['cds_1']])])
    match = "".join(["|" if c1 == c2 else " " for c1, c2 in zip(comb_0, comb_1)])

    line_length = 100
    line = int(0)
    out = ""
    while line * line_length < align_print.get_alignment_length():
        r = range(line * line_length, (line + 1) * line_length)
        # print first sequence
        out = out + comb_0[r.start:r.stop] + "\n"
        out = out + match[r.start:r.stop] + "\n"
        out = out + comb_1[r.start:r.stop] + "\n"
        out = out + "\n"
        # print()
        line = line + 1

    # print one gene record
    fig = plt.figure(figsize=(8, 8.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.text(0.05, 0.95, out,
            horizontalalignment='left',
            verticalalignment='top',
            family='monospace',
            fontsize=8,
            transform=ax.transAxes)
    ax.set_axis_off()
    fig.savefig(os.path.join("..", "figures", "human_to_chimp.png"))
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    print('finished!')
