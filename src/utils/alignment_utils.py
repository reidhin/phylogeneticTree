"""
Alignment utils
"""
import os
import pandas as pd
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align.Applications import MuscleCommandline
from Bio.Align import MultipleSeqAlignment
from io import StringIO
from Bio import AlignIO


def get_distance_dataframe(alignments, trans_dict=None):
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
    return df


def align_sequences(input_file: str, output_file: str = "alignment.fasta") -> MultipleSeqAlignment:
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
    muscle_exe = os.path.join('..', 'muscle3.8.31_i86linux64')

    # define the command line for muscle
    muscle_cline = MuscleCommandline(muscle_exe, input=input_file)

    # use 2 iterations; when sequences are far apart, the attempt to reach a more finer alignment leads to an error
    muscle_cline.maxiters = 2

    # report the final command line
    print(muscle_cline)

    # execute the command
    stdout, stderr = muscle_cline()

    # save for later faster processing or testing
    with open(output_file, "w") as alignment_file:
        alignment_file.write(stdout)

    # return the aligned sequences
    return AlignIO.read(StringIO(stdout), "fasta")
