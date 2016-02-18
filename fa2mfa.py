#!/usr/bin/env python
"""
Convert multiple fasta files to a single mfa file.

USAGE:
fa2mfa.py --fadir file.gbk --prefix ecoli_k12

NOTE:


AUTHOR:
Andreas Sjödin
andreas.sjodin@gmail.com

Version 0.1

"""

import os
from argparse import ArgumentParser
from datetime import datetime

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid


def get_filepaths(directory, filetypes):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """

    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            if filename.endswith(filetypes):
                # Join the two strings in order to form the full filepath.
                filepath = os.path.join(root, filename)
                file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.


def get_fastaseq(path):
    """
    This function will recursively check for content in files
     and only include unique sequence records.
    """
    seq_dict = []
    seguid_dict = []

    files = get_filepaths(path, allowedTypes)
    for file in files:
        print("")
        print("Found", len(files), "records in", file)
        for seq_record in SeqIO.parse(file, "fasta"):
            if seguid(seq_record.seq) not in seguid_dict:
                seguid_dict.append(seguid(seq_record.seq))
                print(seq_record.id, "is added to the list")
                seq_dict.append(seq_record)
            else:
                print(seq_record.id, "is already in list")
    return seq_dict


if __name__ == '__main__':
    description = "Combine multiple fasta files to a single multi-fasta file without duplicates"

    parser = ArgumentParser(description=description)
    parser.add_argument("--fadir", required=True, help="FASTA directory")
    parser.add_argument("--prefix", required=True, help="prefix")
    args = parser.parse_args()

    t0 = datetime.now()
    allowedTypes = ('fna', 'fa', 'fasta', 'fsa')

    path = os.path.abspath(args.fadir)
    print("#############################################################")
    print("## Starting conversion \n")
    print("Looking for files in folder: ", path)

    fasta_seq = get_fastaseq(path)
    outputfile = args.prefix + ".fasta"
    SeqIO.write(fasta_seq, outputfile, "fasta")

    dt = datetime.now() - t0
    print(dt)

