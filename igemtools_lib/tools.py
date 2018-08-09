from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from itertools import chain

import random

BASE_LOOKUP = {
    1: "A",
    2: "T",
    3: "C",
    4: "G",
    "A": 1,
    "T": 2,
    "C": 3,
    "G": 4
}

def read_fa(file):
    seqs = []
    with open(file,  "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record.alphabet = generic_dna
            seqs.append(record)
    return seqs

def is_compatible(seqs, working_bases):
    """
    Determine compatibility

    :param working_bases:
    :return: bool
    """
    all_seqs = []
    all_seqs.append([s.seq for s in working_bases])
    all_seqs.append([s.seq for s in seqs])
    all_seqs = chain(all_seqs)

    # Check for palindromes
    for s in working_bases:
        if is_palindrome(s):
            print("is palindrome!")
            return False

    # Check for compatibility5, -4, -1, -.1
    if not is_specific(all_seqs):
        return False

    return True

def is_palindrome(s, l=6):
    """
    Determine if sequence contains a palindrome
    Split task up into parts and parallelize
    Treat as binary tree
    :param s:
    :return:
    """

    nr = [
        BASE_LOOKUP[b]
        for b in s.seq
    ]
    nr_str = [
        str(n)
        for n in nr
    ]
    nr_str_rev = nr_str[::-1]
    words = np.array([
        int(''.join(nr_str[r:r+l]))
        for r in range(len(nr_str)-l+1)
    ])
    words_rev = np.array([
        int(''.join(nr_str_rev[r:r+l]))
        for r in range(len(nr_str_rev)-l+1)
    ])

    mat1 = np.repeat(words, len(words), axis=0)
    mat1 = mat1.reshape((len(words), len(words)))
    mat2 = np.repeat(words_rev, len(words_rev), axis=0)
    mat2 = mat2.reshape((len(words_rev), len(words_rev)))
    mat2 = np.transpose(mat2)

    subtr = np.subtract(mat1, mat2)
    truth = np.where(subtr == 0, True, False)
    if np.sum(truth) > 0:
        return False
    else:
        return True



def is_specific(all_seqs):
    for s in all_seqs:
        for ss in all_seqs:
            if not s == ss:
                aln = pairwise2.align.globalms(s[0], ss[0], 3, -1, -2, -1, score_only=True)
                if aln > 0:
                    return False
    return True

def generate_bases(gc=0.5, l=60, r=3):
    working_bases = []
    for i in range(l):

        # Prevent repeats
        b =  BASE_LOOKUP[random.randint(1, 4)]
        if len(working_bases) >= r:
            while working_bases[-1*r].count(b) >= r:
                print("found same, {}, {}, {}".format(working_bases[-1*r], b, r))
                b = BASE_LOOKUP[random.randint(1, 4)]
        working_bases.append(b)

    if 0.40 < cg_content(working_bases) < 0.60:
        s = SeqRecord(Seq(''.join(working_bases)))
        s.alphabet = generic_dna
        return s
    else:
        return generate_bases(gc, l)



def cg_content(b):
    cg = float((b.count("C") + b.count("G")) / len(b))
    return cg

def write_fa(file, seqs):
    with open(file, "w") as output_handle:
        for s in seqs:
            SeqIO.write(s, output_handle, "fasta")
