#!/usr/bin/env python3

import io
from Bio import SeqIO

FASTA = """
>mule
TAATACCCCGGATATATGTCCTCACATAGTTCGAGGTCGAGAAAAATGAC
TTCCCACCAAGTGGACTCAGCTCGAGTAAACGCCAACGATACGTCCATTA
GGTGTGTGCCgcaactagtcggacccgttgtgacggaaacaggtccccgc
caagtcacacgggcatgtcatggacTCTCGATCGTTCATCGCCTTCTTGG
GTACCGCAGCCGCAATTAAGCCGTGTCTTCTTCCCCCTTCAAACGGGAAT
CGTGTCGACTTCTTAGGAGCAGNNNNNNNNNNCTAACTCCAGAG
>donkey
TAATACCCCGGATATATGTCTTAACATAGTTCCAGGTCGAGAAGAATGAC
TTGCCACCAAGTGGACTCAGATTCAGTCAACGCGAACGATAAGTCCATTA
GGTGTGTACCgcaactagtgggacccgttgtgacggaaacaggtcaccgc
caagtcacacgtgcatgtcatgtacTCTCGATCGTTTATCGCCTTCTTGG
GTACCGCAGCCGAAATTAAGCCGTGTCTTCTTCCCACTTCAAACGGGAAT
CGTGTCGACTTTACAGGAACAGNNNNNNNNNNATAACGCCAGAG
"""

def baseasint(b):
    if 'A' == b:
        return 0
    if 'C' == b:
        return 1
    if 'G' == b:
        return 2
    if 'T' == b:
        return 3
    if 'N' == b:
        return 4

with io.StringIO(FASTA) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        print(record.id)
        encoded = list(map(baseasint, record.seq.upper()))

encoded


#
# Test bigwig file
#
import pyBigWig
bw = pyBigWig.open("../Data/DNASE/fold_coverage_wiggles/DNASE.H1-hESC.fc.signal.bigwig")
bw
bw.chroms()
bw.header()
# Mean value in range
bw.stats("chr1", 1000000, 1003000)
# Max value in range
bw.stats("chr1", 1000000, 1003000, type="max")
# Max value in chromosome
bw.stats("chr1", type="max")
# All values in range
bw.values("chr1", 1000000, 1003000)
bw.values("chr1", 1000000, 1000030)
bw.values("chr1", 0, 3)
