from tools import *

seqs = read_fa("test.fa")
for s in seqs:
    print(s.seq.comple)