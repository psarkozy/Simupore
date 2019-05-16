import math
import random

def generate_k_mers(k):
    k_mers = ['']
    for i in range(k):
        newkmers = []
        for pmer in k_mers:
            for base in 'ACGT':
                newkmers.append(pmer + base)
        k_mers = newkmers
    return k_mers

def load_kmer_mean_spreads(kmer_mean_spread_file):
    kmer_mean_spread = {}
    for line in open(kmer_mean_spread_file).readlines()[1:]:
        lsp = line.split('\t', 6)
        kmer_mean_spread[lsp[0]] = (float(lsp[1]), float(lsp[2]))
    return kmer_mean_spread


def revcmp(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def load_kmer_lookahead_times(lookaheadfile):

    kmer_lookahead_dt ={}
    for line in open(lookaheadfile).readlines():
        line = line.split()
        kmer_lookahead_dt[line[0]] = float(line[1])
    return kmer_lookahead_dt
    #for kmer in generate_k_mers(4):
    #    kmer_lookahead_dt[kmer] = 1.1
    #return

def generate_readid(): #f2c7ef17-b081-444c-ad41-895b2c5729f7
    out = []
    for i in [8,4,4,4,12]:
        part = ''
        for j in range(i):
            part+= random.choice('012345679abcdef')
        out.append(part)
    return '-'.join(out)