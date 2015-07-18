# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 06:30:51 2015

@author: Temporary Account
"""
# Coursera Algorithms for DNA sequencing Homework #2
# William Wan

# import previously written Boyer-Moore preprocessing algorithm written
import bm_preproc as bm

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
    
def naive(p, t):
    occurrences = []
    comp = 0 # initialize # of character comparisons
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            comp += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, comp, i+1 # returns indicies of matches, # character comparisons, # of alignments tried

# copied from class practical
def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching """
    i = 0
    occurrences = []
    comp = 0 # initialize # of character comparisons
    aligns = 0 # initializes # of alignments tried
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            comp += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
        aligns += 1
    return occurrences, comp, aligns


# test cases
t = 'ATGCCTTGCA'
p = 'GCCT'
p_bm = bm.BoyerMoore(p, alphabet='ACGT')
t1 = boyer_moore(p, p_bm, t)
# output of  - naive(p, t) is ([2], 10, 7) - seems correct



# main body for quiz problems

genome = readGenome('chr1.GRCh38.excerpt.fasta')

# Question #1 - How many alignments does the naive exact matching algorithm try 
# when matching the string 
p1 = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
q1 = naive(p1, genome)

# 1st attempt ([56922], 984143, 799954)

# Question #2 - How many character comparisons does the naive exact matching 
# algorithm try when matching the string 
p2 = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
q2 = naive(p2, genome)

# 1st attempt ([56922], 984143, 799954)

# Question #3 - How many alignments does Boyer-Moore try 
p3 = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
p3_bm = bm.BoyerMoore(p3, alphabet='ACGT')
q3 = boyer_moore(p3, p3_bm, genome)

# 1st attempt ([56922], 165191, 127974)









