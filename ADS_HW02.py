# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 06:30:51 2015

@author: Temporary Account
"""
# Coursera Algorithms for DNA sequencing Homework #2
# William Wan

# import previously written Boyer-Moore preprocessing algorithm written
import bm_preproc as bm
import bisect
import sys

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

class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer
    
    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """

    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq

    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def queryIndex(p, t, index):
    k = index.k # k-mer size used to create the index
    offsets = []
    for i in index.query(p):
        if p[k:] == t[i+k:i+len(p)]:  # verify that rest of P matches
            offsets.append(i)
    return offsets

def queryIndex_approx(p, t, index, mm): # modification of queryIndex function to use pigeon hole principle and kmers in p
# checks before and after kmer index hits for up to mm mismatches
# returns correct indices, but does too many matches    
    occurrences = []
    hits = index.query(p) # initialize matches with 1st kmer in p
    uniques = hits # matches indexed to where p would start in t
    hits_raw = index.query(p) # non-unique matches
    k = index.k
    for i in range(1, len(p)-k): # look for all k-mers in p that also match in t
        indexTemp = index.query(p[i:]) # create temporary array of hits for each k-mer in p

        for j in indexTemp:
#            print(not (j-i) in uniques)
            if (j-i) not in uniques: # see if current hit referenced to start of p has been found already
                hits += [j]
                uniques += [j-i] # append unique match list

        hits_raw += indexTemp # list of all matches 

#    hits = list(hits)
    for i in hits: # naive matching for before and after exact matches
        mismatches = 0
        exact = t[i:i+k] # this is the k-mer from p that matched t 
        p_ind = p.find(exact) # start of exact in p

        for j in range(p_ind): # before kmer
            if p[j] != t[i-p_ind+j]:
                mismatches += 1
                if mismatches > mm:
                    break

        for j in range(p_ind+k, len(p)): # after kmer
            if p[j] !=  t[i+j-p_ind]:
                mismatches += 1
                if mismatches > mm:
                    break
        if mismatches <= mm:
            occurrences.append(i-p_ind) # adding p_ind gets us closer
    return occurrences, hits, hits_raw, uniques
#   sum(uniques) < sum(hits) since uniques is list of matches WRT start of P

# test cases
t = 'ATGCCTTGCA'
p = 'GCCT'
t2 = 'ATGCTAGCTACGTTACGT'
p2 = 'GCT'
p_bm = bm.BoyerMoore(p, alphabet='ACGT')
p2_bm = bm.BoyerMoore(p2, alphabet='ACGT')
test1 = boyer_moore(p, p_bm, t)
test2 = boyer_moore(p2, p2_bm, t2)
test2n = naive(p2, t2)
# output of  - naive(p, t) is ([2], 10, 7) - seems correct

t = 'ACTTGGAGATCTTTGAGGCTAGGTATTCGGGATCGAAGCTCATTTCGGGGATCGATTACGATATGGTGGGTATTCGGGA'
p = 'GGTATTCGGGA'
index3 = Index(t, 3)
print(queryIndex(p, t, index3))
print(queryIndex_approx(p, t, index3, 0))
# end test cases



# main body for quiz problems

genome = readGenome('chr1.GRCh38.excerpt.fasta')

# Question #1 - How many alignments does the naive exact matching algorithm try 
# when matching the string 
p1 = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
q1 = naive(p1, genome)

# 1st attempt ([56922], 984143, 799954) - 799954 is correct

# Question #2 - How many character comparisons does the naive exact matching 
# algorithm try when matching the string 
p2 = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
q2 = naive(p2, genome)

# 1st attempt ([56922], 984143, 799954) - 984143 is correct

# Question #3 - How many alignments does Boyer-Moore try 
p3 = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
p3_bm = bm.BoyerMoore(p3, alphabet='ACGT')
q3 = boyer_moore(p3, p3_bm, genome)

# 1st attempt ([56922], 165191, 127974) - 127974 is correct


# Question #4 - 
# Write a function that, given a length-24 pattern P and given an Index object built 
# on 8-mers, finds all approximate occurrences of P within T with up to 2 mismatches. 
# Insertions and deletions are not allowed.

# How many times does the string GGCGCGGTGGCTCACGCCTGTAAT 
#( derived from human Alu sequences) occur with up to 2 substitutions in the excerpt of 
# human chromosome 1? Hint: multiple index hits might direct you to the same match 
# multiple times; be careful not to count a match more than once.
genomeIndex = Index(genome, 8)
p4 = 'GGCGCGGTGGCTCACGCCTGTAAT'

offsetsExact = queryIndex(p4, genome, genomeIndex)
offsets2mm = queryIndex_approx(p4, genome, genomeIndex, 2) # 19 items based on Q6

# no more than 13 times based on k-mer indexing
# 5 exact matches using queryIndex function

# 1st attempt 10 hits is incorrect - 13 exact kmer hits returned - need to modify function to return hits with 2 mismatches
# 2nd attempt 19 hits based on function from Q6 (correct answer)


# Question #5 -
p5 = 'GGCGCGGTGGCTCACGCCTGTAAT'

q5 = genomeIndex.query(p5)
len(queryIndex_approx(p4, genome, genomeIndex, 0)[1])

# 1st attempt len(q5) ] 13 is incorrect
# 2nd attempt used pigeon-hole principle to assume that exact kmer matches mean no more than 2 matches before or after the kmer
len(queryIndex_approx(p4, genome, genomeIndex, 0)[1]) # gives 446 as total number of hits returned is also incorrect

# 117 is also wrong even though function returns correct number of occurrences


# Question #6

def approximate_ssmatch(p, t, index, mm):
    occurrences = []
    hits = []
#    p_bm = bm.BoyerMoore(p)
    for i in range(int(len(p)/index.k)): # go through hits for differents subsequences of p; need to confirm this is sufficient
        hits += (index.query(p[i:]))
        for hit in index.query(p[i:]): # go through hits for each matching subsequence
            indStart = hit-i # possible occurence is left of subsequence by i
            if not indStart in occurrences:
                mismatches = 0
                for j, mer in enumerate(p): # we want the mer in p as well as an index to look up t
                    if mer != t[indStart+j]:
                        mismatches += 1
    #                    print('mismatch')
                        if mismatches > mm:
                            break
                if mismatches <=mm:
                    occurrences.append(indStart)
    return occurrences, hits


p6 = 'GGCGCGGTGGCTCACGCCTGTAAT'
genomeSubSeqIndex = SubseqIndex(genome, 8, 3)

q6 = approximate_ssmatch(p6, genome, genomeSubSeqIndex, 2)
# 19 offsets matched p6 with up to 2 mismatches
q6Answer = len(set(q6[1])) # 79 hits is correct


