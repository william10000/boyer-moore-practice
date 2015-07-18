# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 06:30:51 2015

@author: Temporary Account
"""
# Coursera Algorithms for DNA sequencing Homework #1
# William Wan

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
    
def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def naive(p, t):
    occurrences = []
    
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def naive_with_rc(p, t):
    occurrences = []
    p_rc = reverseComplement(p)
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters in original string
                match = False
                break
        if match == False: # check reverse complement if forward seq doesn't match
            match = True
            for j in range(len(p)):  # loop over characters
                if t[i+j] != p_rc[j]:  # compare characters in rev. complement
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

genome = readGenome('lambda_virus.fa')

q1 = naive_with_rc('AGGT', genome)
q2 = naive_with_rc('TTAA', genome)

q3 = naive_with_rc('ACTAAGT', genome)
q4 = naive_with_rc('AGTCGA', genome)

q5 = naive_2mm('TTCAAGCC', genome)
q6 = naive_2mm('AGGAGGTT', genome)

#1 q1 - q6 are probably correct

################ question 7
def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def phred33ToQ(qual):
    return ord(qual) - 33

def createHist(qualities):
    # Create a histogram of quality scores
    hist = [0]*50
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist

def findGCByPos(reads):
    ''' Find the GC ratio at each position in the read '''

    # Keep track of the number of G/C bases and the total number of bases at each position
    gc = [0] * 100
    totals = [0] * 100

    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
    
    # Divide G/C counts by total counts to get the average at each position
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc

def totQualByPos(quals):
    totals = [0] * 100
    for reads in quals:
        for i in range(len(reads)):
            totals[i] += phred33ToQ(reads[i])
    return totals

seqs, quals = readFastq('ERR037900_1.first1000.fastq')

h = createHist(quals)
gc = findGCByPos(seqs)
posQuals = totQualByPos(quals)

# print(h)

# Plot the histograms

import matplotlib.pyplot as plt
#plt.plot(range(len(h)), h)
#plt.show()

plt.plot(range(len(gc)), gc)
plt.show()

plt.plot(range(len(posQuals)), posQuals)
plt.show()

q7 = gc
q72 = posQuals

print('Answers are: ')
print('Q1: ', len(q1))
print('Q2: ', len(q2))
print('Q3: ', q3[0])
print('Q4: ', q4[0])
print('Q5: ', len(q5))
print('Q6: ', q6[0])
print('Q7: ', gc.index(min(gc)))

