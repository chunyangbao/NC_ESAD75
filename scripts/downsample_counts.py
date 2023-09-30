#!/usr/bin/env python3

import os
import sys
import random

def downsample(counts, prop):
    ds_counts = []
    for count in counts:
        d = sum([int(random.random() <= prop) for read in range(count)])
        ds_counts.append(str(d))
    return ds_counts

# Usage
if len(sys.argv) != 3:
    exit("Usage:\n\t" + sys.argv[0] + " colname_1,colname_2,colname_3 <input.tsv> output.tsv\n\nInput was:\n\tA tab-separated values (TSV) file")

# argv
c = sys.argv[1] # column name, e.g. 'participant_id,Mutect2_passed_combined_vcf,Mutect2_passed_combined_vcf_index'
p = sys.argv[2] 

c = c.split(',')
p = float(p)

random.seed(10)

# main
l = sys.stdin.readline()
l = l.strip('\n').split('\t')
h = dict((map(reversed, enumerate(l)))) # column name:indices
ci = [h[cx] for cx in c] # column name indices
sys.stdout.write('\t'.join(l) + '\n') # print the 1st line

for l in sys.stdin:
    l = l.strip('\n').split('\t')
    lc = [int(l[cix]) for cix in ci]
    dc = downsample(lc, p)
    for x,y in zip(ci, dc):
        l[x] = y
    sys.stdout.write('\t'.join(l) + '\n')

### The End ###
