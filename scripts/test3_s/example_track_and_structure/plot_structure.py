import sys
import re
from visual import *

foldfile = open(sys.argv[1], 'r')
sequence = "AGCTGGGTTTCCCGATT"
structure = "....(((...)))...."
values = [0.8]*4+[0.2]*3+[0.9]*3+[0.2]*3+['NULL']*4

title = ''
title_m = re.compile('>.*>(\S+)') 
lines = foldfile.readlines()
title = title_m.match(lines[0].strip()).groups()[0] if title_m.match(lines[0].strip()) else None
highlight_start = int(lines[0].split('|')[-1].split('.')[0])
highlight_end = int(lines[0].split('|')[-1].split('.')[1])
sequence = lines[1].strip()
structure = lines[2].strip()
assert len(sequence) == len(structure)


values = []
try:
    shapefile = open(sys.argv[2], 'r')
    for l in shapefile:
        l = l.strip()
        values.append(l)
    if len(values) != len(sequence):
        values.extend(['NULL'] * (len(sequence) - len(values)))
except IndexError:
    print "no values provided!"
    values = ['NULL'] * len(sequence)
# Plot_RNAStructure_Shape(seq=sequence, ss=structure, shape_list=values, mode='fill', title=title,VARNAProg='/Users/gongjing/Downloads/VARNAv3-93-src.jar')
Plot_RNAStructure_highlight(seq=sequence, ss=structure, hg_base_list=[], mode='fill', correctT=True, 
	highlight_region=[(highlight_start,highlight_end)], title=sys.argv[1].split('/')[-1], wait=True, VARNAProg="/Users/gongjing/Downloads/VARNAv3-93-src.jar",
	outputfile=sys.argv[1].replace('.txt', '.plot.svg'))

