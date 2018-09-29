#-*- coding:utf-8 -*-


class DuplexGroup(object):
    __slots__ = ['support', 'reads', 'chr1', 'strand1', 'start1', 'end1', 'chr2', 'strand2', 'start2', 'end2', 'collapsedTo', 'cov']
    def __init__(self):
        self.support = -1
        self.reads = -1
        self.chr1 = -1
        self.strand1 = -1
        self.start1 = -1
        self.end1 = -1
        self.chr2 = -1
        self.strand2 = -1
        self.start2 = -1
        self.end2 = -1
        self.collapsedTo = -1
        self.cov = -1



def checkOverlapDG(dg_item1, dg_item2, maxGap=10, maxTotal=30, checkReads=False):
    overlap = 0
    if dg_item1.end1 + maxGap < dg_item2.start1 or \
        dg_item1.chr2 != dg_item2.chr2 or dg_item1.strand2 != dg_item2.strand2 or \
        dg_item1.chr1 != dg_item2.chr1 or dg_item1.strand1 != dg_item2.strand1:
        overlap = -1
    else:
        gap2 = dg_item2.start2 - dg_item1.end2 if dg_item2.start2 > dg_item1.start2 else dg_item1.start2 - dg_item2.end2
        total1 = max(dg_item2.end1, dg_item1.end1) - dg_item1.start1
        total2 = max(dg_item1.end2, dg_item2.end2) - min(dg_item1.start2, dg_item2.start2)
        if total1 < maxTotal and total2 < maxTotal and gap2 < maxGap:
            overlap = 1
    return overlap


def mergeDuplexGroup(dg_item1, dg_item2):
    dg_item1.reads += ';'+dg_item2.reads
    dg_item1.support += dg_item2.support
    dg_item1.start1 = min(dg_item1.start1, dg_item2.start1)
    dg_item1.end1 = max(dg_item1.end1, dg_item2.end1)
    dg_item1.start2 = min(dg_item1.start2, dg_item2.start2)
    dg_item1.end2 = max(dg_item1.end2, dg_item2.end2)


def collapseDuplexGroup(file1, file2, maxGap=10, maxTotal=30):
    firstPossible = 0
    IN = open(file1)
    OUT = open(file2, 'w')
    duplexGroup = []
    line = IN.readline()
    # start to collapse
    iii = 0
    while line:
        (chr1, start1, end1, strand1, chr2, start2, end2, strand2, reads, support) = line.strip().split()
        dg = DuplexGroup();
        dg.chr1 = chr1;  dg.start1 = int(start1);  dg.end1 = int(end1);  dg.strand1 = strand1;  
        dg.chr2 = chr2;  dg.start2 = int(start2);  dg.end2 = int(end2);  dg.strand2 = strand2;
        dg.reads = reads; dg.support = int(support);
        duplexGroup.append( dg )
        iii += 1
        lastDGoverlapped = 0
        for idx in range(firstPossible, len(duplexGroup) - 1):
            overlapped = checkOverlapDG(duplexGroup[idx], dg, maxGap=maxGap, maxTotal=maxTotal)
            if overlapped == -1:
                if not lastDGoverlapped:
                    firstPossible = idx + 1
            elif overlapped > 0:
                lastDGoverlapped = 1
                mergeDuplexGroup(duplexGroup[idx], dg)
                #dg.collapsedTo = 22
                del duplexGroup[-1]
                break
            else:
                lastDGoverlapped = 1
        line = IN.readline()
    for dg in duplexGroup:
        print >>OUT, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (dg.chr1, dg.start1, dg.end1, dg.strand1, \
            dg.chr2, dg.start2, dg.end2, dg.strand2, dg.reads, dg.support)


import sys
collapseDuplexGroup(sys.argv[1], sys.argv[2], maxGap=10, maxTotal=30)

