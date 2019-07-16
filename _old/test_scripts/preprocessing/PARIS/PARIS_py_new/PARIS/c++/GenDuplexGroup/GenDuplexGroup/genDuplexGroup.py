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


def genDuplexGroup(duplexGroupBedFileSortedFile, duplex_group_name, minOverlap=5, multipleDG=False):
    # process BED
    #duplexGroupBedFileSortedFile = "/Users/lee/Desktop/paris-simple/c++/data/dupleGroup.txt"
    #duplex_group_name = "/Users/lee/Desktop/paris-simple/c++/data/find_dupleGroup_python.txt";
    duplexCount = 0
    lineCount = 0; firstPossible = 0
    BED = open(duplexGroupBedFileSortedFile)
    line = BED.readline()
    duplexGroup = {}
    while line:
        if line.startswith('#'):
            pass
        else:
            lineCount += 1
            ( chr1, start1, end1, id1, score1, strand1, chr2, start2, end2, id2, score2, strand2 ) = line.strip().split()
            lastDGoverlapped = 0
            nonOverlapped = 1
            for idx in range(firstPossible, duplexCount):
                overlapped = checkOverlap( duplexGroup[idx], chr1, strand1, int(start1), int(end1), chr2, strand2, int(start2), int(end2) )
                if overlapped >= minOverlap:
                    nonOverlapped = 0
                    lastDGoverlapped = 1
                    addRead2DuplexGroup ( duplexGroup[idx], id1, int(start1), int(end1), int(start2), int(end2) )
                    #if read[int(id1)].clique == -1:
                    #    read[int(id1)].clique = `idx`
                    #else:
                    #    read[int(id1)].clique += ";"+`idx`
                    if not multipleDG: break
                elif overlapped == -1:
                    ## the start of this read is beyond the firstPossible duplex, so no need to check firstPossible in the future
                    if not lastDGoverlapped:
                        firstPossible = idx + 1
                else:
                    lastDGoverlapped = 1
            if nonOverlapped:
                newDuplexGroup ( duplexGroup, duplexCount, id1, chr1, strand1, int(start1), int(end1), chr2, strand2, int(start2), int(end2) )
                #read[int(id1)].clique = `duplexCount`
                duplexCount += 1
        line = BED.readline()
    BED.close()
    OUT = open(duplex_group_name, 'w')
    for dgID in duplexGroup:
        print >>OUT, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
            (duplexGroup[dgID].chr1, duplexGroup[dgID].start1, duplexGroup[dgID].end1, duplexGroup[dgID].strand1, \
            duplexGroup[dgID].chr2, duplexGroup[dgID].start2, duplexGroup[dgID].end2, duplexGroup[dgID].strand2, \
            duplexGroup[dgID].reads, duplexGroup[dgID].support )
    OUT.close()


def checkOverlap(dg_item, chr1, strand1, start1, end1, chr2, strand2, start2, end2):
    overlap = 0
    if dg_item.end1 < start1 or dg_item.chr2 != chr2 or dg_item.strand2 != strand2 or \
        dg_item.chr1 != chr1 or dg_item.strand1 != strand1:
        overlap = -1
    elif dg_item.start1 <= end1 and dg_item.start2 <= end2 and start2 <= dg_item.end2:
        # dg_item have overlaps with read
        #overlap1 = ( end1 if dg_item.end1 > end1 else dg_item.end1 ) - ( dg_item.start1 if dg_item.start1 > start1 else start1 )   
        #overlap2 = ( end2 if dg_item.end2 > end2 else dg_item.end2 ) - ( dg_item.start2 if dg_item.start2 > start2 else start2 )
        overlap1 = min(end1, dg_item.end1) - max(start1, dg_item.start1)
        overlap2 = min(end2, dg_item.end2) - max(start2, dg_item.start2)
        overlap = min(overlap1, overlap2)
    return overlap

def addRead2DuplexGroup( dg_item, id1, start1, end1, start2, end2 ):
    dg_item.support += 1
    dg_item.reads += ';' + id1
    dg_item.start1 = max(start1, dg_item.start1) 
    dg_item.end1 = min(end1, dg_item.end1) 
    dg_item.start2 = max(start2, dg_item.start2)
    dg_item.end2 = min(end2, dg_item.end2) 

def newDuplexGroup( duplexGroup, idx, read_id, chr1, strand1, start1, end1, chr2, strand2, start2, end2 ):
    duplexGroup[idx] = DuplexGroup()
    duplexGroup[idx].support = 1;
    duplexGroup[idx].reads = read_id
    duplexGroup[idx].chr1 = chr1
    duplexGroup[idx].strand1 = strand1
    duplexGroup[idx].start1 = start1
    duplexGroup[idx].end1 = end1
    duplexGroup[idx].chr2 = chr2
    duplexGroup[idx].strand2 = strand2
    duplexGroup[idx].start2 = start2
    duplexGroup[idx].end2 = end2

import sys

genDuplexGroup(sys.argv[1], sys.argv[2], minOverlap=5, multipleDG=False)












