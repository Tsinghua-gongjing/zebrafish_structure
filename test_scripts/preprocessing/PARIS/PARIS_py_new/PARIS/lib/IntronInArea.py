
""" 
Usage: Get all introns overlapped with a given interval

Example:
    intronObj = intronInAreaClass(os.environ.get('HOME')+"/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed")
    intronObj.searchIntron('chrY', 25590513-10, 25592116+10, '-')
    [[25590514, 25592115]]

Coordination System:
    Input: [1-based-start, 1-based-end]
    Output: [1-based-intron-start, 1-based-intron-end]

"""

import sys, os, datetime

def now():
    "Return Current Time"
    return str(datetime.datetime.now())

def parseIntronInterval(genomeCoorFileName):
    """ Parse All Introns From genomeCoorFile

    """
    def parseIntron(exonString):
        exons = [ (int(exon.split('-')[0]), int(exon.split('-')[1])) for exon in  exonString.split(',') ]
        exons.sort(key=lambda x: x[0])
        intronList = []
        idx = 0
        while idx < len(exons) - 1:
            intron_start = int(exons[idx][1])+1
            intron_end = int(exons[idx+1][0])-1
            if intron_start > intron_end:
                print >>sys.stderr, "Error: Not A Valid Intron: %s" % (exonString, )
                return []
            intronList.append( [intron_start, intron_end] )
            idx += 1
        return intronList
    intron = {}
    IN = open(genomeCoorFileName)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        (Chr, Start, End, Strand, geneID, transID, geneType, Exon) = arr[:8]
        if Chr not in intron: 
            intron[Chr] = {}
            intron[Chr]['+'] = []
            intron[Chr]['-'] = []
        intronList = parseIntron(Exon)
        if intronList:
            intron[Chr][Strand] += intronList
        line = IN.readline()
    return intron


def binize(intron_Container, bw=100000):
    """ binize introns to speed up searching
        
    """
    def statisticPsuedoChrSize(intron_Container):
        chr_size = {}
        for Chr in intron_Container:
            last_1000_intron = intron_Container[Chr]['+'][-500:]+intron_Container[Chr]['-'][-500:]
            max_coor = 0
            for intron in last_1000_intron:
                if max_coor < max(intron):
                    max_coor = max(intron)
            chr_size[Chr] = max_coor + 1000
        return chr_size
    print >>sys.stderr,"Binize intron to speed up searching.\n\t%s" % (now(), )
    chr_size = statisticPsuedoChrSize(intron_Container)
    Bin = {}
    for Chr in chr_size.keys():
        #print >>sys.stderr, "Now Process Chromosome: %s.\n\t%s" % (Chr, now())
        Bin[Chr] = {}; Bin[Chr]['+'] = {}; Bin[Chr]['-'] = {}
        count = int( chr_size[Chr]/bw ) + 1
        for idx in range(count):
            Bin[Chr]['+'][idx] = []; Bin[Chr]['-'][idx] = []
        for Strand in ('+', '-'):
            intron_idx = 0
            for intron_idx in range(len(intron_Container[Chr][Strand])):
                Start = int( int(intron_Container[Chr][Strand][intron_idx][0]) / bw )
                End = int( int(intron_Container[Chr][Strand][intron_idx][1]) / bw )
                for bin_idx in range(Start, End+1):
                    Bin[Chr][Strand][bin_idx].append( intron_idx )
    return Bin

def searchIntronInterval(intron_Container, Bin, Chr, Start, End, Strand, bw=100000, verbose=True):
    if Chr not in Bin:
        if verbose: print >>sys.stderr,"Chromosome %s not in genomeCoor File" % (Chr, )
        raise KeyError
    idxBin = int(Start/bw)
    overlapIntrons = []
    while idxBin <= int(End/bw):
        if idxBin > Bin[Chr][Strand].keys()[-1]:
            break
        for intron_idx in Bin[Chr][Strand][idxBin]:
            intron = intron_Container[Chr][Strand][intron_idx]
            if intron[0] <= End and Start <= intron[1]:
                overlapIntrons.append( intron )
        idxBin += 1
    return overlapIntrons


class intronInAreaClass(object):
    """ Parse Introns in a region with a genomeCoor File
    Example:
        intronObj = intronInAreaClass(os.environ.get('HOME')+"/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed")
        intronObj.searchIntron('chrY', 25590513-10, 25592116+10, '-')
    """
    def __init__(self, genomeCoorFileName):
        self.genomeCoorFileName = genomeCoorFileName
        self.intron_Container = parseIntronInterval(genomeCoorFileName=genomeCoorFileName)
        self.Bin = binize(self.intron_Container, bw=100000)
    def searchIntron(self, Chr, Start, End, Strand, verbose=True):
        return searchIntronInterval(self.intron_Container, self.Bin, Chr=Chr, Start=Start, End=End, Strand=Strand, verbose=verbose)



