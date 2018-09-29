#-*- coding:utf-8 -*-

from Environment import *

# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 备份
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
"""

import copy
read_bak = copy.deepcopy(read)
readmap_bak = copy.deepcopy(readmap)

# 恢复
read = copy.deepcopy(read_bak)
readmap = copy.deepcopy(readmap_bak)



# 备份2
_uniqRead = GLOBAL.readUniqCount

# 恢复2
_readIDNames = read.keys()
_count = 0
for readID in _readIDNames:
    if readID > _uniqRead:
        _count += 1
        name = read[readID].name
        del readmap[name]
        del read[readID]

print 'delete %d items' % (_count, )
GLOBAL.readUniqCount = _uniqRead
"""

def genPairClusterFromJunctionFile( junctionFile, read, readmap, task_name, removeRedundancy=False, minOverhang=2, useLongestPairing=False, intronFlanking=3, verbose=True ):
    """
    Example:
        read = {}
        readmap = {}
        task_name = 'try'
        junctionFile = "/Share/home/zhangqf8/lipan/paris/debug_test/rawData/Chimeric.out.junction"
        genPairClusterFromJunctionFile( junctionFile, read, readmap, task_name, removeRedundancy=False)
    No multi-map in junction file
    """
    duplexGroupBed = GLOBAL.outDir+"/tmp."+task_name+".duplexGroup.bed"
    sortedJunctionFile = GLOBAL.outDir+"/tmp."+task_name+".sorted.junction"
    encodeJunctionFile = GLOBAL.outDir+"/tmp."+task_name+".encode.junction"
    # encodeJunctionFile
    if removeRedundancy:
        if sortSupportParallel():
            #import multiprocessing
            #cores = multiprocessing.cpu_count()
            cmd = "sort --parallel=8  -k1,1 -k4,4 -k3,3 -k6,6 -k11,11n -k13,13n -k12,12 -k14,14 %s >> %s"
        else:
            cmd = "sort  -k1,1 -k4,4 -k3,3 -k6,6 -k11,11n -k13,13n -k12,12 -k14,14 %s >> %s"
        runCMD( cmd % (junctionFile, sortedJunctionFile) )
        uniqJunction(sortedJunctionFile, encodeJunctionFile, read, readmap, task_name)
    else:
        encodeJunction(junctionFile, encodeJunctionFile, read, readmap, task_name)
    junctionFile = encodeJunctionFile
    # define a Alert Object
    JUNCLines = countFileLines(junctionFile, skipHead=['@'])
    Alerter = ProcessAlert(total=JUNCLines, stepsToAlert=2000, OUT=GLOBAL.LOG, careRecentTimes=5)
    # parse DG from junction
    lineCount = 0; validCount = 0
    JUNC = open(junctionFile); 
    DG = open(duplexGroupBed, 'a')
    print >>GLOBAL.LOG, 'read junction file %s...\n\t%s' % (junctionFile, now())
    line = JUNC.readline()
    while line:
        if line.startswith('#'):
            pass
        else:
            lineCount += 1
            Alerter.alert()
            duplexStemLine, duplexIntervalLine = genPairClusterFromOneJunction( line, read, minOverhang=minOverhang, useLongestPairing=useLongestPairing, intronFlanking=intronFlanking, verbose=verbose )
            if duplexStemLine:
                print >>DG, duplexStemLine
                validCount += 1
        line = JUNC.readline()
    Alerter.finish()
    JUNC.close()
    print >>GLOBAL.LOG, "%d lines read from junction file: %s." % (lineCount, junctionFile)
    print >>GLOBAL.LOG, "among which %d lines generate supports for base pairing.\n\tTime: %s" % (validCount, now())
    # return 
    return lineCount, validCount

def uniqJunction(sortedJunctionFile, encodeJunctionFile, read, readmap, task_name):
    readMapFile = GLOBAL.outDir+"/tmp."+task_name+".readmap.txt"
    # define a Alert Object
    JUNCLines = countFileLines(sortedJunctionFile, skipHead=['@'])
    Alerter = ProcessAlert(total=JUNCLines, stepsToAlert=100000, OUT=GLOBAL.LOG, careRecentTimes=3)
    # open files
    JUNC = open(sortedJunctionFile); 
    OUT = open(encodeJunctionFile, 'w')
    MAP = open(readMapFile, 'a')
    showTime('read junction file: %s...\n\t%s' % (sortedJunctionFile, now()))
    # for recording the information of last line
    chrName1 = ""; strand1 = ""; pos1 = 0; cigar1 = "";
    chrName2 = ""; strand2 = ""; pos2 = 0; cigar2 = "";
    lineCount = 0; uniqCount = 0
    line = JUNC.readline()
    while line:
        if line.startswith('#'):
            pass
        elif line.startswith('@'):
            print >>OUT, line.strip()
        else:
            lineCount += 1
            Alerter.alert()
            data = line.strip().split("\t")
            if int(data[10]) != pos1 or int(data[12]) != pos2 or \
                data[11] != cigar1 or data[13] != cigar2 or \
                data[0] != chrName1 or data[3] != chrName2 or \
                data[2] != strand1 or data[5] != strand2:
                    uniqCount += 1
                    # record
                    pos1 = int(data[10]); pos2 = int(data[12])
                    cigar1 = data[11]; cigar2 = data[13]
                    chrName1 = data[0]; chrName2 = data[3]
                    strand1 = data[2]; strand2 = data[5]
                    # add to read
                    GLOBAL.readUniqCount += 1
                    read[GLOBAL.readUniqCount] = READ()
                    read[GLOBAL.readUniqCount].count = 1
                    read[GLOBAL.readUniqCount].collapsedFrom = GLOBAL.readUniqCount
                    read[GLOBAL.readUniqCount].collapsedTo = GLOBAL.readUniqCount
                    read[GLOBAL.readUniqCount].name = data[9]
                    readmap[ data[9] ] = GLOBAL.readUniqCount
                    # write to file
                    print >>MAP, data[9]+"\t"+str(GLOBAL.readUniqCount)
                    data[9] = str(GLOBAL.readUniqCount)
                    print >>OUT, '\t'.join(data)
            else:
                print >>MAP, data[9]+'\t'+str(GLOBAL.readUniqCount)
        line = JUNC.readline()
    Alerter.finish()
    JUNC.close(); MAP.close(); OUT.close()
    print >>GLOBAL.LOG, "%d lines read from juntion file %s." % (lineCount, sortedJunctionFile)
    print >>GLOBAL.LOG, "among which %d lines are unique.\n\tTime: %s" % (uniqCount, now())


def encodeJunction(junctionFile, encodeJunctionFile, read, readmap, task_name):
    readMapFile = GLOBAL.outDir+"/tmp."+task_name+".readmap.txt"
    # define a Alert Object
    JUNCLines = countFileLines(junctionFile, skipHead=['@'])
    Alerter = ProcessAlert(total=JUNCLines, stepsToAlert=100000, OUT=GLOBAL.LOG, careRecentTimes=3)
    # open files
    JUNC = open(junctionFile); 
    OUT = open(encodeJunctionFile, 'w')
    MAP = open(readMapFile, 'a')
    showTime('read junction file: %s...\n\t%s' % (junctionFile, now()))
    lineCount = 0
    line = JUNC.readline()
    while line:
        if line.startswith('#'):
            pass
        elif line.startswith('@'):
            print >>OUT, line.strip()
        else:
            lineCount += 1
            Alerter.alert()
            data = line.strip().split("\t")
            GLOBAL.readUniqCount += 1
            read[GLOBAL.readUniqCount] = READ()
            read[GLOBAL.readUniqCount].count = 1
            read[GLOBAL.readUniqCount].collapsedFrom = GLOBAL.readUniqCount
            read[GLOBAL.readUniqCount].collapsedTo = GLOBAL.readUniqCount
            read[GLOBAL.readUniqCount].name = data[9]
            readmap[ data[9] ] = GLOBAL.readUniqCount
            # write to file
            print >>MAP, data[9]+"\t"+str(GLOBAL.readUniqCount)
            data[9] = str(GLOBAL.readUniqCount)
            print >>OUT, '\t'.join(data)
        line = JUNC.readline()
    Alerter.finish()
    JUNC.close(); MAP.close(); OUT.close()
    print >>GLOBAL.LOG, "%d lines read from junction file %s.\n\tTime: %s" % (lineCount, junctionFile, now())



def genPairClusterFromOneJunction(line, read, minOverhang=2, intronFlanking=3, useLongestPairing=False, verbose=False):
    """ Parse a Line from junction file
    Example:
        junctionLine = "chr11  19992294    +   chr11   20017298    +   0   0   0   1   19992271    23M36S  20017299    23S36M"
        read = {1:READ()}
        print genPairClusterFromOneJunction(junctionLine, read)[0]
    Return:
        chr11   19992271        19992294        1       1       +       chr11   20017299        20017335        1       1       +
        chr11   19992271        20017335        1       1       +
    """
    data = line.strip().split()
    cigar = ''; isChiastic = 0;
    if data[0] != data[3] or data[2] != data[5]:
        isChiastic = 2
        cigar = data[11] + "|" + data[13]
    else:
        if data[2] == '+':
            if int(data[1]) > int(data[4]):
                isChiastic = 1; 
                cigar = getNewCigar(isChiastic=1, strand='+', start1=int(data[10]), start2=int(data[12]), frag1Cigar=data[11], frag2Cigar=data[13])
            else:
                # can be writen in a line
                cigar = getNewCigar(isChiastic=0, strand='+', start1=int(data[10]), start2=int(data[12]), frag1Cigar=data[11], frag2Cigar=data[13])
        elif data[2] == '-':
            if int(data[1]) < int(data[4]):
                isChiastic = 1; 
                cigar = getNewCigar(isChiastic=1, strand='-', start1=int(data[10]), start2=int(data[12]), frag1Cigar=data[11], frag2Cigar=data[13])
            else:
                # can be writen in a line
                cigar = getNewCigar(isChiastic=0, strand='-', start1=int(data[10]), start2=int(data[12]), frag1Cigar=data[11], frag2Cigar=data[13])
        else:
            print >>GLOBAL.ERR, "Unexpected Junction Line: \n\t%s" % (line, )
            return None, None
    if not cigar:
        print >>GLOBAL.ERR, "Skip line of inapproprieate alignment: \n\t%s" % (line, )
        return None, None
    alignment, pair1s, pair1e, pair2s, pair2e = getJuncPair( chr1=data[0], strand1=data[2], donor=int(data[1]), pos1=int(data[10]), cigar1=data[11], 
                                                            chr2=data[3], strand2=data[5], acceptor=int(data[4]), pos2=int(data[12]), cigar2=data[13], 
                                                            intronFlanking=intronFlanking, minOverhang=minOverhang, useLongestPairing=useLongestPairing,
                                                            verbose=verbose )
    if not alignment:
        print >>GLOBAL.ERR, "No valid duplex: \n\t%s" % (line, )
        return None, None
    # modified 2017-7-7
    #if (data[0] == data[3]) and (data[2] == data[5]) and (pair1s < pair2e) and (pair2s < pair1e):
    if (pair1s < pair2e) and (pair2s < pair1e):
        print >>GLOBAL.ERR, "inapprorieate duplex: \n\t%s" % (line, )
        return None, None
    readID = int(data[9])
    read[ readID ].cigar = `isChiastic`+':'+cigar
    # left arm
    read[ readID ].left = Arm()
    read[ readID ].left.chr = data[0]
    read[ readID ].left.strand = data[2]
    # right arm
    read[ readID ].right = Arm()
    read[ readID ].right.chr = data[3]
    read[ readID ].right.strand = data[5]
    if data[0] == data[3] and data[2] == data[5]:
        if pair1s <= pair2s:
            read[ readID ].left.start = pair1s
            read[ readID ].left.end = pair1e
            read[ readID ].right.start = pair2s
            read[ readID ].right.end = pair2e
        else:
            read[ readID ].right.start = pair1s
            read[ readID ].right.end = pair1e
            read[ readID ].left.start = pair2s
            read[ readID ].left.end = pair2e
    else:
        read[ readID ].left.start = pair1s
        read[ readID ].left.end = pair1e
        read[ readID ].right.start = pair2s
        read[ readID ].right.end = pair2e
    # insure from small to large
    if read[ readID ].left.chr > read[ readID ].right.chr:
        # exchange two elemets
        read[ readID ].left, read[ readID ].right = read[ readID ].right, read[ readID ].left
    elif read[ readID ].left.chr == read[ readID ].right.chr and read[ readID ].left.strand > read[ readID ].right.strand:
        # exchange two elemets
        read[ readID ].left, read[ readID ].right = read[ readID ].right, read[ readID ].left
    # organize
    stemBed = '\t'.join([ read[ readID ].left.chr, `read[ readID ].left.start`, `read[ readID ].left.end`, data[9], '1', read[ readID ].left.strand ])
    stemBed += '\t'+'\t'.join([ read[ readID ].right.chr, `read[ readID ].right.start`, `read[ readID ].right.end`, data[9], '1', read[ readID ].right.strand ])
    intervalBed = stemBed
    # organize
    """
    stemBed = '\t'.join([ data[0], `pair1s`, `pair1e`, data[9], '1', data[2] ])
    stemBed += '\t'+'\t'.join([ data[3], `pair2s`, `pair2e`, data[9], '1', data[5] ])
    intervalBed = stemBed
    if data[0] == data[3] and data[2] == data[5]:
        if pair1s < pair2s: 
            intervalBed = "\t".join([ data[0], `pair1s`, `pair2e`, data[9], '1', data[2] ])
        else:
            intervalBed = "\t".join([ data[0], `pair2s`, `pair1e`, data[9], '1', data[2] ])
    """
    #print >>GLOBAL.ERR, stemBed + "\t" + alignment
    return stemBed, intervalBed


def getJuncPair( chr1, strand1, donor, pos1, cigar1, 
                chr2, strand2, acceptor, pos2, cigar2,
                intronFlanking=3, minOverhang=False, useLongestPairing=False, verbose=True ):
    """ Get Both Arm of the Reads
    Example:
        getJuncPair( chr1='chrX', strand1='-', donor=143884355, pos1=143884356, cigar1='26S20M1S', 
                        chr2='chr21', strand2='+', acceptor=8215620, pos2=8215621, cigar2='21S26M',
                        intronFlanking=3, minOverhang=False, useLongestPairing=True )
        #Check useLongestPairing=False
        GLOBAL.Genome.fetch('chrX', 143884356-1, 143884376, '-')
        GLOBAL.Genome.fetch('chr21', 8215621-1, 8215647, '+')
        #Check useLongestPairing=True
        GLOBAL.Genome.fetch('chrX', 143884357-1, 143884370, '-')
        GLOBAL.Genome.fetch('chr21', 8215623-1, 8215636, '+')
    Return:
        1-based-start-frag1, 1-based-end_frag1, 1-based-start-frag2, 1-based-end-frag2
    """
    isIntron = 0
    if chr1 == chr2 and strand1 == strand2:
        if strand1 == '+' and donor < acceptor: 
            isIntron = checkJuncIntron( chr1, '+', donor, acceptor, verbose=verbose )
        elif strand1 == '-' and donor > acceptor: 
            isIntron = checkJuncIntron( chr1, '-', acceptor, donor, verbose=verbose )
    if isIntron:
        print >>GLOBAL.ERR, "Skip read that align to a intron!"
        return None, None, None, None, None
    try:
        frag1 = getOneFragment(chr1, strand1, pos1, cigar1)
        frag2 = getOneFragment(chr2, strand2, pos2, cigar2)
    except KeyError:
        return None, None, None, None, None
    #print frag1, frag2
    if minOverhang and ( len(frag1) < minOverhang or len(frag2) < minOverhang ):
        print >>GLOBAL.ERR, "Skip read that does not contain enough overhang in either end!"
        return None, None, None, None, None
    if useLongestPairing:
        frag2 = frag2[::-1]
        maxScore, alignment, intv1s, intv1e, intv2s, intv2e = localAlignment( frag1, frag2 )
        if not maxScore:
            print >>GLOBAL.ERR, "Skip read that cannot be paired!"
            return None, None, None, None, None
        #print sys.stderr, "relative intervals: ", 
        if strand1 == '+':
            intv1s += pos1 - 1
            intv1e += pos1
        else:
            tmpPos = len(frag1) - intv1s + pos1
            intv1s = len(frag1) - intv1e + pos1
            intv1e = tmpPos + 1
        if strand2 == '+':
            tmpPos = len(frag2) - intv2s + pos2
            intv2s = len(frag2) - intv2e + pos2
            intv2e = tmpPos + 1
        else:
            intv2s += pos2 - 1
            intv2e += pos2
        return alignment, intv1s, intv1e, intv2s, intv2e
    else:
        return 1, pos1, pos1+len(frag1), pos2, pos2+len(frag2)


def checkJuncIntron(chrName, strand, pos1, pos2, intronFlanking=3, verbose=True):
    """ Check if the chrName-strand-pos1-pos2 is the intron region
    Example:
        checkJuncIntron(chrName='chr1', strand='+', pos1=12057, pos2=12179)
        checkJuncIntron(chrName='chr1', strand='+', pos1=12057-10, pos2=12179)
    """
    overlapped_introns = []
    if GLOBAL.Intron:
        try:
            overlapped_introns = GLOBAL.Intron.searchIntron(chrName, pos1-intronFlanking, pos2+intronFlanking, strand)
        except KeyError:
            if verbose: print >>GLOBAL.ERR, "Chromosome %s not find in Genome, Skip it" % (chrName, )
            overlapped_introns = []
    for interval in overlapped_introns:
        if abs( interval[0]-pos1 ) < intronFlanking and abs( interval[1] - pos2 ) < intronFlanking:
            return True
    return False


def getOneFragment(chrName, strand, pos, cigar):
    """ infer the sequence of Match part 
    Example:
        getOneFragment(chrName='chr17', strand='-', pos=31886010, cigar='23S33M1S')
        getOneFragment(chrName='chr11', strand='+', pos=20017299, cigar='23S36M')

        GLOBAL.Genome.fetch('chr17', 31886010-1, 31886010+33, '-')
        GLOBAL.Genome.fetch('chr11', 20017299-1, 20017299+36, '+')
    """
    largetsMatch = parseCigar(cigar, getLargestM=True)
    frag = GLOBAL.Genome.fetch(chrName, pos-1, pos+largetsMatch-1, strand)
    return frag


def getNewCigar(isChiastic, strand, start1, start2, frag1Cigar, frag2Cigar):
    """ infer a new cigar from junction reads
    Examples:
        getNewCigar(isChiastic=1, strand='+', start1=3375624, start2=3375110, frag1Cigar='36M35S', frag2Cigar='36S35M')
        # 35M479N36M
        getNewCigar(isChiastic=1, strand='-', start1=55509122, start2=55509750, frag1Cigar='20S27M', frag2Cigar='20M27S')
        # 27M601N20M

        getNewCigar(isChiastic=0, strand='+', start1=19992271, start2=20017299, frag1Cigar='23M36S', frag2Cigar='23S36M')
        # 23M25005N36M
        getNewCigar(isChiastic=0, strand='-', start1=31886010, start2=31880567, frag1Cigar='23S33M1S', frag2Cigar='23M34S')
        # 23M5420N33M1S

    Example of unapprorieat:
        chr6    32978286    +   chr6    32978220    +   0   0   1   4304    32978252    2S34M36S    32978221    36S25M11S
        getNewCigar(isChiastic=1, strand='+', start1=32978252, start2=32978221, frag1Cigar='2S34M36S', frag2Cigar='36S25M11S')
    """
    match1, matchSize1 = parseCigar(frag1Cigar)
    match2, matchSize2 = parseCigar(frag2Cigar)
    # infer cigar
    cigar = ''
    newReadFrag1 = ''; newReadFrag2 = ''
    if not isChiastic:
        if strand == '+':
            lenN = start2 - start1
            for idx in range(len(match1)-1):
                cigar += str(matchSize1[idx])+str(match1[idx])
                if match1[idx] in ('M', '=', 'X', 'D', 'N'):
                    lenN -= matchSize1[idx]
            if lenN < 0:
                cigar = 0
            else:
                cigar += str(lenN)+'N'
                for idx in range(1, len(match2)):
                    cigar += str(matchSize2[idx])+match2[idx]
        elif strand == '-':
            lenN = start1 - start2
            for idx in range(len(match2)-1):
                cigar += str(matchSize2[idx])+match2[idx]
                if match2[idx] in ('M', '=', 'X', 'D', 'N'):
                    lenN -= matchSize2[idx]
            if lenN < 0:
                cigar = 0
            else:
                cigar += `lenN`+'N'
                for idx in range(1, len(match1)):
                    cigar += `matchSize1[idx]` + match1[idx]
    else:
        if strand == '+':
            lenN = start1 - start2
            for idx in range(1, len(match2)):
                # Why
                if match2[idx] == 'S':
                    cigar += `matchSize2[idx]` + 'M'
                else:
                    cigar += `matchSize2[idx]` + match2[idx]
                if match2[idx] in ('M', '=', 'X', 'D', 'N', 'S'):
                    lenN -= matchSize2[idx]
            # Why
            if match1[0] == 'S':
                lenN -= matchSize1[0]
            if lenN < 0:
                cigar = 0
            else:
                cigar += `lenN` + 'N'
                for idx in range(len(match1) - 1):
                    if match1[idx] == 'S':
                        cigar += `matchSize1[idx]` + 'M'
                    else:
                        cigar += `matchSize1[idx]` + match1[idx]
        if strand == '-':
            lenN = start2 - start1
            for idx in range(1, len(match1)):
                # Why
                if match1[idx] == 'S':
                    cigar += `matchSize1[idx]` + 'M'
                else:
                    cigar += `matchSize1[idx]` + match1[idx]
                if match1[idx] in ('M', '=', 'X', 'D', 'N', 'S'):
                    lenN -= matchSize1[idx]
            # Why
            if match2[0] == 'S':
                lenN -= matchSize2[0]
            if lenN < 0:
                cigar = 0
            else:
                cigar += `lenN`+'N'
                for idx in range(len(match2)-1):
                    if match2[idx] == 'S':
                        cigar += `matchSize2[idx]` + 'M'
                    else:
                        cigar += `matchSize2[idx]` + match2[idx]
    return cigar


#cigar = getNewCigar(isChiastic=1, strand='-', start1=26436106, start2=205616212, frag1Cigar="21S17M1S", frag2Cigar="21M18S")






