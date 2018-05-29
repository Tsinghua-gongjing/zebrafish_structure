#-*- coding:utf-8 -*-


from Environment import *

def genPairClusterFromSamFile( samFile, read, readmap, task_name, removeRedundancy=False, minOverhang=2, useLongestPairing=False, preserveMultimap=False, intronFlanking=3, verbose=True ):
    """
    Example:
        read = {}
        readmap = {}
        task_name = 'try'
        samFile = "/Share/home/zhangqf8/lipan/paris/debug_test/rawData/Aligned.out.sort.sam"
        genPairClusterFromSamFile( samFile, read, readmap, task_name, removeRedundancy=False)
    """
    duplexGroupBed = GLOBAL.outDir+"/tmp."+task_name+".duplexGroup.bed"
    sortedSamFile = GLOBAL.outDir+"/tmp."+task_name+".sorted.sam"
    encodedSamFile = GLOBAL.outDir+"/tmp."+task_name+".encode.sam"
    # encodeSAMFile
    if removeRedundancy:
        cmd_1 = "grep \"^@\" %s > %s"
        if sortSupportParallel():
            #import multiprocessing
            #cores = multiprocessing.cpu_count()
            cmd_2 = "grep -v \"^@\" %s | sort --parallel=8 -k3,3 -k4,4n -k6,6 -k10,10 >> %s"
        else:
            cmd_2 = "grep -v \"^@\" %s | sort -k3,3 -k4,4n -k6,6 -k10,10 >> %s"
        runCMD( cmd_1 % (samFile, sortedSamFile) )
        runCMD( cmd_2 % (samFile, sortedSamFile) )
        uniqSam(sortedSamFile, encodedSamFile, read, readmap, preserveMultimap=preserveMultimap, task_name=task_name)
    else:
        encodeSam(samFile, encodedSamFile, read, readmap, preserveMultimap=preserveMultimap, task_name=task_name)
    samFile = encodedSamFile
    # define a Alert Object
    SAMLines = countFileLines(samFile, skipHead=['@'])
    Alerter = ProcessAlert(total=SAMLines, stepsToAlert=100000, OUT=GLOBAL.LOG)
    # parse DG from sam
    lineCount=0; validCount = 0
    SAM = open(samFile)
    DG = open(duplexGroupBed, 'w')
    print >>GLOBAL.LOG, 'read sam file %s...\n\t%s' % (samFile, now())
    line = SAM.readline()
    while line:
        if line.startswith('#'):
            pass
        elif line.startswith('@'):
            pass
        else:
            lineCount += 1
            Alerter.alert()
            duplexStemLine, duplexIntervalLine = genPairClusterFromSamLine( line, read, minOverhang=minOverhang, useLongestPairing=useLongestPairing, preserveMultimap=preserveMultimap, intronFlanking=intronFlanking, verbose=verbose)
            if duplexStemLine:
                print >>DG, duplexStemLine
                validCount += 1
        line = SAM.readline()
    Alerter.finish()
    SAM.close()
    print >>GLOBAL.LOG, "%d lines read from sam file: %s." % (lineCount, samFile)
    print >>GLOBAL.LOG, "among which %d lines generate supports for base pairing.\n\tTime: %s" % (validCount, now())
    # return 
    flush()
    return lineCount, validCount


def uniqSam(sortedSamFile, encodedSamFile, read, readmap, task_name, preserveMultimap=False):
    readMapFile = GLOBAL.outDir+"/tmp."+task_name+".readmap.txt"
    # define a Alert Object
    SAMLines = countFileLines(sortedSamFile, skipHead=['@'])
    Alerter = ProcessAlert(total=SAMLines, stepsToAlert=100000, OUT=GLOBAL.LOG)
    # open files
    SAM = open(sortedSamFile);
    OUT = open(encodedSamFile, 'w')
    MAP = open(readMapFile, 'w')
    showTime('read sam file: %s...\n\t%s' % (sortedSamFile, now()))
    # for recording the information of last line
    chrName = ""; pos = 0; cigar = "";
    lineCount = 0; uniqCount = 0
    line = SAM.readline()
    while line:
        if line.startswith('#'):
            pass
        elif line.startswith('@'):
            print >>OUT, line.strip()
        else:
            lineCount += 1
            Alerter.alert()
            data = line.strip().split("\t")
            # replace SRR3404926.24301833 with SRR3404926.24301833_1 1 is the multi-map index
            read_id = improved_readID(line)
            if not preserveMultimap and bool(int(data[1]) & 2304):
                line = SAM.readline()
                continue
            if data[3] != pos or data[2] != chrName or data[5] != cigar:
                uniqCount += 1
                pos = data[3]; chrName = data[2]; cigar = data[5];
                GLOBAL.readUniqCount += 1
                read[GLOBAL.readUniqCount] = READ()
                read[GLOBAL.readUniqCount].count = 1
                read[GLOBAL.readUniqCount].collapsedFrom = GLOBAL.readUniqCount
                read[GLOBAL.readUniqCount].collapsedTo = GLOBAL.readUniqCount
                read[GLOBAL.readUniqCount].name = read_id
                readmap[ read_id ] = GLOBAL.readUniqCount
                # write to file
                print >>MAP, read_id+"\t"+str(GLOBAL.readUniqCount)
                data[0] = str(GLOBAL.readUniqCount)
                print >>OUT, '\t'.join(data)
            else:
                print >>MAP, read_id+'\t'+str(GLOBAL.readUniqCount)
        line = SAM.readline()
    Alerter.finish()
    SAM.close(); MAP.close(); OUT.close()
    print >>GLOBAL.LOG, "%d lines read from sam file %s." % (lineCount, sortedSamFile)
    print >>GLOBAL.LOG, "among which %d lines are unique.\n\tTime: %s" % (uniqCount, now())
    flush()


def encodeSam(samFile, encodedSamFile, read, readmap, task_name, preserveMultimap=False):
    readMapFile = GLOBAL.outDir+"/tmp."+task_name+".readmap.txt"
    # define a Alert Object
    SAMLines = countFileLines(samFile, skipHead=['@'])
    Alerter = ProcessAlert(total=SAMLines, stepsToAlert=100000, OUT=GLOBAL.LOG)
    # open files
    SAM = open(samFile)
    OUT = open(encodedSamFile, 'w')
    MAP = open(readMapFile, 'w')
    showTime('read sam file: %s...\n\t%s' % (samFile, now()))
    lineCount = 0
    line = SAM.readline()
    while line:
        if line.startswith('#'):
            pass
        elif line.startswith('@'):
            print >>OUT, line.strip()
        else:
            lineCount += 1
            Alerter.alert()
            data = line.strip().split("\t")
            # replace SRR3404926.24301833 with SRR3404926.24301833_1 1 is the multi-map index
            read_id = improved_readID(line)
            if not preserveMultimap and bool(int(data[1]) & 2304):
                line = SAM.readline()
                continue
            GLOBAL.readUniqCount += 1
            read[GLOBAL.readUniqCount] = READ()
            read[GLOBAL.readUniqCount].count = 1
            read[GLOBAL.readUniqCount].collapsedFrom = GLOBAL.readUniqCount
            read[GLOBAL.readUniqCount].collapsedTo = GLOBAL.readUniqCount
            read[GLOBAL.readUniqCount].name = read_id
            # multi-map make it cann't work
            readmap[ read_id ] = GLOBAL.readUniqCount
            # write to file
            print >>MAP, read_id+"\t"+str(GLOBAL.readUniqCount)
            data[0] = str(GLOBAL.readUniqCount)
            print >>OUT, '\t'.join(data)
        line = SAM.readline()
    Alerter.finish()
    SAM.close(); MAP.close(); OUT.close()
    print >>GLOBAL.LOG, "%d lines read from sam file %s.\n\tTime: %s" % (lineCount, samFile, now())
    flush()

#genPairClusterFromSamFile( samFile, read, readmap, task_name, removeRedundancy=False, sorted=False )

def genPairClusterFromSamLine(line, read, minOverhang=2, intronFlanking=3, useLongestPairing=False, preserveMultimap=False, verbose=False):
    """ Parse a Line from sam file
    Example:
        samLine = "1 256 chr21   8209800 0   25M189522N30M2S *   0   0   AATACATGCCGACGGGCGCTGACCCCAGTCGGTCCTGAGAGATGGGCGAGCGCCGAG   IHIIIIIGGIHHHHHIIIIIHHHIIIIIDGHIHHIIHIGHIIHHIIIIIIHHIIIII   NH:i:10 HI:i:8  AS:i:47"
        read = {1:READ()}
        print genPairClusterFromSamLine(samLine, read)
    Return:
        chr21   8209800 8209824 1       1       +       chr21   8399347 8399376 1       1       +
        chr21   8209800 8399376 1       1       +
    """
    data = line.strip().split()
    if not preserveMultimap and bool(int(data[1]) & 2304):
        return None, None
    if 'N' not in data[5] and 'D' not in data[5]:
        if verbose: print >>GLOBAL.ERR, "Skip read that does not contain cleavage!\n\t%s" % (line, )
        return None, None
    strand = '-' if bool(int(data[1]) & 16) else '+'
    alignment, pair1s, pair1e, pair2s, pair2e = getSamPair( chrName=data[2], strand=strand, pos=int(data[3]), cigar=data[5], minOverhang=minOverhang, useLongestPairing=useLongestPairing, verbose=verbose );
    if not alignment:
        if verbose: print >>GLOBAL.ERR, "No valid duplex: \n\t%s" % (line, )
        return None, None
    if pair1s < pair2e and pair2s < pair1e:
        if verbose: print >>GLOBAL.ERR, 'inapproprieate duplex: \n\t%s' % (line, )
        return None, None
    # left arm
    read[ int(data[0]) ].left = Arm()
    read[ int(data[0]) ].left.chr = data[2]
    read[ int(data[0]) ].left.strand = strand
    read[ int(data[0]) ].left.start = pair1s
    read[ int(data[0]) ].left.end = pair1e
    # right arm
    read[ int(data[0]) ].right = Arm()
    read[ int(data[0]) ].right.chr = data[2]
    read[ int(data[0]) ].right.strand = strand
    read[ int(data[0]) ].right.start = pair2s
    read[ int(data[0]) ].right.end = pair2e
    # organize
    stemBed = '\t'.join([ data[2], `pair1s`, `pair1e`, data[0], '1', strand ])
    stemBed += '\t'+'\t'.join([ data[2], `pair2s`, `pair2e`, data[0], '1', strand ])
    intervalBed = '\t'.join( [data[2], `pair1s`, `pair2e`, data[0], '1', strand] )
    #print >>GLOBAL.ERR, stemBed + "\t" + alignment
    return stemBed, intervalBed


def getSamPair( chrName, strand, pos, cigar, intronFlanking=3, minOverhang=False, useLongestPairing=False, verbose=True ):
    """ Get Both Arm of the Reads
    Example:
        getSamPair( chrName='chr21', strand='+', pos=8209800, cigar='25M189522N30M2S', intronFlanking=3, minOverhang=False, useLongestPairing=False )
        # Check
        Genome.fetch('chr21', 8399347-1, 8399376, '+')
    Return:
        1-based-start-frag1, 1-based-end_frag1, 1-based-start-frag2, 1-based-end-frag2
    """
    match, matchSize = parseCigar( cigar )
    try:
        frag1, frag2, frag1Pos, frag2Pos = getBothFragment(chrName, strand, pos, match, matchSize, intronFlanking=intronFlanking, verbose=verbose)
    except KeyError:
        if verbose: print >>GLOBAL.ERR, "Chromosome %s not find in Genome, Skip it" % (chrName, )
        return None, None, None, None, None
    if not frag1 or not frag2:
        if verbose: print >>GLOBAL.ERR, "Skip read that aligns to an intron!\n\t%s\t%s\t%d\t%s" % (chrName, strand, pos, cigar)
        return None, None, None, None, None
    elif minOverhang and ( len(frag1) < minOverhang or len(frag2) < minOverhang ):
        if verbose: print >>GLOBAL.ERR, "Skip read that does not contain enough overhang in either end!"
        return None, None, None, None, None
    if useLongestPairing:
        frag2 = frag2[::-1]
        maxScore, alignment, intv1s, intv1e, intv2s, intv2e = localAlignment( frag1, frag2 )
        if not maxScore:
            if verbose: print >>GLOBAL.ERR, "Skip read that cannot be paired! %s %s %s %s" % (chrName, strand, pos, cigar)
            return None, None, None, None, None
        #print sys.stderr, "relative intervals: ", 
        if strand == '+':
            intv1s += frag1Pos - 1
            intv1e += frag1Pos
            tmpPos = len(frag2) - intv2s + frag2Pos
            intv2s = len(frag2) - intv2e + frag2Pos
            intv2e = tmpPos + 1
        else:
            tmpPos = len(frag1) - intv1s + frag1Pos
            intv1s = len(frag1) - intv1e + frag1Pos
            intv1e = tmpPos + 1
            intv2s += frag2Pos - 1
            intv2e += frag2Pos
        return alignment, intv1s, intv1e, intv2s, intv2e
    else:
        return 1, frag1Pos, frag1Pos+len(frag1), frag2Pos, frag2Pos+len(frag2)

# getSamPair( chrName='chr21', strand='+', pos=8209800, cigar='25M189522N30M2S', intronFlanking=3, minOverhang=10, useLongestPairing=True )


def getBothFragment(chrName, strand, pos, match, matchSize, intronFlanking=3, verbose=True):
    """ Get Both Seq From Genome
    Example:
        Genome = SeqFunc.seqClass('/Share/home/zhangqf8/lipan/paris/debug_test/annotation/human-virus.fa')
        match, matchSize = parseCigar( cigar='25M50722N30M2S' )
        getBothFragment(chrName='chr21', strand='+', pos=8392835, match=match, matchSize=matchSize)
        # check
        Genome.fetch('chr21', 8392835-1, 8392835+25-1, '+')
        Genome.fetch('chr21', 8392835+25+50722-1, 8392835+25+50722-1+30, '+')
    Return:
        frag1, frag2, 1-based-start-frag1, 1-based-start-frag2
    """
    frag1 = frag2 = ""
    maxNonMatch = maxND(chrName, strand, pos, match, matchSize, intronFlanking=intronFlanking, verbose=verbose)
    if not maxNonMatch: return frag1, frag2, 0, 0
    # infer genome position
    pos1 = pos
    for idx in range(maxNonMatch+1):
        if match[idx] in ('S', 'H', 'I', 'P'):
            pass
        elif match[idx] in ('M', '=', 'X'):
            frag1 += GLOBAL.Genome.fetch(chrName, pos-1, pos+matchSize[idx]-1, strand)
            pos += matchSize[idx]
        elif match[idx] == 'D':
            frag1 += GLOBAL.Genome.fetch(chrName, pos-1, pos+matchSize[idx]-1, strand) 
            pos += matchSize[idx]
        elif match[idx] == 'N':
            # the max N we find in most cases，special: 20M2000N20M3000N20M
            pos += matchSize[idx]
    pos2 = pos
    for idx in range(maxNonMatch+1, len(matchSize)):
        if match[idx] in ('S', 'H', 'I', 'P'):
            pass
        elif match[idx] in ('M', '=', 'X'):
            frag2 += GLOBAL.Genome.fetch(chrName, pos-1, pos+matchSize[idx]-1, strand)
            pos += matchSize[idx]
        elif match[idx] == 'D':
            frag2 += GLOBAL.Genome.fetch(chrName, pos-1, pos+matchSize[idx]-1, strand)
            pos += matchSize[idx]
        elif match[idx] == 'N':
            pos += matchSize[idx]
    return frag1, frag2, pos1, pos2


def maxND(chrName, strand, pos, match, matchSize, intronFlanking=3, verbose=True):
    """ Get the max N or D index from cigar
    Example:
        Intron = intronInAreaClass(os.environ.get('HOME')+"/lipan/DYNAMIC/GTF/hg38.genomeCoor.bed")
        # not in Intron
        match, matchSize = parseCigar( cigar='25M50722N30M2S' )
        maxND(chrName='chr21', strand='+', pos=8392835, match=match, matchSize=matchSize, intronFlanking=3)
        # in Intron
        match, matchSize = parseCigar( cigar='10M386N10M' )
        maxND(chrName='chr1', strand='+', pos=12218, match=match, matchSize=matchSize, intronFlanking=3)
    Return:
        >= 1 means a valid index
        == 0 means not found
    """
    totalLen = sum(matchSize)
    overlapped_introns = []
    if GLOBAL.Intron:
        try:
            #pass
            overlapped_introns = GLOBAL.Intron.searchIntron(chrName, pos-intronFlanking, pos+totalLen+intronFlanking, strand, verbose=verbose)
        except KeyError:
            "without intron information of this chromosome"
            if verbose: print >>GLOBAL.ERR, "Chromosome %s not find in Genome, Skip it" % (chrName, )
            overlapped_introns = []
    maxMatch = 0; tmpMax = 0;
    oldPos = pos
    for idx in range(len(matchSize)):
        if match[idx] in ('S', 'H', 'I', 'P'):
            pass
        elif match[idx] in ('M', '=', 'X'):
            pos += matchSize[idx]
        elif match[idx] in ('D', 'N'):
            oldPos = pos
            pos += matchSize[idx]
            isIntron = 0
            for interval in overlapped_introns:
                if abs( interval[0] - oldPos ) < intronFlanking and abs( interval[1] - pos ) < intronFlanking:
                    isIntron = 1
                    break
            if not isIntron and matchSize[idx] > tmpMax:
                tmpMax = matchSize[idx]
                maxMatch = idx
    return maxMatch




