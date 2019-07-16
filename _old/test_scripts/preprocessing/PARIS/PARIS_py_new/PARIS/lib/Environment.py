#!/Share/home/zhangqf/usr/anaconda/bin/python
#-*- coding:utf-8 -*-

import re, sys, os, getopt, time, datetime, commands, random
import seq as SeqFunc
from IntronInArea import *
import Memory

class GLOBAL():
    outDir = ""
    readUniqCount = 0
    Intron = None
    Genome = None
    LOG = sys.stdout
    ERR = sys.stderr
    @staticmethod
    def init_GLOBAL():
        GLOBAL.outDir = ""
        GLOBAL.readUniqCount = 0
        GLOBAL.Intron = None
        GLOBAL.Genome = None
        GLOBAL.LOG = sys.stdout
        GLOBAL.ERR = sys.stderr


"""
import seq as SeqFunc
from IntronInArea import *

GLOBAL.Genome = SeqFunc.seqClass(GenomeFile)
GLOBAL.Intron = intronInAreaClass(GenomeCoorFile)
"""

class READ(object):
    __slots__ = ['count', 'collapsedTo', 'collapsedFrom', 'name', 'left', 'right', 'cigar', 'clique', 'cluster', 'ngTag']
    def __init__(self):
        self.count = -1
        self.collapsedTo = -1
        self.collapsedFrom = -1
        self.name = -1   # readID: SRR3404926.37587503
        self.left = -1
        self.right = -1
        self.cigar = -1
        self.clique = -1 # duplex group
        self.cluster = -1
        self.ngTag = -1

class Arm(object):
    __slots__ = ['chr', 'strand', 'start', 'end']
    def __init__(self):
        self.chr = -1
        self.strand = -1
        self.start = -1
        self.end = -1

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

class ProcessAlert(object):
    """ A Class To Record The Prosedure of the Current Task
    Example:
        alertH = ProcessAlert(total=283029, stepsToAlert=10000, OUT=sys.stdout)
        for i in range(283029):
            alertH.alert()
        alertH.finish()
    Exception:
        If alert called times exceed the total, UnboundLocalError will be raised
    """
    def __init__(self, total, stepsToAlert, OUT, careRecentTimes=10):
        self.total = total
        self.careRecentTimes = careRecentTimes
        self.recent = []
        self.current = 0
        self.lastTime = datetime.datetime.now()
        self.startTime = datetime.datetime.now()
        self.stepsToAlert = stepsToAlert
        self.totalSteps = total / stepsToAlert
        self.OUT = OUT
    @staticmethod
    def microSeconds2Human(miSecond):
        """ Convert miscro Seconds into Human Read Format
            microSeconds2Human(2832382382)
            # '47.21 mins'
        """
        unit = ['useconds', 'mseconds', 'seconds', 'mins', 'hours', 'days', 'weeks', 'months', 'years']
        relation = [1000, 1000, 60, 60, 24, 7, 4, 12]
        idx = 0
        while idx < 8 and miSecond > relation[idx]:
            miSecond = 1.0 * miSecond / relation[idx]
            idx += 1
        return "%.2f %s" % (miSecond, unit[idx])
    @staticmethod
    def timeInterval(time_first, time_next):
        return (time_next - time_first).total_seconds() * 1000 * 1000
    def alert(self, other=''):
        if self.current > self.total:
            del self
        self.current += 1
        if self.current % self.stepsToAlert == 0:
            runnedPercent = 1.0 * self.current / self.total
            leftSteps = ( 1 - runnedPercent ) * self.totalSteps
            nowTime = datetime.datetime.now()
            try:
                recentInterval = ProcessAlert.timeInterval(self.lastTime, nowTime)  #( nowTime - self.lastTime ).microseconds
                speed = 1000000.0 * self.stepsToAlert / recentInterval
                self.lastTime = nowTime
            except TypeError:
                self.startTime = nowTime
                self.lastTime = nowTime
                return
            self.recent.append(recentInterval)
            if len(self.recent) > self.careRecentTimes: self.recent = self.recent[1:]
            left_time = ""
            if len(self.recent) == self.careRecentTimes:
                each_step_time = sum(self.recent)/len(self.recent)
                left_time = leftSteps * each_step_time
                left_time = ProcessAlert.microSeconds2Human(left_time)
                print >>self.OUT, "Have read %d lines(%.2f%%)...estimate %s left. %s\n\tCurrent Memory: %s; Speed: %.2f lines/second\n\t%s" % (self.current, 100.0*self.current/self.total, left_time, other, Memory.memory(), speed, now())
            else:
                print >>self.OUT, "Have read %d lines(%.2f%%)... %s\n\tCurrent Memory: %s; Speed: %.2f lines/second\n\t%s" % (self.current, 100.0*self.current/self.total, other, Memory.memory(), speed, now())
        flush()
    def finish(self, other=''):
        nowTime = datetime.datetime.now()
        try:
            total_time = ProcessAlert.timeInterval(self.startTime, nowTime)
        except TypeError:
            pass
        ave_speed = 1000000.0 * self.current / total_time
        huamn_total_time = ProcessAlert.microSeconds2Human(total_time)
        print >>self.OUT, "Task finished seccessfully: total time: %s; average speed: %.2f lines/second. %s\n\t%s" % (total_time, ave_speed, other, now())
        del self



def sortObjDictByAttr(rawDict, keys=[]):
    """ Sort Object Dict By its Attributes
    Example:
        # Class Object is values in dict
        class Hey(object):
            __slots__ = ['chr', 'strand', 'start', 'end']
            def __init__(self):
                import random
                self.chr = ['chr1', 'chr2', 'chr3', 'chrX'][random.randint(0, 3)]
                self.strand = ['+', '-'][random.randint(0, 1)]
                self.start = random.randint(0, 5)
                self.end = random.randint(0, 5)
        # beautiful print dict
        def printDict(Dict):
            for idx in Dict:
                print '=======%s======' % (idx, )
                for Key in ['chr', 'strand', 'start', 'end']:
                    print Key, getattr(Dict[idx], Key)
        #sort dict
        example_dict = {1: Hey(), 2: Hey(), 3: Hey(), 4: Hey()}
        printDict(example_dict)
        sortObjDictByAttr(example_dict, keys=['chr', 'strand', 'start', 'end'])
    """
    dictItems = rawDict.items()
    SortLambda = lambda x: [getattr(x[1], keys[idx]) for idx in range(len(keys))]
    dictItems.sort(key=SortLambda)
    return [ item[0] for item in dictItems ]

def defined(Dict, Key):
    """ judge if the Key in the Dict
    Example:
        Dict = {1:2, 2:9}
        defined(Dict, Key=2)
        defined(Dict, Key=3)
    """
    try:
        Dict[Key]
    except KeyError:
        return False
    else:
        return True

def now():
    "Return Current Time"
    return str(datetime.datetime.now())

def showTime(Doing):
    print >>GLOBAL.LOG, '\nStart to %s\n\t%s' % (Doing, now())

def Complain_Shell_Error(ShellCommand, ShellReturnCode, ShellShowInfo):
    print >>GLOBAL.ERR, "# *******************************"
    print >>GLOBAL.ERR, "An Unexpected ShellCommand Occurred: \n\t%s" % (ShellCommand, )
    print >>GLOBAL.ERR, "Returned Code: %s" % (ShellReturnCode, )
    print >>GLOBAL.ERR, "Showed Info: \n\t%s" % (ShellShowInfo, )
    print >>GLOBAL.ERR, "Process is forced to stop now... \n\t%s " % (now(), )
    print >>GLOBAL.ERR, "# *******************************"
    flush()
    raise Exception("Shell Error")

def runCMD(cmd):
    print >>GLOBAL.LOG, "\t==================Shell Command====================="
    print >>GLOBAL.LOG, "\t"+cmd
    print >>GLOBAL.LOG, "\t===================================================="
    flush()
    return_code, return_info = commands.getstatusoutput( cmd )
    if return_code != 0:
        Complain_Shell_Error(ShellCommand=cmd, ShellReturnCode=return_code, ShellShowInfo=return_info)
    print >>GLOBAL.ERR, return_info
    flush()

def countFileLines(fileName, skipHead=[]):
    """ Count File Lines
    Example:
        countFileLines(fileName="/150T/zhangqf/tmp_lp/PARIS/766_72/Aligned.out.sam", skipHead=['@'])
    """
    import commands
    if not os.path.exists(fileName):
        print >>STD.ERR, "Error: countFileLines ==> %s doesn't exists " % (fileName, )
        raise IOError
    grepLine = ""
    for skip_head in skipHead:
        grepLine += "("+skip_head+")|"
    if grepLine:
        grepLine = "grep -Pv \"^"+grepLine[:-1]+"\""
        CMD = "{} {} | {}".format(grepLine, fileName, "wc -l")
    else:
        CMD = "wc -l "+fileName
    Lines = int(commands.getstatusoutput( CMD )[1].split()[0])
    return Lines


def OpenFile(fileName, mode='r'):
    try:
        FILE_HANDLE = open(fileName, mode)
    except IOError, e:
        print >>GLOBAL.ERR, "Open File {} failed. Error information: \n\t{}".format(fileName, e)
        raise IOError
    flush()
    return FILE_HANDLE

def flush():
    GLOBAL.LOG.flush()
    GLOBAL.ERR.flush()

def sortSupportParallel():
    return_code, return_info = commands.getstatusoutput("sort --help")
    if "--parallel" in return_info:
        return True
    else:
        return False

def sam_attr(samLine):
    """ Parse Attributes from SAM file lines
    Example:
        samLine = "SRR3404926.24301833 256 chr21   8209800 0   25M189522N30M2S *   0   0   AATACATGCCGACGGGCGCTGACCCCAGTCGGTCCTGAGAGATGGGCGAGCGCCGAG   IHIIIIIGGIHHHHHIIIIIHHHIIIIIDGHIHHIIHIGHIIHHIIIIIIHHIIIII   NH:i:10 HI:i:8  AS:i:47"
        samAttr = sam_attr(samLine)
    """
    data = samLine.strip().split()
    attrList = data[11:]
    attrDict = {}
    for attr in attrList:
        attr_items = attr.split(':')
        attrDict[ attr_items[0] ] = int(attr_items[2]) if attr_items[2].isdigit() else attr_items[2]
    return attrDict

def improved_readID(samLine):
    """ replace SRR3404926.24301833 with SRR3404926.24301833_1 1 is the multi-map index
    Example:
        samLine = "SRR3404926.24301833 256 chr21   8209800 0   25M189522N30M2S *   0   0   AATACATGCCGACGGGCGCTGACCCCAGTCGGTCCTGAGAGATGGGCGAGCGCCGAG   IHIIIIIGGIHHHHHIIIIIHHHIIIIIDGHIHHIIHIGHIIHHIIIIIIHHIIIII   NH:i:10 HI:i:8  AS:i:47"
        improved_readID(samLine)
    """
    samAttr = sam_attr(samLine)
    read_id = samLine.split()[0] + '_' + `samAttr['HI']`
    return read_id


def localAlignment(seq1, seq2, AlignmentParameters=None):
    """ Local Align 2 Sequences
    Example:
        localAlignment('ATCGATCGATCGATCGATGCTAGA', 'GTAGCGTAGCATCGATCGATGCAGTATGCTA')
    Return:
        score, aligment, [1-based-start, 1-based-end], [1-based-start, 1-based-end]
    """
    if not AlignmentParameters:
        AlignmentParameters = {
            "matchEnergy": {
                'AA': -999, 'AC': -999, 'AG': -999, 'AT': 1, 'AN': -3,
                'CA': -999, 'CC': -999, 'CG': 1, 'CT': -999, 'CN': -3,
                'GA': -999, 'GC': 1, 'GG': -999, 'GT': 1, 'GN': -3,
                'TA': 1, 'TC': -999, 'TG': 1, 'TT': -999, 'TN': -3,
                'NA': -3, 'NC': -3, 'NG': -3, 'NT': -3, 'NN': 1,
                },
            "gapPenalty": -1
        }
    # initialzation
    seq1_len = len(seq1)
    seq2_len = len(seq2)
    matrix = []
    for i in range(seq2_len):
        matrix.append( [] )
        for j in range(seq1_len):
            matrix[i].append( {} )
            matrix[i][j]['score'] = 0
            matrix[i][j]['pointer'] = 'none'
    # fill
    maxScore = 0; maxI = 0; maxJ = 0;
    for i in range(1, seq2_len):
        for j in range(1, seq1_len):
            # calculate match score
            pair = seq1[j-1] + seq2[i-1]
            diagonalScore = matrix[i-1][j-1]['score'] + AlignmentParameters['matchEnergy'][pair]
            # calculate gap score
            upScore = matrix[i-1][j]['score'] + AlignmentParameters['gapPenalty']
            leftScore = matrix[i][j-1]['score'] + AlignmentParameters['gapPenalty']
            if diagonalScore <= 0 and upScore <= 0 and leftScore <= 0:
                matrix[i][j]['score'] = 0
                matrix[i][j]['pointer'] = 'none'
                continue
            # choose best score
            if diagonalScore >= upScore:
                if diagonalScore >= leftScore:
                    matrix[i][j]['score'] = diagonalScore
                    matrix[i][j]['pointer'] = 'diagonal'
                else:
                    matrix[i][j]['score'] = leftScore
                    matrix[i][j]['pointer'] = 'left'
            else:
                if upScore >= leftScore:
                    matrix[i][j]['score'] = upScore
                    matrix[i][j]['pointer'] = 'up'
                else:
                    matrix[i][j]['score'] = leftScore
                    matrix[i][j]['pointer'] = 'left'
            if matrix[i][j]['score'] > maxScore:
                maxI = i
                maxJ = j
                maxScore = matrix[i][j]['score']
    # trace-back
    align1 = align2 = ''
    i = maxI; j = maxJ;
    while 1:
        if matrix[i][j]['pointer'] == 'none': break
        if matrix[i][j]['pointer'] == 'diagonal':
            align1 += seq1[j-1]
            align2 += seq2[i-1]
            i -= 1; j -= 1;
        elif matrix[i][j]['pointer'] == 'left':
            align1 += seq1[j-1]
            align2 += '-'
            j -= 1
        elif matrix[i][j]['pointer'] == 'up':
            align1 += '-'
            align2 += seq2[i-1]
            i -= 1
    align = align1[::-1] + ':' + align2
    return maxScore, align, j+1, maxJ, i+1, maxI



def parseCigar(cigar, getLargestM=False, getLeadingS=False, getMatchLen=False):
    """ Split Cigar And Get Some Info from Cigar
    Example:
        parseCigar('20S53M181747N14M') # (['S', 'M', 'N', 'M'], [20, 53, 181747, 14])
        parseCigar('20S53M181747N14M', getLargestM=True) # 53
        parseCigar('20S53M181747N14M', getLeadingS=True) # 20
    """
    def remove(List, elem):
        return filter(lambda a: a != elem, List)
    matchSize = [ int(a) for a in remove(re.split('[MIDNSHP=X]', cigar), '') ]
    match = remove(re.split('\d', cigar), '')
    if getLargestM:
        largestM = 0
        for idx in range(len(match)):
            if match[idx] == 'M':
                if matchSize[idx] > largestM:
                    largestM = matchSize[idx]
        return largestM
    elif getLeadingS:
        if match[0] != 'S':
            #print >>GLOBAL.ERR, "Warning: unexpected CIGAR String: %s" % (cigar)
            return 0
        else:
            return matchSize[0]
    elif getMatchLen:
        pass
    else:
        return match, matchSize




