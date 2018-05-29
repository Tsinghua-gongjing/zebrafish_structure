#-*- coding:utf-8 -*-

from Environment import *

def genDuplexGroup(duplexGroup, read, task_name, minOverlap=5, multipleDG=False):
    duplexGroupBedFile = GLOBAL.outDir+"/tmp."+task_name+'.duplexGroup.bed'
    duplexGroupBedFileSortedFile = duplexGroupBedFile + '.sorted'
    # sort file
    if sortSupportParallel():
        #import multiprocessing
        #cores = multiprocessing.cpu_count()
        cmd = "sort --parallel=8 -k1,1 -k6,6 -k7,7 -k12,12 -k2,2n -k3,3n -k8,8n -k9,9n %s -o %s"
    else:
        cmd = "sort -k1,1 -k6,6 -k7,7 -k12,12 -k2,2n -k3,3n -k8,8n -k9,9n %s -o %s"
    runCMD( cmd % (duplexGroupBedFile, duplexGroupBedFileSortedFile) )
    # define a Alert Object
    dgBEDLines = countFileLines(duplexGroupBedFileSortedFile)
    Alerter = ProcessAlert(total=dgBEDLines, stepsToAlert=10000, OUT=GLOBAL.LOG)
    # process BED
    duplexCount = 0
    lineCount = 0; firstPossible = 0
    BED = open(duplexGroupBedFileSortedFile)
    line = BED.readline()
    while line:
        if line.startswith('#'):
            pass
        else:
            lineCount += 1
            Alerter.alert("firstPossible: %d, duplexCount: %d" % (firstPossible, duplexCount))
            ( chr1, start1, end1, id1, score1, strand1, chr2, start2, end2, id2, score2, strand2 ) = line.strip().split()
            lastDGoverlapped = 0
            nonOverlapped = 1
            for idx in range(firstPossible, duplexCount):
                overlapped = checkOverlap( duplexGroup[idx], chr1, strand1, int(start1), int(end1), chr2, strand2, int(start2), int(end2) )
                if overlapped >= minOverlap:
                    nonOverlapped = 0
                    lastDGoverlapped = 1
                    addRead2DuplexGroup ( duplexGroup[idx], id1, int(start1), int(end1), int(start2), int(end2) )
                    if read[int(id1)].clique == -1:
                        read[int(id1)].clique = `idx`
                    else:
                        read[int(id1)].clique += ";"+`idx`
                    if not multipleDG: break
                elif overlapped == -1:
                    ## the start of this read is beyond the firstPossible duplex, so no need to check firstPossible in the future
                    if not lastDGoverlapped:
                        firstPossible = idx + 1
                else:
                    lastDGoverlapped = 1
            if nonOverlapped:
                newDuplexGroup ( duplexGroup, duplexCount, id1, chr1, strand1, int(start1), int(end1), chr2, strand2, int(start2), int(end2) )
                read[int(id1)].clique = `duplexCount`
                duplexCount += 1
        line = BED.readline()
    Alerter.finish()
    BED.close()
    print >>GLOBAL.LOG, "Finally %d duplex group generated from file %s\n\t%s" % (duplexCount, duplexGroupBedFile, now())
    return duplexCount

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




# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 调试
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
"""

duplexGroup = {}
genDuplexGroup(duplexGroup, read, task_name=parameters['task'], minOverlap=parameters['minOverlap'], multipleDG=parameters['multipleDG'])

def beatifulPrintDuplexGroup(dg_Item, title=''):
    print "%22s%s%s" % ("===========", title,"===========")
    for Key in ['chr1', 'strand1', 'start1', 'end1', 'chr2', 'strand2', 'start2', 'end2', 'support', 'reads', 'collapsedTo']:
        print "%22s\t%s" % (Key, getattr(dg_Item, Key))

beatifulPrintDuplexGroup(duplexGroup[4])

for dg_key in duplexGroup:
    KeyStr = getattr(duplexGroup[dg_key], 'reads')
    Support = getattr(duplexGroup[dg_key], 'support')
    if KeyStr.count(';') + 1 != Support:
        print 'Unexpected Error'
    if '50056' in KeyStr:
        beatifulPrintDuplexGroup(duplexGroup[dg_key], dg_key)

support = []
for dg_key in duplexGroup:
    KeyStr = getattr(duplexGroup[dg_key], 'reads')
    Support = getattr(duplexGroup[dg_key], 'support')
    support.append(Support)

support.sort()
for i in range(1, 10):
    print "%d ==> %d" % (i, support.count(i))


# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 恢复
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-

for readID in read:
    read[readID].clique = -1

duplexGroup = {}

"""




















