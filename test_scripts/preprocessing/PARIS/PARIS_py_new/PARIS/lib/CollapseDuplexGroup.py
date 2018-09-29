#-*- coding:utf-8 -*-

from Environment import *

def collapseDuplexGroup(duplexGroup, read, maxGap=10, maxTotal=30):
    merged_dgCount = 0
    firstPossible = 0
    dgArray = []
    # define a Alert Object
    Alerter = ProcessAlert(total=len(duplexGroup), stepsToAlert=1000, OUT=GLOBAL.LOG)
    # start to collapse
    for dg in sortObjDictByAttr(duplexGroup, ['chr1', 'strand1', 'chr2', 'strand2', 'start1', 'end1', 'start2', 'end2']):
        dgArray.append(dg)
        lastDGoverlapped = 0
        for idx in range(firstPossible, len(dgArray)-1):
            overlapped = checkOverlapDG(duplexGroup[dgArray[idx]], duplexGroup[dg], read=read, maxGap=maxGap, maxTotal=maxTotal)
            if overlapped == -1:
                if not lastDGoverlapped:
                    firstPossible = idx + 1
            elif overlapped > 0:
                lastDGoverlapped = 1
                mergeDuplexGroup(duplexGroup[dgArray[idx]], duplexGroup[dg])
                duplexGroup[dg].collapsedTo = dgArray[idx]
                merged_dgCount += 1
                break
            else:
                lastDGoverlapped = 1
        Alerter.alert("merged dg: %s" % (merged_dgCount, ))
    Alerter.finish()

def checkOverlapDG(dg_item1, dg_item2, read, maxGap=10, maxTotal=30, checkReads=False):
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
    if overlap == 1 and checkReads:
        # a unmatured option
        read1s = dg_item1.reads.split(';')
        read2s = dg_item2.reads.split(';')
        minOverlap1 = 0
        minOverlap2 = 0
        overlapCount1 = 0
        for r1 in read1s:
            overlapCount2 = 0
            for r2 in read2s:
                overlapCount2 += readsOverlap( read, int(r1), int(r2) )
                if overlapCount2 > minOverlap2:
                    overlapCount1 += 1
                    break
            if overlapCount1 > minOverlap1:
                return 1
    return overlap

def readsOverlap(read, rid1, rid2):
    overlap = read[rid1].left.start < read[rid2].left.end and read[rid1].left.end > read[rid2].left.start and read[rid1].right.start < read[rid2].right.end and read[rid1].right.end > read[rid2].right.start
    return 1 if overlap else 0

def mergeDuplexGroup(dg_item1, dg_item2):
    dg_item1.reads += ';'+dg_item2.reads
    dg_item1.start1 = min(dg_item1.start1, dg_item2.start1)
    dg_item1.end1 = max(dg_item1.end1, dg_item2.end1)
    dg_item1.start2 = min(dg_item1.start2, dg_item2.start2)
    dg_item1.end2 = max(dg_item1.end2, dg_item2.end2)





# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 调试
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-

"""

collapseDuplexGroup(duplexGroup, read, maxGap=10, maxTotal=30)

xx = sortObjDictByAttr(duplexGroup, ['chr1', 'strand1', 'chr2', 'strand2', 'start1', 'end1', 'start2', 'end2'])
xx[:10]
[4707, 4708, 4709, 4710, 4711, 4712, 4713, 4714, 4715, 4716]
beatifulPrintDuplexGroup(duplexGroup[xx[4732-2]])


collapsedTo = 0
for dg_key in duplexGroup:
    dg_collapsedTo = getattr(duplexGroup[dg_key], 'collapsedTo')
    if dg_collapsedTo != -1:
        collapsedTo += 1
        beatifulPrintDuplexGroup(duplexGroup[dg_key], dg_key)
        print ''
        beatifulPrintDuplexGroup(duplexGroup[dg_collapsedTo], dg_collapsedTo)


beatifulPrintDuplexGroup(duplexGroup[27], 27)
beatifulPrintDuplexGroup(duplexGroup[28], 28)
beatifulPrintDuplexGroup(duplexGroup[29], 29)

"""

# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 检查被collapsedTo的dg会不会自己也collapsedTo
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-



