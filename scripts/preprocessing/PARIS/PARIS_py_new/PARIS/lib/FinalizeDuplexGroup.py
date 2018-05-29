#-*- coding:utf-8 -*-

from Environment import *

def finalizeDuplexGroup(duplexGroup, read, readmap, task_name):
    readClusterBed = GLOBAL.outDir+"/tmp.%s.readCluster.bed" % (task_name, )
    RC = open(readClusterBed, 'w')
    # define a Alert Object
    Alerter = ProcessAlert(total=len(duplexGroup), stepsToAlert=100000, OUT=GLOBAL.LOG)
    for dg in duplexGroup.keys():
        Alerter.alert()
        if duplexGroup[dg].collapsedTo != -1: continue
        reads = sorted([int(a) for a in duplexGroup[dg].reads.split(';')])
        reads = list(set(reads))  # remove repeat
        duplexGroup[dg].support = 0
        for idx in range(len(reads)):
            duplexGroup[dg].support += 1
            read[reads[idx]].clique += ';'+`dg`
            names = read[reads[idx]].name.split(';')
            for name in names:
                if readmap[name] > 0:
                    if read[reads[idx]].left.chr == read[reads[idx]].right.chr and read[reads[idx]].left.strand == read[reads[idx]].right.strand:
                        print >>RC, '\t'.join([read[reads[idx]].left.chr, `read[reads[idx]].left.start`, `read[reads[idx]].right.end`, `readmap[ name ]`, '1', read[reads[idx]].left.strand ])
                    else:
                        print >>RC, '\t'.join([read[reads[idx]].left.chr, `read[reads[idx]].left.start`, `read[reads[idx]].left.end`, `readmap[ name ]`, '1', read[reads[idx]].left.strand ])
                        print >>RC, '\t'.join([read[reads[idx]].right.chr, `read[reads[idx]].right.start`, `read[reads[idx]].right.end`, `readmap[ name ]`, '2', read[reads[idx]].right.strand ])
                    readmap[name] = 0 - readmap[name]
    Alerter.finish()
    RC.close()


# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 调试
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-

"""

finalizeDuplexGroup(duplexGroup, read, readmap, task_name=parameters['task'])


for dg_key in duplexGroup:
    KeyStr = getattr(duplexGroup[dg_key], 'reads')
    Support = getattr(duplexGroup[dg_key], 'support')
    if KeyStr.count(';') + 1 != Support:
        print 'Unexpected Error'

xx = readmap.values()
sum([1 for i in xx if i < 0])
sum([1 for i in xx if i > 0])




"""