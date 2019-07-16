
#-*- coding:utf-8 -*-

from Environment import *

def loadDuplexGroup(duplexGroup, read, readmap, task_name):
    collapsedDgFile = GLOBAL.outDir+"/tmp.%s.collapsed_duplexGroup.bed" % (task_name, )
    # output readCluster.bed file
    readClusterBedFile = GLOBAL.outDir+"/tmp.%s.readCluster.bed" % (task_name, )
    CLUSTER = open(readClusterBedFile, 'w')
    # count File Lines
    totalLines = countFileLines(collapsedDgFile)
    # Define an Alert
    Alerter = ProcessAlert(total=totalLines, stepsToAlert=100000, OUT=GLOBAL.LOG)
    # read collapsed duplex group file
    COLLAPSE = open(collapsedDgFile)
    idx = 0
    line = COLLAPSE.readline()
    while line:
        idx += 1
        Alerter.alert()
        (chr1, start1, end1, strand1, chr2, start2, end2, strand2, supportReads, support) = line.strip().split()
        dg = DuplexGroup()
        # uniq reads
        supportReads = sorted(list(set(supportReads.split(';'))), key=lambda x: int(x))
        dg.support = len(supportReads); dg.reads = ";".join(supportReads) 
        dg.chr1 = chr1; dg.start1 = int(start1); dg.end1 = int(end1); dg.strand1 = strand1
        dg.chr2 = chr2; dg.start2 = int(start2); dg.end2 = int(end2); dg.strand2 = strand2
        # add to duplexGroup
        duplexGroup[idx] = dg
        # tag reads's clique
        for read_id in supportReads:
            currentRead = read[int(read_id)]
            if currentRead.clique == -1: currentRead.clique = `idx`
            else: currentRead.clique += ';'+`idx`
            read_names = currentRead.name.split(';')
            for name in read_names:
                if readmap[name] > 0:
                    if currentRead.left.chr == currentRead.right.chr and currentRead.left.strand == currentRead.right.strand:
                        print >>CLUSTER, '\t'.join( [currentRead.left.chr, `currentRead.left.start`, `currentRead.right.end`, `readmap[ name ]`, '1', currentRead.left.strand ] )
                    else:
                        print >>CLUSTER, '\t'.join([currentRead.left.chr, `currentRead.left.start`, `currentRead.left.end`, `readmap[ name ]`, '1', currentRead.left.strand ])
                        print >>CLUSTER, '\t'.join([currentRead.right.chr, `currentRead.right.start`, `currentRead.right.end`, `readmap[ name ]`, '2', currentRead.right.strand ])
                    readmap[name] = 0 - readmap[name]
        line = COLLAPSE.readline()
    Alerter.finish()
    CLUSTER.close()
    COLLAPSE.close()

