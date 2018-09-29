#-*- coding:utf-8 -*-

from Environment import *

def nonOverlappingTag(read, task_name):
    readClusterBedFile = GLOBAL.outDir+"/tmp.%s.readCluster.bed" % (task_name, )
    sortedReadClusterBedFile = readClusterBedFile + ".sorted";
    uniqReadClusterBedFile = readClusterBedFile + ".uniq";
    readClusterFile = readClusterBedFile + ".cluster";
    # sort file
    if sortSupportParallel():
        cmd = "sort --parallel=8 -k1,1 -k6,6 -k2,2n -k3,3n %s -o %s"
    else:
        cmd = "sort -k1,1 -k6,6 -k2,2n -k3,3n %s -o %s"
    runCMD( cmd % (readClusterBedFile, sortedReadClusterBedFile) )
    uniqBed ( sortedReadClusterBedFile, uniqReadClusterBedFile)
    cmd_cluster = "bedtools cluster -i %s -s > %s"
    runCMD( cmd_cluster % (uniqReadClusterBedFile, readClusterFile) )
    # generate proper tags for reads in read
    totalLines = countFileLines(readClusterFile)
    Alerter = ProcessAlert(total=totalLines, stepsToAlert=10000, OUT=GLOBAL.LOG)
    CLUSTER = open(readClusterFile)
    line = CLUSTER.readline()
    ngTag = 0; clusterID = 0
    while line:
        Alerter.alert()
        data = line.strip().split()
        reads = data[3].split(';')
        if clusterID == int(data[6]):
            for read_id in reads:
                read[int(read_id)].cluster = clusterID
                ngTag += 1
                read[int(read_id)].ngTag = ngTag
        else:
            clusterID = int(data[6])
            ngTag = 0
            for read_id in reads:
                read[int(read_id)].cluster = clusterID
                ngTag += 1
                read[int(read_id)].ngTag = ngTag
        line = CLUSTER.readline()
    Alerter.finish()
    CLUSTER.close()

def uniqBed(inputBed, outputBed):
    bedPos = ""; tag = ""; bedStrand = ""; count = 0
    IN = open(inputBed); OUT = open(outputBed, 'w')
    line = IN.readline()
    while line:
        data = line.strip().split()
        tmpPos = "\t".join([data[0], data[1], data[2]])
        if tmpPos != bedPos or data[5] != bedStrand:
            if bedPos:
                print >>OUT, "\t".join([bedPos, tag, `count`, bedStrand])
            bedPos = tmpPos
            bedStrand = data[5]
            tag = data[3]
            count = int(data[4])
        else:
            tag += ";" + data[3]
            count += int(data[4])
        if bedPos:
            print >>OUT, "\t".join([bedPos, tag, `count`, bedStrand])
        line = IN.readline()
    IN.close()
    OUT.close()




