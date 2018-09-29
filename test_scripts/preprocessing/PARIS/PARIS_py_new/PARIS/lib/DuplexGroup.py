#-*- coding:utf-8 -*-

from Environment import *

def gDuplexGroup(task_name, duplexBin, logFile, errFile, minOverlap=5, multipleDG=False, maxGap=10, maxTotal=30):
    """
        Using C++ Code to Generate and Collapse Duplex Group
    """
    readPairFile = GLOBAL.outDir+"/tmp."+task_name+'.duplexGroup.bed'
    sortedReadPairFile = readPairFile + '.sorted'
    dgFile = GLOBAL.outDir+"/tmp."+task_name+'.gen_duplexGroup.bed'
    sortedDgFile = GLOBAL.outDir+"/tmp."+task_name+'.gen_duplexGroup.bed.sorted'
    collapsedDgFile = GLOBAL.outDir+"/tmp."+task_name+'.collapsed_duplexGroup.bed'
    # sort file
    if sortSupportParallel():
        cmd = "sort --parallel=8 -k1,1 -k6,6 -k7,7 -k12,12 -k2,2n -k3,3n -k8,8n -k9,9n %s -o %s"
    else:
        cmd = "sort -k1,1 -k6,6 -k7,7 -k12,12 -k2,2n -k3,3n -k8,8n -k9,9n %s -o %s"
    runCMD( cmd % (readPairFile, sortedReadPairFile) )
    # generate duplex group
    cmd_genDG = "%s GenDuplexGroup --log %s --error %s --inReadPairFile %s --outDGFile %s --minOverlap %s --multipleDG %s"
    runCMD(cmd_genDG % (duplexBin, logFile, errFile, sortedReadPairFile, dgFile, minOverlap, 'yes' if multipleDG else 'no'))
    # sort file
    if sortSupportParallel():
        cmd = "sort --parallel=8 -k1,1 -k4,4 -k5,5 -k8,8 -k2,2n -k3,3n -k6,6n -k7,7n %s -o %s"
    else:
        cmd = "sort -k1,1 -k4,4 -k5,5 -k8,8 -k2,2n -k3,3n -k6,6n -k7,7n %s -o %s"
    runCMD( cmd % (dgFile, sortedDgFile) )
    # collapse duplex group
    cmd_collapseDG = "%s CollapseDuplexGroup --log %s --error %s --inDGFile %s --outCollapseFile %s --maxGap %s --maxDGOverhang %s"
    runCMD(cmd_collapseDG % (duplexBin, logFile, errFile, sortedDgFile, collapsedDgFile, maxGap, maxTotal))

