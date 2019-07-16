#-*- coding:utf-8 -*-

from Environment import *

def old_finalizeReads(read, readmap, duplexGroup):
    dgCount = 0
    # define a Alert Object
    Alerter = ProcessAlert(total=len(readmap), stepsToAlert=100000, OUT=GLOBAL.LOG)
    for name in readmap:
        Alerter.alert()
        if readmap[name] > 0: continue
        readID = 0 - readmap[name]
        if read[readID].clique == -1: continue
        cliques = sorted( [int(a) for a in read[readID].clique.split(';') ] )
        cliques = list(set(cliques)) # remove repeat
        read[ readID ].clique = ';'.join([ `clique` for clique in cliques if duplexGroup[clique].collapsedTo == -1 ])
    Alerter.finish()

def finalizeReads(read, readmap, duplexGroup):
    dgCount = 0
    # define a Alert Object
    Alerter = ProcessAlert(total=len(readmap), stepsToAlert=100000, OUT=GLOBAL.LOG)
    for name in readmap:
        Alerter.alert()
        if readmap[name] > 0: continue
        readID = 0 - readmap[name]
        if read[readID].clique == -1: continue
        cliques = sorted( list(set(read[readID].clique.split(';'))), key=lambda x: int(x) )
        read[ readID ].clique = ';'.join(cliques)
    Alerter.finish()

# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 调试
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-

"""

finalizeReads(read, readmap, duplexGroup)


for read_Key in read:
    if read[ read_Key ].clique != -1 and ':' in read[ read_Key ].clique:
        print read_Key, read[ read_Key ].clique

"""

