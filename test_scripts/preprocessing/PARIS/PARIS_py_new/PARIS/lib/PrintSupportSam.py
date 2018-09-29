#-*- coding:utf-8 -*-

from Environment import *

def printSupportSam(outputFile, allSupportSam, duplexGroup, read, readmap, outputReadFile=False, minSupport=2):
    if outputReadFile:
        READ = open(outputReadFile, 'w')
    OUT = open(outputFile, 'w')
    samFiles = allSupportSam.split(':')
    lineCount = 0; validCount = 0; noDefineInReadmap = 0; unParsedJunctionLine = 0; outputSamLine = 0
    firstSam = True
    # define a Alert Object
    #samFiles.reverse()
    totalSamLines = sum([countFileLines(samFile, skipHead=['@']) for samFile in samFiles])
    Alerter = ProcessAlert(total=totalSamLines, stepsToAlert=100000, OUT=GLOBAL.LOG)
    last_read_id = ""; last_line = ""
    for samFile in samFiles:
        SAM = open(samFile)
        print >>GLOBAL.LOG, "check for supporting reads from sam file %s...\n\t%s" % (samFile, now())
        line = SAM.readline()
        while line:
            if line.startswith('#'):
                pass
            elif line.startswith('@'):
                if firstSam:
                    print >>OUT, line.strip()
            else:
                firstSam = False
                lineCount += 1
                Alerter.alert("valid line: %d" % (validCount, ))
                data = line.strip().split()
                improved_read_id = improved_readID(line)
                if defined(readmap, data[0]):
                    read_id = data[0]
                elif defined(readmap, improved_read_id):
                    read_id = improved_read_id
                else:
                    print >>GLOBAL.ERR, "Warning: sam line not defined in readmap. Possible reason: preserveMulti or removeRedundancy is on. File: %s\n\t%s" % (samFile, line.strip())
                    noDefineInReadmap += 1
                    line = SAM.readline()
                    continue
                if readmap[read_id] > 0:
                    line = SAM.readline()
                    continue
                readID = 0 - readmap[read_id]
                realReadID = read[readID].collapsedTo if read[readID].collapsedTo != -1 else readID
                if read[realReadID].cigar != -1:
                    # only junction reads have cigar attr
                    isChiastic, cigar = read[realReadID].cigar.split(':')
                    isChiastic = int(isChiastic)
                    if isChiastic == 0:
                        # preserve: [========]~~~~~~~~~(sense: 0, antisense: 272) remove: ~~~~~~~~~[========](sense: 256, antisense: 16)
                        #if int(data[1]) == 256 or int(data[1]) == 16:
                        # Modified 2017-7-7
                        leadingSeq = parseCigar(data[5], getLeadingS=True)
                        if leadingSeq >= 15:
                            line = SAM.readline()
                            continue
                        data[1] = `int(data[1]) & 3839`
                        data[5] = cigar
                    elif isChiastic == 1:
                        # remove: [========]~~~~~~~(sense: 0, antisense: 272)   preserve: ~~~~~~~~[=========](sense: 256, antisense: 16)
                        # Modified 2017-7-7
                        if data[0] == last_read_id:
                            last_data = last_line.strip().split()
                            last_leadingS = parseCigar(last_data[5], getLeadingS=True)
                            this_leadingS = parseCigar(data[5], getLeadingS=True)
                            if last_leadingS > this_leadingS:
                                data = last_data
                                #line = last_data
                        else:
                            last_line = line.strip()
                            last_read_id = line.strip().split()[0]
                            line = SAM.readline()
                            continue
                        #if int(data[1]) == 0 or int(data[1]) == 272:
                        #    line = SAM.readline()
                        #    continue
                        data[1] = `int(data[1]) & 3839`
                        data[9] = reverseRead( data[5], data[9] )
                        data[10] = reverseRead( data[5], data[10] )
                        if not data[9] or not data[10]:
                            unParsedJunctionLine += 1
                            print >>GLOBAL.ERR, "Unexpected Line:\n\tlast line: %s\n\tthis line: %s\n\tlast leading S: %s this leading S: %s" % (last_line, line, last_leadingS, this_leadingS)
                            line = SAM.readline()
                            continue
                        data[5] = cigar
                    elif isChiastic == 2:
                        # reads map to different chromosomes or strands
                        # a reads may output to outputReadFile twice
                        pass
                    preString = "\t".join(data)+'\tXG:i:'+str(isChiastic)
                    #OUT.writelines()
                else:
                    preString = line.strip()+"\tXG:i:0"
                    #OUT.writelines()
                dgIDs = read[realReadID].clique.split(';')
                willGoOn = False
                for dgID in dgIDs:
                    if duplexGroup[int(dgID)].support >= minSupport:
                        willGoOn = True
                        break
                if not willGoOn:
                    line = SAM.readline()
                    continue
                Tag = "\tDG:i:" + read[realReadID].clique
                if read[realReadID].ngTag != -1: Tag += "\tNG:i:"+`read[realReadID].ngTag` # print >>OUT, "\tNG:i:"+`read[realReadID].ngTag`
                outputSamLine += 1
                print >>OUT, preString+Tag
                try:
                    # not all read have left/right attr
                    if outputReadFile:
                        print >>READ, "\t".join([ read[realReadID].left.chr, `read[realReadID].left.start`, `read[realReadID].left.end`, read[realReadID].name, `realReadID`, read[realReadID].left.strand ])
                        print >>READ, "\t".join([ read[realReadID].right.chr, `read[realReadID].right.start`, `read[realReadID].right.end`, read[realReadID].name, `realReadID`, read[realReadID].right.strand ])
                    validCount += 1
                except AttributeError:
                    print >>GLOBAL.ERR, "%d information is not found in read.\n\t%s" % (realReadID, line.strip())
            line = SAM.readline()
        SAM.close()
        print >>GLOBAL.LOG, "in total %d lines read from sam files." % (lineCount, )
        print >>GLOBAL.LOG, "\tamong which %d lines generate supports for valid base pairing. %d lines undefined in readmap\n\tTime: %s" % (validCount, noDefineInReadmap, now()) 
    Alerter.finish("outputSamLine: %s; unParsedJunctionLine: %s" % (outputSamLine, unParsedJunctionLine))
    OUT.close()
    if outputReadFile: READ.close()
    return lineCount, validCount


def reverseRead(cigar, seq):
    reverseSeq = ""
    leadingSeq = parseCigar(cigar, getLeadingS=True)
    if not leadingSeq:
        return None
    reverseSeq = seq[leadingSeq:]
    reverseSeq += seq[:leadingSeq]
    return reverseSeq





# =-==-=-=-=-=-=--=-=-=-=-=-=-=-
# 调试
# =-==-=-=-=-=-=--=-=-=-=-=-=-=-

"""

printSupportSam(supportSamFile, allSupportSam, read, readmap, outputReadFile=supportReadFile)

xx = readmap.values()
sum([i for i in readmap if readmap[i] < 0])
sum([1 for i in xx if i > 0])

"""



