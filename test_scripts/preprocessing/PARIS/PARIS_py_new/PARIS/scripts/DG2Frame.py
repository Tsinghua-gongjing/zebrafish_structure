#-*- coding:utf-8 -*- 
#!env python

import re, sys, os, getopt, time, datetime
import genome2Gene, genome2Trans, ParseTrans

Usage = """
## --------------------------------------
Parse Gene and Transcript of DG Regions
 
Command:
%s -i DGFile -g GFF3/GTF -o output_prefix --genomeCoor genomeCoorFile -s [gencode|ncbi]

# what it is:
 -i                 DG file, multiple file seperated by :, equal to -o
 -g                 genome annotation
 -o                 output file prefix, , multiple file seperated by :, equal to -i 
 -s                 <gencode/ncbi> data source
 --genomeCoor       genomeCoor.txt file produced by parseBedFromGTF.py
# Attention
  Strongly recommend to set --removeTransVersion no --removeGeneVersion no --removeScaffold no
  to produce genomeCoorFile

""" % (sys.argv[0], )

def init():
    params = { }
    opts, args = getopt.getopt(sys.argv[1:], 'hg:o:s:i:', ['genomeCoor='])
    for op, value in opts:
        if op == '-h':
            print Usage;
            sys.exit(-1)
        elif op == '-i':
            mutipleDG = value.split(':')
            params['dg'] = [ os.path.abspath(dgFile) for dgFile in mutipleDG]
        elif op == '-g':
            params['gtf'] = os.path.abspath(value)
        elif op == '-o':
            multipleOut = value.split(':')
            params['output'] = [ os.path.abspath(out).rstrip(".") for out in multipleOut ]
        elif op == '-s':
            if value == 'ncbi':
                params['source'] = 'NCBI'
            elif value == 'gencode':
                params['source'] = 'Gencode'
            else:
                print "Error: -s <gencode/ncbi>"
                sys.exit(-1)
            #assert value in ('ncbi', 'gencode'), 'Error: -t <gencode/ncbi>'
        elif op == '--genomeCoor':
            params['genomeCoor'] = os.path.abspath(value)
    if 'dg' not in params:
        print 'Error: confirm your -i'
        print Usage
        sys.exit(-1)
    if 'gtf' not in params or not os.path.isfile(params['gtf']):
        print 'Error: confirm your -g'
        print Usage
        sys.exit(-1)
    if 'output' not in params:
        print 'Error: specify -o'
        print Usage
        sys.exit(-1)
    if 'source' not in params:
        print 'Error: specify -s'
        print Usage
        sys.exit(-1)
    if 'genomeCoor' not in params or not os.path.isfile(params['genomeCoor']):
        print 'Error: confirm your --genomeCoor'
        print Usage
        sys.exit(-1)
    if len(params['output']) != len(params['dg']):
        print "Error: -i numbers of input file must equal to -o out files"
        sys.exit(-1)
    return params



#sys.argv = ['bin.py', '-i', '/150T/zhangqf/tmp_lp/PARIS/tanglei/59_24/59_24_DG', '-o', '/150T/zhangqf/tmp_lp/PARIS/tanglei/59_24/analysis/59_24', '--genomeCoor', '/Share/home/zhangqf8/lipan/paris/human_virus/GRCh38.genomeCoor.bed', '-s', 'gencode', '-g', '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Gencode/GRCh38.gtf']

def main():
    params = init()
    print "DG ==> Frame..."
    for (inFile, outFile) in zip(params['dg'], params['output']):
        DG2Frame(inFile, outFile+".frame")
    g2g = genome2Gene.genome2GeneClass(params['gtf'], source=params['source'])
    g2t = genome2Trans.genome2TransClass(params['gtf'], source=params['source'])
    parseTrans = ParseTrans.ParseTransClass(params['genomeCoor'])
    print "Frame ==> gFrame..."
    for outFile in params['output']:
        Frame2GeneFrame(outFile+".frame", outFile+".gframe", g2g)
    print "Frame ==> tFrame..."
    for outFile in params['output']:
        Frame2TransFrame(outFile+".frame", outFile+".tframe", g2t, parseTrans)





def DG2Frame(dfFile, frameFile):
    import re
    pattern = "Group\s(\d+).*position\s([\w.\d]+)\(([+-])\):(\d+)-(\d+)\|([\w.\d]+)\(([+-])\):(\d+)-(\d+).*support\s(\d+).*left\s(\d+).*right\s(\d+).*score\s([\d\.]*\d)"
    DG = open(dfFile)
    FRAME = open(frameFile, 'w')
    print >>FRAME, "#Group\tlchr\tlstrand\tlstart\tlend\trchr\trstrand\trstart\trend\tsupport\tlcount\trcount\tscore"
    line = DG.readline()
    while line:
        if line.startswith('Group'):
            data = re.findall(pattern, line)
            print >>FRAME, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % data[0]
        line = DG.readline()
    FRAME.close()
    DG.close()


def Frame2GeneFrame(frameFile, geneFrameFile, genome2GeneHandle):
    def getGeneString(Chr, Start, End, Strand, genome2GeneHandle):
        Genes = genome2GeneHandle.genomeCoor2GeneCoor(Chr, Start, End, Strand, pureGeneID=False)
        Str = ""
        for gene in Genes:
            Str += "||" if Str else ""
            geneData = gene[3].split('_')
            geneID = geneData[0]; geneName = "_".join(geneData[1:])
            Str += "%s:%s-%s:%s:%s:%d-%d" % (gene[0], gene[1], gene[2], geneID, geneName, gene[4], gene[5])
        return Str
    import datetime
    FRAME = open(frameFile)
    GFRAME = open(geneFrameFile, 'w')
    print >>GFRAME, "# From: %s; Date: %s; Annotation: %s" % (frameFile, str(datetime.datetime.now()), genome2GeneHandle.anno_file)
    line = FRAME.readline()
    while line:
        if line.startswith('#'):
            print >>GFRAME, line.strip() + "\tlgenes\trgenes"
            pass
        else:
            data = line.strip().split()
            lString = rString = ""
            try:
                lString = getGeneString(data[1], int(data[3]), int(data[4]), data[2], genome2GeneHandle)
            except KeyError:
                pass
            try:
                rString = getGeneString(data[5], int(data[7]), int(data[8]), data[6], genome2GeneHandle)
            except KeyError:
                pass
            lString = lString if lString else "NULL"
            rString = rString if rString else "NULL"
            print >>GFRAME, line.strip() + "\t" + lString + "\t" + rString
        line = FRAME.readline()
    GFRAME.close()
    FRAME.close()


def Frame2TransFrame(frameFile, transFrameFile, genome2TransHandle, parseTrans):
    def getTransString(Chr, Start, End, Strand, genome2TransHandle, parseTrans):
        Trans = genome2TransHandle.genomeCoor2TransCoor(Chr, Start, End, Strand, pureTransID=False)
        Str = ""
        for trans in Trans:
            Str += "||" if Str else ""
            #print trans
            Str += "%s:%s-%s:%s:%s-%s" % tuple(trans)
            transID = trans[3]; transStart = trans[4]; transEnd = trans[5]
            try:
                transFeature = parseTrans.getTransFeature(transID=transID, showAttr=False)
                Str += ":%s:%s:%s-%s:%s" % (transFeature['gene_name'], transFeature['gene_type'], transFeature['cds_start'], transFeature['cds_end'], transFeature['trans_len'])
            except KeyError:
                Str += ":NULL:NULL:NULL-NULL:NULL"
        return Str
    import datetime
    FRAME = open(frameFile)
    RFRAME = open(transFrameFile, 'w')
    print >>RFRAME, "# From: %s; Date: %s; Annotation: %s" % (frameFile, str(datetime.datetime.now()), genome2TransHandle.anno_file)
    line = FRAME.readline()
    while line:
        if line.startswith('#'):
            print >>RFRAME, line.strip() + "\tltrans\trtrans"
            pass
        else:
            data = line.strip().split()
            lString = rString = ""
            try:
                lString = getTransString(data[1], int(data[3]), int(data[4]), data[2], genome2TransHandle, parseTrans)
            except KeyError:
                pass
            try:
                rString = getTransString(data[5], int(data[7]), int(data[8]), data[6], genome2TransHandle, parseTrans)
            except KeyError:
                pass
            lString = lString if lString else "NULL"
            rString = rString if rString else "NULL"
            print >>RFRAME, line.strip() + "\t" + lString + "\t" + rString
        line = FRAME.readline()
    RFRAME.close()
    FRAME.close()


main()


