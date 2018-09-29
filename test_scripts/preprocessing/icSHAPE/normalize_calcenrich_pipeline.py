import sys
from optparse import OptionParser
import pandas as pd
import subprocess as sp
import os

def normalizeRT(filename, stdoutFile, stderrFile, winSize=200, step=30):
    normalizer = '/Share/home/zhangqf7/jinsong_zhang/icSHAPE/icSHAPE-master/scripts/normalizeRTfile_bywindow.pl'
    output = filename + "_norm"
    stdout = open(stdoutFile, "w")
    stderr = open(stderrFile, "w")
    args = ['perl', normalizer, "-d", "32", "-l", "32", "-m", "mean:vigintile2", 
            "-o", output, "-i", filename, "-w", str(winSize), "-s", str(step)]
    print >>stderr, "commandLine: %s" % " ".join(args)
    process = sp.Popen(args, stdout=stdout, stderr=stderr)
    #stdout.close()
    #stderr.close()
    return process, stdout, stderr


def calcEnrich(bgFile, fgFile, outputFile, stdoutFile, stderrFile):
    calculator = '/Share/home/zhangqf7/jinsong_zhang/icSHAPE/icSHAPE-master/scripts/calcEnrich_test.pl'
    args = ['perl', calculator, "-f", fgFile, "-b", bgFile, "-o", outputFile, 
            "-w", "factor5:scaling1", "-y", "10", "-x", "0.25", "-e", 'complex']
    stdout = open(stdoutFile, "w")
    stderr = open(stderrFile, "w")
    print >>stderr, "commandLine: %s" % " ".join(args)
    sp.call(args, stdout=stdout, stderr=stderr)
    stdout.close()
    stderr.close()


def filterEnrich(tmpFile, output, bgFile, fgFile, stdoutFile, stderrFile, hitCoverage=2, Basedensity=200):
    Filter = '/Share/home/zhangqf7/jinsong_zhang/icSHAPE/icSHAPE-master/scripts/filterEnrich_test.pl'
    args = ['perl', Filter, "-i", tmpFile, "-o", output, "-f", fgFile, "-b", bgFile, 
            "-s", "5", "-e", "30", "-T", str(hitCoverage), "-t", str(Basedensity)]
    stdout = open(stdoutFile, "w")
    stderr = open(stderrFile, "w")
    print >>stderr, "commandLine: %s" % " ".join(args)
    sp.call(args, stdout=stdout, stderr=stderr)
    stdout.close()
    stderr.close()


def main():
    usage = "Usage: %prog [option] -b backgroundRT -f foregroundRT -p prefix"
    parser = OptionParser(usage)
    parser.add_option("-f", "--foreground", dest="fgFile", action="store") #required
    parser.add_option("-b", "--background", dest="bgFile", action="store") #required
    parser.add_option("-o", "--outputDir", dest="outputDir", action="store")
    parser.add_option("-p", "--prefix", dest="prefix", action="store") #required
    parser.add_option("-w", "--winSize", dest="winSize", action="store", type="int",
                       default=200, help="specify window size to normalize RT.[default: %default]")
    parser.add_option("-s", "--step", dest="step", action="store", type="int",
                      default=30, help="specify step to normalize RT.[default: %default]")
    parser.add_option("-T", dest="hitCoverage", action="store", type="int",
                      default=2, help="specify the hit coverage threshold to filter enrich score.[default: %default]")
    parser.add_option("-t", dest="Basedensity", action="store", type="int",
                      default=200, help="specify the basedensity threshold to filter enrich score.[default: %default]")
    (option, args) = parser.parse_args()
    if (not option.bgFile) or (not option.fgFile) or (not option.prefix):
        print "please specify background and foreground RT file"
        sys.exit(1) 
    bgFile = option.bgFile
    fgFile = option.fgFile
    hitCoverage = option.hitCoverage
    Basedensity = option.Basedensity
    winSize = option.winSize
    step = option.step
    if not option.outputDir:
        outputDir = os.getcwd()
    else:
        outputDir = option.outputDir
    prefix = option.prefix
    bgNormFile = os.path.join(outputDir, bgFile + "_norm.win%s.stp%s" % (winSize, step))
    fgNormFile = os.path.join(outputDir, fgFile + "_norm.win%s.stp%s" % (winSize, step))
    bgProcess, bgout, bgerr = normalizeRT(bgFile, stdoutFile=os.path.join(outputDir, prefix + ".bg.out"), stderrFile=os.path.join(outputDir, prefix + ".bg.err"), winSize=winSize, step=step)
    fgProcess, fgout, fgerr = normalizeRT(fgFile, stdoutFile=os.path.join(outputDir, prefix + ".fg.out"), stderrFile=os.path.join(outputDir, prefix + ".fg.err"), winSize=winSize, step=step)
    bgProcess.wait()
    print "background RT finished"
    fgProcess.wait()
    print "foreground RT finished"
    bgout.close()
    bgerr.close()
    fgout.close()
    fgerr.close()
    tmpFile = os.path.join(outputDir, prefix + ".icshape.w%s.s%s.tmp.out" % (winSize, step))
    calcEnrich(bgNormFile, fgNormFile, tmpFile, stdoutFile=os.path.join(outputDir, prefix + ".calcenrich.out"), stderrFile=os.path.join(outputDir, prefix + '.calcenrich.err'))
    FinalEnrichFile = os.path.join(outputDir, prefix + ".icshape.w%s.s%s.T%s.t%s.out" % (winSize, step, hitCoverage, Basedensity))
    filterEnrich(tmpFile, FinalEnrichFile, bgFile, fgFile, stdoutFile=os.path.join(outputDir, prefix + ".filter.out"), stderrFile=os.path.join(outputDir, prefix + ".filter.err"), 
                 hitCoverage=hitCoverage, Basedensity=Basedensity)


if __name__ == '__main__':
    main()
