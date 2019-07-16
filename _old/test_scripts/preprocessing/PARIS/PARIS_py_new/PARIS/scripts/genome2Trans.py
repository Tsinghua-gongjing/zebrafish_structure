#-*- coding:utf-8 -*-


"""
给出一段基因组的区域，返回这个区域中转录本对应的坐标

  使用方法：  
    1. ensembl/gencode转换
        gt = genome2TransClass("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Gencode/GRCh37.gtf")
        gt.genomeCoor2TransCoor('chr2', 200000, 350000, '+', pureTransID=False)
    2. ncbi hg19转换
        ncbi_gt = genome2TransClass("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/ncbi/GRCh37.p13.gff3", source='NCBI')
        ncbi_gt.genomeCoor2TransCoor('chrX', 200000, 350000, '+', pureTransID=False)
    3. ensembl/gencode mm10
        gt = genome2TransClass("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.GRCm38.81.gtf")
        gt.genomeCoor2TransCoor('chr9', 13245741-10000, 13245741+1000000, '-', pureTransID=False)

  两个常用的GTF文件：
    /Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.gtf
    /Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.GRCm38.81.gtf

输入：
  注意这里的输入是左闭右开，0-based坐标； 

返回：
  包含转录本区段的列表，每一个元素格式是这样的：[chrID, genomeStart, genomeEnd, transID, transStart, transEnd]
  其中返回的转录本的起点和终点都是1-based的，且前闭后闭; 基因组的起点和终点也是1-based的，前闭后闭

异常：
  当染色体不在GTF文件中时，发出一个KeyError的异常，注意在使用时进行捕获

Version:
    2017-5-29 使用新定义的GTF/GFF3读取函数
"""


import re, sys, os, getopt, time, datetime

def now():
    "返回当前时间字符串"
    return str(datetime.datetime.now())


def binize(RNA_Container, bw=100000):
    "把基因组上的区域分段，加速查找过程。bw越小，占内存越高，越快"
    def statisticPsuedoChrSize(RNA_Container):
        """获取染色体的最大计数( <= 染色体的长度)
            RNA_Container：{rnaID1 => {chr, start, end, strand...}...}
        """
        ChrLen = {}
        for rnaID in RNA_Container:
            Chr = RNA_Container[rnaID]['chr']
            End = int(RNA_Container[rnaID]['end'])
            try:
                ChrLen[Chr] = End if End > ChrLen[Chr] else ChrLen[Chr]
            except KeyError:
                ChrLen[Chr] = End
        return ChrLen
    print "Binize genome to speed up searching.\n\t%s" % (now(), )
    chr_size = statisticPsuedoChrSize(RNA_Container)
    Bin = {}
    for Chr in chr_size.keys():
        Bin[Chr] = {}; Bin[Chr]['+'] = {}; Bin[Chr]['-'] = {}
        count = int( chr_size[Chr]/bw ) + 1
        for idx in range(count):
            Bin[Chr]['+'][idx] = []; Bin[Chr]['-'][idx] = []
    for rnaID in RNA_Container:
        Chr = RNA_Container[rnaID]['chr']
        Strand = RNA_Container[rnaID]['strand']
        Start = int( int(RNA_Container[rnaID]['start']) / bw )
        End = int( int(RNA_Container[rnaID]['end']) / bw )
        for idx in range(Start, End+1):
            Bin[Chr][Strand][idx].append( rnaID )
    return Bin

def genomeRange2TransCoor(ref_bin, GFF_Container, Chr, Start, End, Strand, pureChrSym=False, bw=100000):
    "左闭右开，0-based"
    idxBin = int(Start/bw)
    overlapRegions = []
    while idxBin <= int(End/bw):
        if idxBin > ref_bin[Chr][Strand].keys()[-1]:
            break
        for rnaID in ref_bin[Chr][Strand][idxBin]:
            if End < int(GFF_Container['RNA'][rnaID]['start']) or Start > int(GFF_Container['RNA'][rnaID]['end']):
                continue
            #for transcriptID in GFF_Container['gene_info'][rnaID]['transcript']:
            overlapRegion = overlapTrans(GFF_Container, rnaID, Strand, Start, End)
            for item in overlapRegion: item.insert(0, Chr)
            if len(overlapRegion) > 0:
                overlapRegions += overlapRegion
        idxBin += 1
    return overlapRegions

def overlapTrans(GFF_Container, rnaID, strand, absStart, absEnd):
    "查找一个基因组区段与某个转录本的交集"
    absStart += 1
    transcript = GFF_Container['RNA'][rnaID]
    numExon = len( GFF_Container['exon'][rnaID] )
    relExonStart = 0; startPosInExon = 0; endPosInExon = 0; absStartInExon = 0; absEndInExon = 0;
    overlapRegions = []
    for idxExon in range(numExon):
        exon = GFF_Container['exon'][rnaID][idxExon]
        exonStart = int(exon['start'])
        exonEnd = int(exon['end'])
        exonLength = exonEnd - exonStart + 1
        if exonStart <= absEnd and exonEnd >= absStart:
            "有交集"
            if strand == '+':
                startPosInExon = absStart - exonStart + 1
                if startPosInExon < 1: startPosInExon = 1
                endPosInExon = absEnd - exonStart + 1
                if endPosInExon > exonLength: endPosInExon = exonLength
            elif strand == '-':
                startPosInExon = exonEnd - absEnd + 1
                if startPosInExon < 1: startPosInExon = 1
                endPosInExon = exonEnd - absStart + 1
                if endPosInExon > exonLength: endPosInExon = exonLength
            relStart = relExonStart + startPosInExon
            relEnd = relExonStart + endPosInExon
            absStartInExon = absStart if absStart >= exonStart else exonStart
            absEndInExon = exonEnd if absEnd >= exonEnd else absEnd
            "add to overlapRegions"
            overlapRegions.append( [absStartInExon, absEndInExon, rnaID, relStart, relEnd] )
        relExonStart += exonLength
    return overlapRegions


"""测试函数
import NCBI_Genome
GFF_File = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/ncbi/GRCh38.p7.gff3'
NCBI_Genome_Anno = NCBI_Genome.NCBI_Genome_Class(GFF_File)
NCBI_GFF_Container = NCBI_Genome_Anno.gff3_container
binize(NCBI_GFF_Container['gene_info'], NCBI_GFF_Container['chr_size'], bw=100000)

ChrLen = statisticPsuedoChrSize(NCBI_GFF_Container['RNA'])

Bin = binize(NCBI_GFF_Container['RNA'], bw=100000)

genomeRange2TransCoor(Bin, NCBI_GFF_Container, 'NC_000001.11', Start=17506, End=17706, Strand='-', bw=100000)
"""

class genome2TransClass(object):
    """把基因组的坐标转成转录组的坐标的类"""
    def __init__(self, anno_file, source='Ensembl', bin_width=100000):
        """
        anno_file: 注释文件，Ensembl/Gencode为GTF文件；NCBI为GFF3文件
        source: 注明注释文件的来源(Ensembl/Gencode/NCBI)
        """
        self.anno_file = anno_file
        self.source = source
        self.bw = bin_width
        if source not in ('Ensembl', 'Gencode', 'NCBI'):
            print 'Error: source should be one of Ensembl/Gencode/NCBI'
            return None
        if source == 'NCBI':
            import NCBI_Genome
            self.GenomeAnno = NCBI_Genome.NCBI_Genome_Class(anno_file)
            self.bin = binize(self.GenomeAnno.gff3_container['RNA'], bw=bin_width)
            self.NC_To_Chr = self.GenomeAnno.build_NC_To_chr_dict()
            self.Chr_To_NC = self.GenomeAnno.build_chr_To_NC_dict()
        if source in ('Ensembl', 'Gencode'):
            import GENCODE_Genome
            self.GenomeAnno = GENCODE_Genome.GENCODE_Genome_Class(anno_file)
            self.bin = binize(self.GenomeAnno.gtf_container['RNA'], bw=bin_width)
        print "genome2TransClass: 输入是0-based，输出是1-based"
    def genomeCoor2TransCoor(self, Chr, Start, End, Strand, pureTransID=False, inChrSym='chr'):
        """函数返回这个基因组坐标对应的转录本的坐标. 
            Chr可以是'chr1'或者'1'
            如果source是NCBI, chrSym是一个有效的参数：chr/NC
        """
        if self.source == 'NCBI':
            if inChrSym not in ('chr', 'NC'):
                print 'Error: inChrSym should be one of chr/NC'
            ChrID = Chr
            if inChrSym == 'chr':
                ChrID = self.Chr_To_NC[Chr]
            raw_transCoorList = genomeRange2TransCoor(self.bin, self.GenomeAnno.gff3_container, ChrID, Start=Start, End=End, Strand=Strand, bw=self.bw)
            ['NC_000001.11', 17606, 17706, 'rna1', 487, 587]
            clean_transCoorList = []
            for (chrID, chr_start, chr_end, rnaID, trans_start, trans_end) in raw_transCoorList:
                RNA = self.GenomeAnno.gff3_container['RNA'][rnaID]
                "rnaID to TransID"
                try:
                    TransID = RNA['transcript_id']
                    if pureTransID: TransID = TransID.strip().split('.')[0]
                except KeyError:
                    if RNA['gbkey'] in ('tRNA', 'rRNA'):
                        TransID = RNA['gene']+'_'+RNA['gbkey']
                    else:
                        continue
                clean_transCoorList.append( [Chr, chr_start, chr_end, TransID, trans_start, trans_end] )
            return clean_transCoorList
        if self.source in ('Ensembl', 'Gencode'):
            raw_transCoorList = genomeRange2TransCoor(self.bin, self.GenomeAnno.gtf_container, Chr=Chr, Start=Start, End=End, Strand=Strand, bw=self.bw)
            clean_transCoorList = []
            for (Chr, chr_start, chr_end, rnaID, trans_start, trans_end) in raw_transCoorList:
                if pureTransID: rnaID = rnaID.strip().split('.')[0]
                clean_transCoorList.append( [Chr, chr_start, chr_end, rnaID, trans_start, trans_end] )
            return clean_transCoorList




