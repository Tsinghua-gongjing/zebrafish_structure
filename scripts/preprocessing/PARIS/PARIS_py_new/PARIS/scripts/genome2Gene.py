#-*- coding:utf-8 -*-


"""
给出一段基因组的区域，返回这个区域中基因对应的坐标

  使用方法：  
    1. ensembl/gencode转换
        gencode_gg = genome2GeneClass("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Gencode/GRCh37.gtf")
        gencode_gg.genomeCoor2GeneCoor('chr1', 332235, 332236, '-', pureGeneID=False)
    2. ncbi hg19转换
        ncbi_gg = genome2GeneClass("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/ncbi/GRCh37.p13.gff3", source='NCBI')
        ncbi_gg.genomeCoor2GeneCoor('chr1', 92222, 133748, '-', pureGeneID=False)
    3. ensembl/gencode mm10
        gencode_gg = genome2GeneClass("/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.GRCm38.81.gtf")
        gencode_gg.genomeCoor2GeneCoor('chr1', 5276106-2, 5276106, '-', pureGeneID=False)

  两个常用的GTF文件：
    /Share/home/zhangqf8/lipan/DYNAMIC/GTF/Homo_sapiens.gtf
    /Share/home/zhangqf8/lipan/DYNAMIC/GTF/Mus_musculus.GRCm38.81.gtf

输入：
  注意这里的输入是左闭右开，0-based坐标； 

返回：
  包含基因区段的列表，每一个元素格式是这样的：[chrID, genomeStart, genomeEnd, geneID, geneStart, geneEnd]
  其中返回的基因的起点和终点都是1-based的，且前闭后闭; 基因组的起点和终点也是1-based的，前闭后闭

异常：
  当染色体不在GTF文件中时，发出一个KeyError的异常，注意在使用时进行捕获

Version:
    2017-6-14 修改genome2Trans得到genome2Gene
"""


import re, sys, os, getopt, time, datetime

def now():
    "返回当前时间字符串"
    return str(datetime.datetime.now())


def binize(Gene_Container, bw=100000):
    "把基因组上的区域分段，加速查找过程。bw越小，占内存越高，越快"
    def statisticPsuedoChrSize(Gene_Container):
        """获取染色体的最大计数( <= 染色体的长度)
            Gene_Container{rnaID1 => {chr, start, end, strand...}...}
        """
        ChrLen = {}
        for geneID in Gene_Container:
            Chr = Gene_Container[geneID]['chr']
            End = int(Gene_Container[geneID]['end'])
            try:
                ChrLen[Chr] = End if End > ChrLen[Chr] else ChrLen[Chr]
            except KeyError:
                ChrLen[Chr] = End
        return ChrLen
    print "Binize genome to speed up searching.\n\t%s" % (now(), )
    chr_size = statisticPsuedoChrSize(Gene_Container)
    Bin = {}
    for Chr in chr_size.keys():
        Bin[Chr] = {}; Bin[Chr]['+'] = {}; Bin[Chr]['-'] = {}
        count = int( chr_size[Chr]/bw ) + 1
        for idx in range(count):
            Bin[Chr]['+'][idx] = []; Bin[Chr]['-'][idx] = []
    for geneID in Gene_Container:
        Chr = Gene_Container[geneID]['chr']
        Strand = Gene_Container[geneID]['strand']
        Start = int( int(Gene_Container[geneID]['start']) / bw )
        End = int( int(Gene_Container[geneID]['end']) / bw )
        for idx in range(Start, End+1):
            Bin[Chr][Strand][idx].append( geneID )
    return Bin

def genomeRange2GeneCoor(ref_bin, GFF_Container, Chr, Start, End, Strand, bw=100000):
    "左闭右开，0-based"
    idxBin = int(Start/bw)
    overlapRegions = []
    checkedGenes = []
    while idxBin <= int(End/bw):
        if idxBin > ref_bin[Chr][Strand].keys()[-1]:
            break
        for geneID in ref_bin[Chr][Strand][idxBin]:
            if geneID not in checkedGenes:
                checkedGenes.append(geneID)
            else:
                continue
            if End < int(GFF_Container['gene'][geneID]['start']) or Start > int(GFF_Container['gene'][geneID]['end']):
                continue
            #for transcriptID in GFF_Container['gene_info'][geneID]['transcript']:
            overlapRegion = overlapGene(GFF_Container, geneID, Strand, Start, End)
            #for item in overlapRegion: item.insert(0, Chr)
            if overlapRegion != None:
                overlapRegion.insert(0, Chr)
                overlapRegions += [overlapRegion]
        idxBin += 1
    return overlapRegions

def overlapGene(GFF_Container, geneID, strand, absStart, absEnd):
    "查找一个基因组区段与基因的交集"
    absStart += 1
    gene = GFF_Container['gene'][geneID]
    geneStart = int(gene['start'])
    geneEnd = int(gene['end'])
    geneLength = geneEnd - geneStart + 1
    if absStart <= geneEnd and geneStart <= absEnd:
        if strand == '+':
            startPosInGene = absStart - geneStart + 1
            if startPosInGene < 1: startPosInGene = 1
            endPosInGene = absEnd - geneStart + 1
            if endPosInGene > geneLength: endPosInGene = geneLength
        elif strand == '-':
            startPosInGene = geneEnd - absEnd + 1
            if startPosInGene < 1: startPosInGene = 1
            endPosInGene = geneEnd - absStart + 1
            if endPosInGene > geneLength: endPosInGene = geneLength
        absStartInExon = absStart if absStart >= geneStart else geneStart
        absEndInExon = geneEnd if absEnd >= geneEnd else absEnd
        return [absStartInExon, absEndInExon, geneID, startPosInGene, endPosInGene]
    return None


"""测试函数
import NCBI_Genome
GFF_File = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/ncbi/GRCh38.p7.gff3'
NCBI_Genome_Anno = NCBI_Genome.NCBI_Genome_Class(GFF_File)
NCBI_GFF_Container = NCBI_Genome_Anno.gff3_container
Bin = binize(NCBI_GFF_Container['gene'], bw=100000)
genomeRange2GeneCoor(Bin, NCBI_GFF_Container, Chr='NC_000001.11', Start=34609, End=36081, Strand='-', bw=100000)


import GENCODE_Genome
GTF_File = '/Share/home/zhangqf8/lipan/DYNAMIC/GTF/Gencode/GRCh38.gtf'
GENCODE_Genome_Anno = GENCODE_Genome.GENCODE_Genome_Class(GTF_File)
GENCODE_GFF_Container = GENCODE_Genome_Anno.gtf_container
Bin = binize(GENCODE_GFF_Container['gene'], bw=100000)
genomeRange2GeneCoor(Bin, GENCODE_GFF_Container, Chr='chr1', Start=14402, End=29570, Strand='-', bw=100000)
"""

class genome2GeneClass(object):
    """把基因组的坐标转成基因的坐标的类"""
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
            self.bin = binize(self.GenomeAnno.gff3_container['gene'], bw=bin_width)
            self.NC_To_Chr = self.GenomeAnno.build_NC_To_chr_dict()
            self.Chr_To_NC = self.GenomeAnno.build_chr_To_NC_dict()
        if source in ('Ensembl', 'Gencode'):
            import GENCODE_Genome
            self.GenomeAnno = GENCODE_Genome.GENCODE_Genome_Class(anno_file)
            self.bin = binize(self.GenomeAnno.gtf_container['gene'], bw=bin_width)
        print "genome2GeneClass: 输入是0-based，输出是1-based"
    def genomeCoor2GeneCoor(self, Chr, Start, End, Strand, pureGeneID=False, inChrSym='chr'):
        """函数返回这个基因组坐标对应的基因的坐标. 
            Chr可以是'chr1'或者'1'
            如果source是NCBI, chrSym是一个有效的参数：chr/NC
        """
        if self.source == 'NCBI':
            if inChrSym not in ('chr', 'NC'):
                print 'Error: inChrSym should be one of chr/NC'
            ChrID = Chr
            if inChrSym == 'chr':
                ChrID = self.Chr_To_NC[Chr]
            raw_geneCoorList = genomeRange2GeneCoor(self.bin, self.GenomeAnno.gff3_container, ChrID, Start=Start, End=End, Strand=Strand, bw=self.bw)
            ['NC_000001.11', 17606, 17706, 'geneID', 487, 587]
            clean_geneCoorList = []
            for (chrID, chr_start, chr_end, geneID, gene_start, gene_end) in raw_geneCoorList:
                gene = self.GenomeAnno.gff3_container['gene'][geneID]
                "geneID to geneID_geneName"
                try:
                    geneName = self.GenomeAnno.gff3_container['gene'][geneID]['gene']
                except KeyError:
                    print 'Error: Unexpect GeneID %s with no gene attribute, Skip' % (geneID, )
                    continue
                clean_geneCoorList.append( [Chr, chr_start, chr_end, geneID+'_'+geneName, gene_start, gene_end] )
            return clean_geneCoorList
        if self.source in ('Ensembl', 'Gencode'):
            raw_geneCoorList = genomeRange2GeneCoor(self.bin, self.GenomeAnno.gtf_container, Chr=Chr, Start=Start, End=End, Strand=Strand, bw=self.bw)
            clean_geneCoorList = []
            for (Chr, chr_start, chr_end, geneID, gene_start, gene_end) in raw_geneCoorList:
                try:
                    geneName = self.GenomeAnno.gtf_container['gene'][geneID]['gene_name']
                except KeyError:
                    print 'Error: Unexpect GeneID %s with no gene_name attribute, Skip' % (geneID, )
                    continue
                if pureGeneID: geneID = geneID.strip().split('.')[0]
                clean_geneCoorList.append( [Chr, chr_start, chr_end, geneID+'_'+geneName, gene_start, gene_end] )
            return clean_geneCoorList
















