from nested_dict import nested_dict
import re, subprocess

def gene_type(raw_type):
        "Convert Raw Gene Type to Our Defined Gene Type"
        valid_gene_type = ('protein_coding', 'pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'misc_RNA', 'rRNA')
        lncRNA_class = ('antisense','lincRNA','processed_transcript','sense_intronic','TEC','sense_overlapping')
        if raw_type == 'protein_coding': return 'mRNA';
        if raw_type in valid_gene_type: return raw_type;
        if re.match('.*pseudogene',raw_type): return 'pseudogene';
        if raw_type in lncRNA_class: return 'lncRNA';
        return 'other'

def DG2Frame(dfFile, frameFile):
    """
    for paris.py output .dg
    Group 14 == position ENSDART00000055085(-):1736-1752|ENSDART00000172018(-):279-303, support 2, left 2, right 3098, score 0.500484.
    """
    # trans_dict = loadTransGtfBed2('/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/refseq_ensembl_merge.tarns.bed.2')
    # trans_dict = loadTransGtfBed2('/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/refseq_ensembl91_merge.tarns.bed.2')
    trans_dict = loadTransGtfBed2('/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/human_zebrafish/human_zebrafish.bed.2')
    import re
    pattern = "Group\s(\d+).*position\s([\w.\d]+)\(([+-])\):(\d+)-(\d+)\|([\w.\d]+)\(([+-])\):(\d+)-(\d+).*support\s(\d+).*left\s(\d+).*right\s(\d+).*score\s([\d\.]*\d)"
    DG = open(dfFile)
    FRAME = open(frameFile, 'w')
    print >>FRAME, "#Group\tlchr\tlstrand\tlstart\tlend\trchr\trstrand\trstart\trend\tsupport\tlcount\trcount\tscore\tltype\trtype"
    BEDPE = open(frameFile.replace('.txt', '.bedpe'), 'w')
    line = DG.readline()
    while line:
        if line.startswith('Group'):
            data = re.findall(pattern, line)
            print data
            if data[0][2] == '+' and data[0][6] == '+':
                ltype = trans_dict[data[0][1]]['type']
                rtype = trans_dict[data[0][5]]['type']
                print >>FRAME, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (data[0]) + '\t%s\t%s'%(ltype, rtype)
                Group, lchr, lstrand, lstart, lend, rchr, rstrand, rstart, rend, support, lcount, rcount, score = data[0]
                print >> BEDPE, '\t'.join([lchr, lstart, lend, rchr, rstart, rend, Group, support, lstrand, rstrand, lcount, rcount, score, ltype, rtype])
        line = DG.readline()
    FRAME.close()
    DG.close()
    BEDPE.close()

def DG2Frame2(dfFile, frameFile):
    """
    for dg2cluster output
    Group 0 == position ENSMUST00000000001.4;ENSMUSG00000000001.4;Gnai3(+):4-30|ENSMUST00000000001.4;ENSMUSG00000000001.4;Gnai3(+):62-73, support 2, left 3, right 2, score 0.583.
    """
    trans_dict = loadTransGtfBed2()
    import re
    pattern = "Group\s(\d+).*position\s([\w.\d;-]+)\(([+-])\):(\d+)-(\d+)\|([\w.\d;-]+)\(([+-])\):(\d+)-(\d+).*support\s(\d+).*left\s(\d+).*right\s(\d+).*score\s([\d\.]*\d)"
    DG = open(dfFile)
    FRAME = open(frameFile, 'w')
    print >>FRAME, "#Group\tlchr\tlstrand\tlstart\tlend\trchr\trstrand\trstart\trend\tsupport\tlcount\trcount\tscore\tltype\trtype"
    line = DG.readline()
    while line:
        if line.startswith('Group'):
            data = re.findall(pattern, line)
            #print data[0]
            if data[0][2] == '+' and data[0][6] == '+':
                ltype = gene_type(trans_dict[data[0][1].split(';')[0]]['type'])
                rtype = gene_type(trans_dict[data[0][5].split(';')[0]]['type'])
                data_ls = list(data[0])
                data_ls[1] = data_ls[1].split(';')[0]
                data_ls[5] = data_ls[5].split(';')[0]
                print >>FRAME, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (tuple(data_ls)) + '\t%s\t%s'%(ltype, rtype)
        line = DG.readline()
    FRAME.close()
    DG.close()

def loadTransGtfBed2(ref_bed='/Share/home/zhangqf7/gongjing/mes/ref/mm10.transCoor.bed.2'):
    H = open(ref_bed)
    line = H.readline()
    trans_dict = nested_dict()
    header_ls = ['tx', 'gene', 'type', 'length', 'utr_5_start', 'utr_5_end', 'cds_start', 'cds_end', 'utr_3_start', 'utr_3_end']
    while line:
        if line.startswith('#'): line = H.readline(); continue
        arr = line.strip().split('\t')
        gene = arr[1].split('=')[0].split()[0]
        for i,j in zip(header_ls, arr):
            trans_dict[arr[0]][i] = j
        line = H.readline()
    H.close()
    print "read: %s, n=%s"%(ref_bed, len(trans_dict))
    return trans_dict.to_dict()

def bedpe_to_bam(bedpe=None):
    if bedpe is None:
        bedpe = '/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-2-allclean/27-DG.bedpe'
    genome_size = '/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/combine_transcriptome_star/chrNameLength.txt'
    bam = bedpe.replace('.bedpe', '.bam')
    sort_bam = bam.replace('.bam', '.sort.bam')
    subprocess.call(["bedtools bedpetobam -i %s -g %s > %s"%(bedpe, genome_size, bam)],shell=True)
    subprocess.call(["samtools sort -o %s %s"%(sort_bam, bam)], shell=True)
    subprocess.call(["samtools index %s"%(sort_bam)], shell=True)

def main():
    #dfFile = '/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-5-rep-combine/27-DG'
    #dfFile = '/Share/home/zhangqf7/gongjing/zebrafish/data/paris/PARIS_2016/mES/8.duplexgroup'
    # dg_ls = ['/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-%s-rep-combine/downsampling2/27-DG'%(i) for i in [1,2,3,4,5]]
    # for dfFile in dg_ls:
    #     frameFile = dfFile + '.txt'
    #     DG2Frame(dfFile=dfFile, frameFile=frameFile)

    #bedpe_to_bam('/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-5-allclean/27-DG.bedpe')

    dfFile = '/Share/home/zhangqf7/gongjing/zebrafish/data/paris/shi-zp-23/27-DG.trimaticAdapter'
    frameFile = dfFile + '.txt'
    DG2Frame(dfFile=dfFile, frameFile=frameFile)

if __name__ == '__main__':
    main()

