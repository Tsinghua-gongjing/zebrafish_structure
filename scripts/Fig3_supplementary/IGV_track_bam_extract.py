import subprocess
import sys, os

def bam_info():
	bam_dict = {'RNAseq_h2':'/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/control_h2-rep2.sorted.bam',
				'RNAseq_h4':'/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/control_h4-rep2.sorted.bam',
				'RNAseq_h6':'/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/control_h6-rep2.sorted.bam',
				'RIP_h2':'/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/hur-h2.5_fragment.sorted.bam',
				'RIP_h4':'/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/hur-h4_fragment.sorted.bam',
				'RIP_h6':'/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/hur-h6_fragment.sorted.bam',
				'RNAseq_MO_h2':'/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/hur_MO-h2_rep2.sorted.bam',
				'RNAseq_MO_h4':'/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/hur_MO-h4_rep2.sorted.bam',
				'RNAseq_MO_h6':'/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/hur_MO-h6_rep2.sorted.bam'}
	return bam_dict

def transcript_ls():
	decay_ls = ['NM_001110127', 'NM_207055', 'NM_212563', 'NM_001013351', 'NM_001001943',
	            'NM_001034978', 'NM_001115057', 'NM_001033746', 'NM_001077537', 'NM_201291']
	stable_ls = ['NM_130941', 'NM_212792', 'NM_001082921', 'NM_200276']

	return ' '.join(decay_ls), ' '.join(stable_ls)

def fetch_transcript_bam():
	if len(sys.argv) == 1:
		transcript = 'NM_130941'
		savedir = '/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/candidate_bam'
	else:
		transcript = sys.argv[1]
		savedir = sys.argv[2]
	bam_dict = bam_info()
	decay, stable = transcript_ls()
	sample_ls = ['RNAseq_h4', 'RNAseq_h6', 'RNAseq_MO_h4', 'RNAseq_MO_h6', 'RIP_h4', 'RIP_h6']
	savefn = 'stable'
	fa = '/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/refseq_ensembl91_merge.tarns.fa'
	for sample in sample_ls:
		if os.path.isdir('%s/%s'%(savedir, savefn)):
			pass
		else:
			os.makedirs('%s/%s'%(savedir, savefn))
		save_bam = '%s/%s/%s.bam'%(savedir, savefn, sample)
		subprocess.call(["samtools view -bS %s %s >%s; samtools index %s; igvtools count -w 1 %s %s %s"%(bam_dict[sample], stable, save_bam, save_bam, save_bam, save_bam.replace('.bam', '.wig'), fa)], shell=True)


def wig_normalize_by_depth(wig=None, mapped_reads_num=None, reads_len=150, effective_genome_size=37421801):
	depth = mapped_reads_num * reads_len / float(effective_genome_size)
	wig_norm = wig.replace('.wig', '.norm.wig')
	WIG_NORM = open(wig_norm, 'w')
	with open(wig, 'r') as WIG:
		for line in WIG:
			line = line.strip()
			if not line or line.startswith(('#', 'track', 'variableStep')):
				print >>WIG_NORM, line.strip()
				continue
			arr = line.split('\t')
			value_norm = float(arr[1]) / depth
			print >>WIG_NORM, '\t'.join(map(str, [arr[0], value_norm]))
	WIG_NORM.close()

	bed_norm = wig_norm.replace('.wig', '.bed')
	tranWig_to_genomeWig_py = '/Share/home/zhangqf7/gongjing/zebrafish/script/coor/GAP.1.1.0/tranWig_to_genomeWig.py'
	gene_coordinate_bed = '/Share/home/zhangqf7/gongjing/zebrafish/data/reference/gtf/danRerV10.refGene.genes.bed'
	subprocess.call(["wig2bed < %s > %s; python %s %s %s"%(wig_norm, bed_norm, tranWig_to_genomeWig_py, gene_coordinate_bed, bed_norm)], shell=True)

if __name__ == '__main__':
	fetch_transcript_bam()

	wig_normalize_by_depth(wig='/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/candidate_bam/decay/RIP_h4.wig', mapped_reads_num=36790087)
	wig_normalize_by_depth(wig='/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/candidate_bam/decay/RIP_h6.wig', mapped_reads_num=24055217)
	wig_normalize_by_depth(wig='/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/candidate_bam/decay/RNAseq_h4.wig', mapped_reads_num=34548669)
	wig_normalize_by_depth(wig='/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/candidate_bam/decay/RNAseq_h6.wig', mapped_reads_num=41536928)
	wig_normalize_by_depth(wig='/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/candidate_bam/decay/RNAseq_MO_h4.wig', mapped_reads_num=40103529)
	wig_normalize_by_depth(wig='/Share/home/zhangqf7/gongjing/zebrafish/data/RIP/candidate_bam/decay/RNAseq_MO_h6.wig', mapped_reads_num=44087896)