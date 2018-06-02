import sys
import re
import structure

def readIc(ic_file):
	"Read icSHAPE Reactivity File"
	icDict = dict()
	IN = open(ic_file)
	line = IN.readline()
	while line:
		arr = line.strip().split()
		icDict[arr[0]] = arr[3:]
		line = IN.readline()
	IN.close()
	return icDict

def read_fasta(fasta):
	seq_dict = {}
	with open(fasta, 'r') as FASTA:
		for line in FASTA:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			if line.startswith('>'):
				seq_id = line.split(' ')[0].replace('>', '')
				seq_dict[seq_id] = ""
			else:
				seq_dict[seq_id] += line.strip()
	return seq_dict

def main(icshape_out=None, fasta=None, tx_txt=None, savedir=None):
	"""
	icshape_out: icshape reactivity profile
	transcript_id
	start, end: 0-based
	fasta: reference 
	"""
	if icshape_out is None:
		icshape_out = '/Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/shield.icshape.w200.s30.T2.t200.out'
	if fasta is None:
		fasta = '/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.fa'
	if tx_txt is None:
		tx_txt = '/Share/home/zhangqf7/gongjing/zebrafish/result/candidate_structure/decay/candidate_decay.txt'
	if savedir is None:
		savedir = '/Share/home/zhangqf7/gongjing/zebrafish/result/candidate_structure/decay'
	print
	print "icshape_out: %s"%(icshape_out)
	"""
	print "transcript_id: %s"%(transcript_id)
	print "start: %s"%(start)
	print "end: %s"%(end)
	"""
	print "fasta: %s"%(fasta)
	print "tx_txt: %s"%(tx_txt)
	print

	icDict = readIc(icshape_out)
	seq_dict = read_fasta(fasta)

	with open(tx_txt, 'r') as TXT:
		for line in TXT:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			arr = line.split('\t')
			transcript_id = arr[0]
			start = int(arr[1])-20
			end = int(arr[2])+20
			seq = seq_dict[transcript_id][int(start)-1:int(end)]
			shape_reactivity = icDict[transcript_id][int(start)-1:int(end)]
			structure_list = structure.predictStructure(Seq=seq, Shape=shape_reactivity, bp_constraint=[], mfe=True, clean=True, si=-0.6, sm=1.8, md=0)
			print "predict structure:", ''.join(seq), structure_list

			sample = icshape_out.split('/')[-1].split('.')[0]
			savefn = '%s/%s.%s.%s.%s.txt'%(savedir, transcript_id, start, end, sample)
			# cross_link_points = structure.search_TT_cross_linking(sequence=seq, ct_list=structure_list)
			highlight_start = 21
			highlight_end = highlight_start + end - start - 40
			with open(savefn, 'w') as SAVEFN:
				print >>SAVEFN, '>%s.%s.%s.%s|%s.%s'%(transcript_id, start, end, sample, highlight_start, highlight_end)
				print >>SAVEFN, seq
				print >>SAVEFN, structure_list

	"""
	seq = seq_dict[transcript_id][int(start):int(end)]
	shape_reactivity = icDict[transcript_id][int(start):int(end)]

	structure_list = structure.predictStructure(Seq=seq, Shape=shape_reactivity, bp_constraint=[], mfe=True, clean=True, si=-0.6, sm=1.8, md=0)
	print "predict structure:", ''.join(seq), structure_list
	"""


if __name__ == '__main__':
	# main(icshape_out=sys.argv[1], transcript_id=sys.argv[2], start=sys.argv[3], end=sys.argv[4], fasta=sys.argv[5])
	main()

	"""
	python structure_predict.py /Share/home/zhangqf7/gongjing/zebrafish/data/icSHAPE_final_out_new_win/1cell.icshape.w200.s30.T2.t200.out NM_001123007 100 130 
	/Share/home/zhangqf7/gongjing/zebrafish/data/reference/transcriptome/danRer10.refSeq.transcriptome.fa
	"""