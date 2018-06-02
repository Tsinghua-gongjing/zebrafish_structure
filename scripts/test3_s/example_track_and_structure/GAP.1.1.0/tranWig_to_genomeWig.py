from ParseTrans import *
import sys
from nested_dict import nested_dict

def bed_coordinate_conversion(genome_coordinate_bed=None, trans_bed=None, genome_bed=None, write_wig=1):
	if genome_bed is None:
		genome_bed = trans_bed.replace('.bed', '.genome.bed')
	
	print "ref genome bed: %s"%(genome_coordinate_bed)
	print "input trans bed: %s"%(trans_bed)
	print "output genome bed: %s"%(genome_bed)

	Parser = ParseTransClass(genomeCoorBedFile = genome_coordinate_bed)
	GENOME_BED = open(genome_bed, 'w')
	convert_dict = nested_dict(2, list) # chr:['start'],['end'],['score']
	with open(trans_bed, 'r') as TRANS_BED:
		for line in TRANS_BED:
			line = line.strip()
			if not line or line.startswith(('#', 'track')):
				continue
			arr = line.split('\t')
			trans_id = arr[0]
			start = int(arr[1])+1 # need 1-based
			end = int(arr[2])
			score = float(arr[4])
			# print trans_id,start
			convert_ls = Parser.transCoor2geneCoor(trans_id, start, end) # output also 1-based
			"""
			>>> Parser.transCoor2geneCoor("ENSDART00000168451", 1, 1)
			[['1', 18716, 18716, 'ENSDARG00000102097', 1, 1]]
			"""
			"""
			for i in convert_ls:
				arr[0] = i[0]
				arr[1] = i[1]
				arr[2] = i[2]
				arr[3] = i[3]
				print >>GENOME_BED, '\t'.join(map(str, arr))
			"""
			arr[0] = convert_ls[0][0]
			arr[1] = convert_ls[0][1]-1
			arr[2] = convert_ls[0][2]
			arr[3] = ','.join([i[3] for i in convert_ls])
			print >>GENOME_BED, '\t'.join(map(str, arr))
			convert_dict[arr[0]]['start'].append(arr[1])
			convert_dict[arr[0]]['end'].append(arr[2])
			convert_dict[arr[0]]['score'].append(arr[4])
	GENOME_BED.close()

	if write_wig:
		wig = trans_bed.replace('.bed', '.genome.wig')
		with open(wig, 'w') as WIG:
			print >>WIG, 'track type=wiggle_0'
			for i,j in convert_dict.items():
				print >>WIG, 'variableStep chrom=%s span=1'%(i)
				start_ls, score_ls = zip(*sorted(zip(j['start'], j['score'])))
				for start,score in zip(start_ls, score_ls):
					print >>WIG, '\t'.join(map(str, [start+1, score]))

if __name__ == '__main__':
	bed_coordinate_conversion(genome_coordinate_bed=sys.argv[1], trans_bed=sys.argv[2])