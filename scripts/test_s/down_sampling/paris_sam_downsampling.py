import os, sys
import random

ChimericSam = sys.argv[2]
Chinericjunction = sys.argv[3]
AlignedSam = sys.argv[1]
outputDir = sys.argv[4]
if not os.path.exists(outputDir):
    os.makedirs(outputDir)
p = float(sys.argv[5])
count = sys.argv[6] # 
outChimericSam = os.path.join(outputDir, os.path.basename(ChimericSam).replace("sam", "ds.sam"))
outChimericJunction = os.path.join(outputDir, os.path.basename(Chinericjunction).replace("junction", "ds.junction"))
outAlignedSam = os.path.join(outputDir, os.path.basename(AlignedSam).replace("sam", "ds.sam"))

random.seed(1234)

random_seed_ls = random.sample(range(1, 10000), 100)

random_seed = random_seed_ls[int(count)-1]
random.seed(random_seed)

print "random seed: %s"%(random_seed)

with open(ChimericSam, 'r') as c_sam, open(Chinericjunction, 'r') as c_jun, open(outChimericSam, 'w') as out_c_sam, open(outChimericJunction, 'w') as out_c_jun:
    l_c_sam = c_sam.readline()
    l_c_jun = c_jun.readline()
    while (l_c_jun and l_c_sam):
        while l_c_sam:
            if l_c_sam.startswith("@"):
                out_c_sam.write(l_c_sam)
                l_c_sam = c_sam.readline()
            else:
                break
        if random.random() < p:
            out_c_jun.write(l_c_jun)
            #junction 1 line : sam 2 line
            out_c_sam.write(l_c_sam)
            l_c_sam = c_sam.readline()
            out_c_sam.write(l_c_sam)
            l_c_jun = c_jun.readline()
            l_c_sam = c_sam.readline()
        else:
            l_c_jun = c_jun.readline()
            l_c_sam = c_sam.readline()
            l_c_sam = c_sam.readline()

with open(AlignedSam, 'r') as a_sam, open(outAlignedSam, 'w') as out_a_sam:
    l_a_sam = a_sam.readline()
    while l_a_sam:
        while l_a_sam:
            if l_a_sam.startswith("@"):
                out_a_sam.write(l_a_sam)
                l_a_sam = a_sam.readline()
            else:
                break
        if random.random() < p:
            out_a_sam.write(l_a_sam)
            l_a_sam = a_sam.readline()
        else:
            l_a_sam = a_sam.readline()
