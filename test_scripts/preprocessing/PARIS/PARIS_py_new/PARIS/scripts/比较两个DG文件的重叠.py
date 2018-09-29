
def readDG(DGFile):
    import re
    DGList = []
    IN = open(DGFile)
    line = IN.readline()
    while line:
        if line.startswith('Group'):
            data = line.strip().split()
            split_items = re.split("[():|]", data[4].strip(','))
            tmp_idd = split_items[:2] + [ int(it) for it in split_items[3].split('-') ] + split_items[4:6] + [ int(it) for it in split_items[7].split('-') ]
            DGList.append( tmp_idd )
        line = IN.readline()
    DGList.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    return DGList

def TwoDGOverLaps(DG1, DG2, similarity=2):
    commonDG = 0
    for dg_item1 in DG1:
        for dg_item2 in DG2:
            if dg_item1[0] == dg_item2[0] and dg_item1[1] == dg_item2[1] and dg_item1[4] == dg_item2[4] and dg_item1[5] == dg_item2[5]:
                if dg_item1[2] < dg_item2[3] and dg_item2[2] < dg_item1[3]:
                    # left arm have overlap
                    if dg_item1[6] < dg_item2[7] and dg_item2[6] < dg_item1[7]:
                        # right arm have overlap
                        left_overlap = abs(dg_item2[2]-dg_item1[2]) + abs(dg_item2[3]-dg_item1[3])
                        right_overlap = abs(dg_item1[6]-dg_item2[6]) + abs(dg_item2[7]-dg_item1[7])
                        if left_overlap + right_overlap <= similarity:
                            commonDG += 1
                            break
    print commonDG

python_dg = readDG("DG_python")
perl_dg = readDG("DG_perl")

TwoDGOverLaps(python_dg, perl_dg, similarity=2)


