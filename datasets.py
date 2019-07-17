import pandas as pd

##dynamic regions
dsr_egg_1cell = 'data/dynamic_region/egg_cell1/window-anno.bed'
dsr_1cell_4cell = 'data/dynamic_region/cell1_cell4/window-anno.bed'
dsr_4cell_64cell = 'data/dynamic_region/cell4_cell64/window-anno.bed'
dsr_64cell_sphere = 'data/dynamic_region/cell64_sphere/window-anno.bed'
dsr_sphere_shield = 'data/dynamic_region/sphere_shield/window-anno.bed'
way2345 = 'data/dynamic_region/way2345/way2345.bed'
way1 = 'data/dynamic_region/way1/way1.bed'

def readdsr(dsr_f):
    df = pd.read_csv(dsr_f, sep="\t", header=None, usecols=[0,1,2,7])
    df.columns = ['TransID', 'start(0-based)', 'end(0-based)', 'region']
    return df

def readhotregion(f):
    df = pd.read_csv(f, sep="\t", header=None)
    df.columns = ['TransID', 'start(0-based)', 'end(0-based)']
    return df

def readAlldsr(Exwriter):
    file_ls = [dsr_egg_1cell, dsr_1cell_4cell, dsr_4cell_64cell, dsr_64cell_sphere, dsr_sphere_shield]
    sample_ls = ['0hpf->0.4hpf', '0.4hpf->1hpf', '1hpf->2hpf', '2hpf->4hpf', '4hpf->6hpf']
    for i, j in zip(file_ls, sample_ls):
        df = readdsr(i)
        df.to_excel(Exwriter, sheet_name=j, index=False)
    return Exwriter

def main():
    writer = pd.ExcelWriter('data/dynamic_region/all_dsr_and_hot_region.xlsx')
    writer = readAlldsr(writer)
    readhotregion(way2345).to_excel(writer, sheet_name='hot_dsr_way2345', index=False)
    readhotregion(way1).to_excel(writer, sheet_name='specific_dsr_way1', index=False)
    writer.close()

main()