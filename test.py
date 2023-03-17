
with open('1_Lung_P_10_11phys.nrrd', 'r+') as f:
    savelines = ''
    for lines in f.readlines():
        savelines += lines

with open('test.nrrd', 'w+') as f:
    f.writelines(savelines)


