import re
# import sys,glob ## command line stuff
import numpy as np
import sys
import related_funs
import gamma
import nrrd


class class_analysis_dos_nrrd:
    def __init__(self, path2doslist):
        self.path2doselist=path2doslist
        self.path2log='./'

    def fun_accumulatedose(self):
        if len(self.path2doselist) < 3:
            errormess = 'Detects wrong input:"dose list at least three"'
            related_funs.writelog(self.path2log, errormess)
            sys.exit()
        datafile1, headfile1 = nrrd.read(self.path2doselist[1])
        for jj in range(1,len(self.path2doselist)-1):
            datafile, headfile = nrrd.read(self.path2doselist[jj])
            datafile1 += datafile
        # nrrd.write(self.path2doselist[-1]+'.nhdr', datafile1, headfile1)
        nrrd.write(self.path2doselist[-1]+'.nhdr', datafile1, headfile1,detached_header=True)
        print(self.path2doselist[-1])
        print("finished writing of: ",self.path2doselist[-1][self.path2doselist[-1].rfind('/')+1:])

    def fun_modifynrrdhed(self,path2dosenrrd):
        findwrongspace=False
        with open(path2dosenrrd,'rb') as f:
            savelines =''
            for lines in f.readlines():
                if ('space directions' in lines):
                    if (', ' in lines or ' ,' in lines):
                        findwrongspace=True
                        mod_line=lines.replace(', ',',')
                        mod_line=mod_line.replace(' ,',',')
                else:
                    mod_line=lines
                savelines+=mod_line+'\n'
        if findwrongspace:
            with open(path2dosenrrd,'wb') as f:
                f.writelines(savelines)


