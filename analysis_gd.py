#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 13:58:34 2022

@author: ysheng
"""
import re
# import sys,glob ## command line stuff
import numpy as np
import sys
import related_funs
import gamma


class class_analysis_gd:
    def __init__(self, patientID, planname, OptMethod,targetnamelist, targetdoselist, oarnamelist, externalname, fractions,
                 savepath, gammaEva, robustevaluation, path2gdlist, nameofgdlist,doseshowninplansgd,Showall,randomrobust,
                 fractionsacc,referecegd,referencefraction,referenceplanneddose):
        self.patientID = patientID
        self.planname = planname
        self.OptMethod= OptMethod
        self.PlanDose=doseshowninplansgd
        self.targetnamelist = targetnamelist

        self.targetdoselist = targetdoselist
        self.oarnamelist = oarnamelist
        self.externalname = externalname
        self.lowerdoseforext = 999  # take lower
        # prescribed dose for calculate of V95 for Ext
        self.ExtV95 = 0
        self.fractions = fractions # fractions to calculate doses
        self.savepath = savepath
        self.gammaEva = gammaEva
        self.robustevaluation = robustevaluation
        self.robusteva = True
        self.Showall = Showall  # true, show result of all.
        if self.robustevaluation == '21':
            self.robust_suffix = ['_nom', '_nx', '_ny', '_nz', '_px', '_py', '_pz',
                                  '_hd', '_hd_nx', '_hd_ny', '_hd_nz', '_hd_px', '_hd_py', '_hd_pz',
                                  '_ld', '_ld_nx', '_ld_ny', '_ld_nz', '_ld_px', '_ld_py', '_ld_pz']
        elif self.robustevaluation == '9':
            self.robust_suffix = ['_nom', '_hd', '_ld', '_nx', '_ny', '_nz', '_px', '_py', '_pz']
        else:
            self.robusteva=False
            self.robust_suffix = ['']

        self.randomrobust=randomrobust
        self.fractionsacc=fractionsacc

        self.senarioname_suffix =[]
        if self.randomrobust is not None and self.fractionsacc is None: #either senario or fxacc. they should be evaluate separately
            calculatedsenarios=self.randomrobust.split(',')
            for jj in range(int(calculatedsenarios[0]),int(calculatedsenarios[1])+1):
                self.senarioname_suffix.append(str(jj))
        elif self.randomrobust is not None and self.fractionsacc is not None:
            calculatedsenarios = self.randomrobust.split(',')
            for jj in range(int(calculatedsenarios[0]),int(calculatedsenarios[1])+1):
                for kk in range(1,int(self.fractionsacc)+1):
                    surffixtemp=str(jj)+'_fxDVH'+'1-'+str(kk)
                    self.senarioname_suffix.append(surffixtemp)
        else:
            self.senarioname_suffix=['']

        self.filesuffixtoattach = []
        if self.robusteva:
            self.filesuffixtoattach = self.robust_suffix
        elif self.randomrobust is not None:
            self.filesuffixtoattach = self.senarioname_suffix

        self.ionType = 'Carbon' # when use this script, please check if your plan name contains P or C that indicates ion type!!!!
        if 'p' in self.planname or 'P' in self.planname:
            self.ionType='Proton'

        self.FileList = path2gdlist
        self.nameofgdlist = nameofgdlist

        if self.savepath == 'None':
            self.savepath = './'
        self.Vxx = [90, 95, 100, 105]
        self.Dxx = [2, 95, 98]
        self.Dcc = [1]
        self.path2log = self.savepath + '00_Doseana_processing.log'
        # self.path2log = '/home/yurii/Sheng/patient_data/00_Doseana_processing.log'

        self.voilist = []
        self.patientIDToW = []
        self.plannameToW = []
        self.VOI_names = []
        self.VOI_volumes = []
        self.VOI_pres_Dose = []
        self.VOI_Parameter = []
        self.VOI_data = []  # 1 dimension/9D/21D
        self.referenceDATAforCompare = []
        self.writelinesinfo = []

        self.referencegd = referecegd  # reference gd file path
        self.referenceVoidata=[]
        self.referencefraction=referencefraction
        self.referenceplanneddose=referenceplanneddose

    def fun_analysis_gd(self):
        if len(self.nameofgdlist) != len(self.FileList):
            errormess = 'Detects wrong input:"gdlist path and given described names are not the same, please check"'
            related_funs.writelog(self.path2log, errormess)
            sys.exit()
        related_funs.writelog(self.path2log, 'Start a new analysis')

        namesuffix=''
        if self.fractionsacc is not None:
            namesuffix='_fxacc'
        elif self.randomrobust is not None:
            namesuffix = '_senario'
        elif self.robusteva:
            namesuffix = '_worst'

        savedata_fildname = self.savepath + self.patientID + '_' + self.planname + '_' + \
                            "_".join(m for m in self.nameofgdlist) + namesuffix+'.txt'
        # write log
        writeloginfo = 'running patient: ' + self.patientID + ' plan: ' + self.planname
        related_funs.writelog(self.path2log, writeloginfo)
        referencedata = True
        if self.referencegd is not None:
            self.referenceVoidata = self.AnalyzeDVHvoirefdata(self.referencegd, '3Dassigned_ref')
        for fileNo in range(0, len(self.FileList)):
            fileToanalysis = self.FileList[fileNo]
            Definednameofdata = self.nameofgdlist[fileNo]
            # do not change order VOI_data first. so self.reference could be available.
            VOI_data = self.AnalyzeDVHvoidata(fileToanalysis, Definednameofdata, referencedata)
            if referencedata: # write some pre information and analysis reference data.
                # patientIDToW,plannameToW,VOI_names,VOI_volumes,VOI_pres_Dose,VOI_Parameter,VOI_data = \
                #     self.AnalyzeDVHReference(fileToanalysis, Definednameofdata)
                patientIDToW, plannameToW, VOI_names, VOI_volumes, VOI_pres_Dose, VOI_TargetdoseTec,\
                    VOI_Parameter,VOI_OptMethod,VOI_ionType = self.WriteDVHPreInfo(fileToanalysis)
                self.writelinesinfo.append(patientIDToW)
                self.writelinesinfo.append(plannameToW)
                self.writelinesinfo.append(VOI_names)
                self.writelinesinfo.append(VOI_volumes)
                self.writelinesinfo.append(VOI_pres_Dose)
                # add SIBH,SIBL,Normal
                self.writelinesinfo.append(VOI_TargetdoseTec)
                self.writelinesinfo.append(VOI_Parameter)
                self.writelinesinfo.append(VOI_OptMethod)
                self.writelinesinfo.append(VOI_ionType)
                if (len(self.referenceVoidata)!=0):
                    if isinstance(self.referenceVoidata[0], list):
                        [self.writelinesinfo.append(i) for i in self.referenceVoidata]
                    else:
                        self.writelinesinfo.append(self.referenceVoidata)
                referencedata = False

            if isinstance(VOI_data[0],list):
                [self.writelinesinfo.append(i) for i in VOI_data]
            else:
                self.writelinesinfo.append(VOI_data)

            # write gamma data to line
            # gammacri = []
            # gammaresult = []
            # writegammadata=''
            # gammacri=['3/3','2/2']
            # if self.gammaEva:
            #     print('detects show gamma, REMIND: all gamma analysis reference was set according to the first gd file in the input.')
            #     if 'phy' in self.FileList[0][self.FileList[0].rfind('/'):]:
            #         path23dnrrd = self.FileList[0][:self.FileList[0].rfind('/')] + '/totalnom_phys.nrrd'
            #     else:
            #         path23dnrrd = self.FileList[0][:self.FileList[0].rfind('/')] + '/totalnom_bio.nrrd'
            #     for path2gd in self.FileList[1:]:
            #         nominname=''
            #         if '3D' in path2gd:
            #             nominname='nom_'
            #         if 'phy' in path2gd[path2gd.rfind('/'):]:
            #             path24Dnrrd = path2gd[:path2gd.rfind('/')] + '/total'+nominname + 'phys.nrrd'
            #         elif 'bio' in path2gd[path2gd.rfind('/'):]:
            #             path24Dnrrd = path2gd[:path2gd.rfind('/')] + '/total'+nominname + 'bio.nrrd'
            #
            #         for cri in gammacri:
            #             try:
            #                 gammaresult.append(self.fun_ana_gamma(path23dnrrd, path24Dnrrd,cri))
            #             except:
            #                 writeloginfo = 'check if reference/compare dose nrrd file exist for patient ' + self.patientID + ' plan ' + self.planname
            #                 related_funs.writelog(self.path2log, writeloginfo)
            #
            #     for i in range(0,len(gammacri)):
            #         writegammadata += self.patientID + ' ' + self.planname + ' - - - 0gamma' + gammacri[i] + ' - '
            #         for j in range(0,int(len(gammaresult)/len(gammacri))):
            #             writegammadata+=str(gammaresult[i+2*j][0])+'% '
            #         writegammadata+='\n'
            #         #
            # save info and analysis data
            with open(savedata_fildname, 'w+') as savefileinfo:
                # savefileinfo.writelines('patientID plan VOI volume pre_dose parameter 3D 4D1 4D2 4D3 ...')
                reverseinfo = list(zip(*self.writelinesinfo))
                # reverseinfo = self.writelinesinfo
                for oneline in reverseinfo:
                    for onedata in oneline:
                        savefileinfo.writelines(str(onedata) + ' ')
                    savefileinfo.write('\n')
    def fun_analysis_refonly(self):
        fileNo=0
        fileToanalysis = self.FileList[fileNo]
        Definednameofdata = self.nameofgdlist[fileNo]
        referencedata=True
        VOI_data = self.AnalyzeDVHvoidata(fileToanalysis, Definednameofdata, referencedata)

        # patientIDToW,plannameToW,VOI_names,VOI_volumes,VOI_pres_Dose,VOI_Parameter,VOI_data = \
        #     self.AnalyzeDVHReference(fileToanalysis, Definednameofdata)
        patientIDToW, plannameToW, VOI_names, VOI_volumes, VOI_pres_Dose, VOI_Parameter = self.WriteDVHPreInfo(
            fileToanalysis)
        self.writelinesinfo.append(patientIDToW)
        self.writelinesinfo.append(plannameToW)
        self.writelinesinfo.append(VOI_names)
        self.writelinesinfo.append(VOI_volumes)
        self.writelinesinfo.append(VOI_pres_Dose)
        self.writelinesinfo.append(VOI_Parameter)
        [self.writelinesinfo.append(i) for i in VOI_data]
        savedata_fildname = self.savepath +'reference_analysis' + '.tx'
        with open(savedata_fildname, 'a+') as savefileinfo:
            # savefileinfo.writelines('patientID plan VOI volume pre_dose parameter 3D 4D1 4D2 4D3 ...')
            reverseinfo = list(zip(*self.writelinesinfo))
            # reverseinfo = self.writelinesinfo
            for oneline in reverseinfo:
                for onedata in oneline:
                    savefileinfo.writelines(str(onedata) + ' ')
                savefileinfo.write('\n')
    def WriteDVHPreInfo(self, fileToanalysis):  # return write data: vol, pdose, parameter.
        # get the voiname, voivolume, voiprescirbeddose, voiparameter info
        patientIDToW = ['00_ID']
        plannameToW = ['00_Planname']
        VOI_names = ['00_OAR/Tname']
        VOI_volumes = ['Volume']
        VOI_pres_Dose = ['Dose/Fractions']
        VOI_Parameter = ['parameter']
        VOI_TargetdoseTec = ['TargetdoseTec'] # SIBL,SIBH,Normal
        VOI_OptMethod = ['OptMethod']
        VOI_ionType=['IonType']
        for targetinfo in range(0, len(self.targetnamelist)):
            targetName = self.targetnamelist[targetinfo]

            dIrrVolcc = -1
            dIrrVolcc = self.getIrrVolByVOI(fileToanalysis + self.filesuffixtoattach[0] + '.dvh.gd',
                                            targetName) / 1000  # mm^3 --> cc
            if dIrrVolcc <= 0.0 and not any(c in targetName for c in 'POZM'):
                writeinfo = "Voi " + targetName + " doesn't exist for file: "+fileToanalysis
                related_funs.writelog(self.path2log, writeinfo)

            for i in range(0, len(self.Vxx) + len(self.Dxx) + len(self.Dcc) + 5):  # vxx dxx dcc mean min max ci hi
                patientIDToW.append(self.patientID)
                plannameToW.append(self.planname)
                VOI_names.append(self.targetnamelist[targetinfo])
                VOI_volumes.append(dIrrVolcc)
                VOI_pres_Dose.append(self.targetdoselist[targetinfo])
                VOI_TargetdoseTec.append(self.fun_determinTargetDoseTech(targetinfo))
                VOI_OptMethod.append(self.OptMethod)
                VOI_ionType.append(self.ionType)
            VOI_Parameter.append('min')
            VOI_Parameter.append('max')
            VOI_Parameter.append('mean')
            VOI_Parameter.append('CI')
            VOI_Parameter.append('HI')
            for j in self.Vxx:
                VOI_Parameter.append('V' + str(j) + '%')
            for m in self.Dxx:
                VOI_Parameter.append('D' + str(m) + '%')
            for n in self.Dcc:
                VOI_Parameter.append('D' + str(n) + 'cc')
        for oarinfo in self.oarnamelist:

            dIrrVolcc = -1
            dIrrVolcc = self.getIrrVolByVOI(fileToanalysis + self.filesuffixtoattach[0] + '.dvh.gd',
                                            oarinfo) / 1000  # mm^3 --> cc
            if dIrrVolcc <= 0.0:
                writeinfo = "Voi " + oarinfo + "doesn't exit for file:"+fileToanalysis
                related_funs.writelog(self.path2log, writeinfo)

            for i in range(0, len(self.Dcc) + 1):  # dcc mean
                patientIDToW.append(self.patientID)
                plannameToW.append(self.planname)
                VOI_names.append(oarinfo)
                VOI_volumes.append(str(dIrrVolcc))
                VOI_pres_Dose.append('/' + str(self.fractions) + 'Fxs/')
                VOI_TargetdoseTec.append(self.fun_determinOarDoseTech())
                VOI_OptMethod.append(self.OptMethod)
                VOI_ionType.append(self.ionType)
            VOI_Parameter.append('mean')
            for n in self.Dcc:
                VOI_Parameter.append('D' + str(n) + 'cc')
        # start get dose DVH info.
        return patientIDToW, plannameToW, VOI_names, VOI_volumes, VOI_pres_Dose, VOI_TargetdoseTec,\
            VOI_Parameter,VOI_OptMethod,VOI_ionType
    def AnalyzeDVHvoidata(self, fileToanalysis, Defineddatanameprefix, referencedata):
        # abs means consider Targrt D95 of original plan as 1, calculate percentage different.
        # return write data: vol, pdose, parameter.
        # start get dose DVH info.
        VOI_data = []
        voidata_worst = ['Worst_' + Defineddatanameprefix]
        voidata_mean = ['Mean_' + Defineddatanameprefix]
        voidata_median = ['Median_' + Defineddatanameprefix]
        voidata_SD = ['SD_' + Defineddatanameprefix]
        ContainsReference = False
        startaveragepoint = 0
        if referencedata:
            ContainsReference = True
            startaveragepoint = 1
        for filesuffex in self.filesuffixtoattach:
            print('<info> --> Starting evaluation: '+fileToanalysis+filesuffex)
            gdfilestoanalysis = fileToanalysis + filesuffex + '.dvh.gd'
            VOI_data.append([Defineddatanameprefix + filesuffex[-3:] ])  # add first row(name). for reference colume, it will change acorrodingly.
            countRefereindex = 0

            if self.fractionsacc is not None:  # calculate target dose for each fraction acc.
                fxacc_fileindex = filesuffex[filesuffex.rfind('-') + 1:]
                self.fractions = float(self.fractionsacc) / float(fxacc_fileindex)

            for targetinfo in range(0, len(self.targetnamelist)):
                targetName = self.targetnamelist[targetinfo]
                targetDose = self.targetdoselist[targetinfo]

                # prepareing data for calcuation of CI
                self.lowerdoseforext = float(targetDose)
                self.ExtV95 = self.getExternalV95(gdfilestoanalysis, self.externalname, 'EXT', self.lowerdoseforext)

                Dmin, Dmax, Dmean, CI, HI, Vxxlist, Dxxlist, Dcclist = self.getDVHMetricsFromFileByVOI(
                    gdfilestoanalysis, targetName,
                    'Target', float(targetDose),
                    self.Vxx, self.Dxx,
                    self.Dcc)
                if ContainsReference:
                    VOI_data[-1][0] = Defineddatanameprefix + filesuffex  # modify the first row(name)
                    pass
                VOI_data[-1].append(str('%.4f' % Dmin))
                VOI_data[-1].append(str('%.4f' % Dmax))
                VOI_data[-1].append(str('%.4f' % Dmean))
                VOI_data[-1].append(str('%.4f' % CI))
                VOI_data[-1].append(str('%.4f' % HI))
                for Vxxinfo in Vxxlist:
                    VOI_data[-1].append(str('%.4f' % Vxxinfo))
                for Dxxinfo in Dxxlist:
                    VOI_data[-1].append(str('%.4f' % Dxxinfo))
                for Dccinfo in Dcclist:
                    VOI_data[-1].append(str('%.4f' % Dccinfo))
            for oarinfo in self.oarnamelist:
                targetDose = max(self.targetdoselist)

                Dmin, Dmax, Dmean, CI, HI, Vxxlist, Dxxlist, Dcclist = self.getDVHMetricsFromFileByVOI(
                    gdfilestoanalysis, oarinfo, 'OAR', float(targetDose), self.Vxx, self.Dxx, self.Dcc)
                VOI_data[-1].append(str('%.4f' % Dmean))
                for Dccinfo in Dcclist:
                    VOI_data[-1].append(str('%.4f' % Dccinfo))
            ContainsReference = False
        #if self.robusteva or (self.randomrobust is not None and self.fractionsacc is None):
        if self.robusteva:
            collectedparameters = self.fun_parameterstobeanalysised()
            countofcollecteddata = 0
            for i in range(1, len(VOI_data[0])):  # calculate worst, mean, median, sd
                floatvoidata = []
                for j in range(startaveragepoint, len(VOI_data)):
                    floatvoidata.append(float(VOI_data[j][i]))
                # floatvoidata_neg=[i for i in floatvoidata if float(i)<0]
                voidata_np = np.array(floatvoidata)
                if collectedparameters[countofcollecteddata] == 0:  # min is worst
                    voidata_worst.append(np.min(voidata_np))
                elif collectedparameters[countofcollecteddata] == 1:  # max is worst
                    voidata_worst.append(np.max(voidata_np))
                countofcollecteddata += 1
                voidata_mean.append(np.mean(voidata_np))
                voidata_median.append(np.median(voidata_np))
                voidata_SD.append(np.std(voidata_np))
            VOI_data.append(voidata_worst)
            VOI_data.append(voidata_mean)
            VOI_data.append(voidata_median)
            VOI_data.append(voidata_SD)
        #if (self.robusteva or (self.randomrobust is not None and self.fractionsacc is None)) and not self.Showall:
        if self.robusteva and not self.Showall:
            if referencedata:
        # if data inculdes reference, return reference and worst
                VOI_data_temp = []
                VOI_data_temp.append(VOI_data[0])
                VOI_data_temp.append(voidata_worst)
                return VOI_data_temp
            else:
                return voidata_worst
        else:
            return VOI_data
    def AnalyzeDVHvoirefdata(self, fileToanalysis, Defineddatanameprefix):
        # abs means consider Targrt D95 of original plan as 1, calculate percentage different.
        # return write data: vol, pdose, parameter.
        # start get dose DVH info.
        temppd=self.PlanDose
        tempfractions=self.fractions
        self.PlanDose=self.referenceplanneddose
        self.fractions=self.referencefraction
        VOI_data = []
        print('<info> --> Starting evaluation: '+fileToanalysis)
        gdfilestoanalysis = fileToanalysis + '.dvh.gd'
        VOI_data.append([Defineddatanameprefix])  # add first row(name). for reference colume, it will change acorrodingly.

        for targetinfo in range(0, len(self.targetnamelist)):
            targetName = self.targetnamelist[targetinfo]
            targetDose = self.targetdoselist[targetinfo]
                # prepareing data for calcuation of CI
            self.lowerdoseforext = float(targetDose)
            self.ExtV95 = self.getExternalV95(gdfilestoanalysis, self.externalname, 'EXT', self.lowerdoseforext)
            Dmin, Dmax, Dmean, CI, HI, Vxxlist, Dxxlist, Dcclist = self.getDVHMetricsFromFileByVOI(
                    gdfilestoanalysis, targetName,
                    'Target', float(targetDose),
                    self.Vxx, self.Dxx,
                    self.Dcc)
            VOI_data[-1].append(str('%.4f' % Dmin))
            VOI_data[-1].append(str('%.4f' % Dmax))
            VOI_data[-1].append(str('%.4f' % Dmean))
            VOI_data[-1].append(str('%.4f' % CI))
            VOI_data[-1].append(str('%.4f' % HI))
            for Vxxinfo in Vxxlist:
                VOI_data[-1].append(str('%.4f' % Vxxinfo))
            for Dxxinfo in Dxxlist:
                VOI_data[-1].append(str('%.4f' % Dxxinfo))
            for Dccinfo in Dcclist:
                VOI_data[-1].append(str('%.4f' % Dccinfo))
        for oarinfo in self.oarnamelist:
            targetDose = max(self.targetdoselist)

            Dmin, Dmax, Dmean, CI, HI, Vxxlist, Dxxlist, Dcclist = self.getDVHMetricsFromFileByVOI(
                gdfilestoanalysis, oarinfo, 'OAR', float(targetDose), self.Vxx, self.Dxx, self.Dcc)
            VOI_data[-1].append(str('%.4f' % Dmean))
            for Dccinfo in Dcclist:
                VOI_data[-1].append(str('%.4f' % Dccinfo))

        self.PlanDose=temppd
        self.fractions=tempfractions

        return VOI_data
    def getDVHMetricsFromFileByVOI(self, filename, voiname, voitype, voidose, voiVxx, voiDxx, voiDcc):

        if voidose != self.PlanDose:
            dTarget_Fraction_Dose = float(voidose) / float(self.fractions)
            coefficient = float(self.PlanDose) / dTarget_Fraction_Dose
        else:
            coefficient = 1

        bGotIrrVol = False
        bVOIinData = False
        dIrrVolcc = -1.0

        Dmin = -1
        Dmax = -1
        Dmean = -1
        CI = -1
        HI = -1
        Vxxlist = []
        Dxxlist = []
        Dcclist = []

        xvalues = []
        yvalues = []

        bVOIinData = True
        bGotIrrVol = False
        dIrrVolcc = self.getIrrVolByVOI(filename, voiname) / 1000  # mm^3 --> cc

        if dIrrVolcc > 0.0:
            bGotIrrVol = True
            # print("Voi: "+voiname+" has volume "+str(dIrrVolcc)+" cc")
            xvalues, yvalues = self.getDVHValuesFromFileByVOI(filename, voiname)

        xvalues = np.array(xvalues)
        yvalues = np.array(yvalues)
        xvalues = xvalues * coefficient

        Dmin = self.Dmin(voiname, filename) * float(self.PlanDose) * float(self.fractions)
        Dmean = self.Dmean(voiname, filename) * float(self.PlanDose) * float(self.fractions)
        Dmax = self.Dmax(voiname, filename) * float(self.PlanDose) * float(self.fractions)

        DVal = 0
        if (voitype == "Target"):
            # finding D95 & D5
            # D5  = self.D_n(5.0,[yvalues],xvalues)[0]
            # D95 = self.D_n(95,[yvalues],xvalues)[0]
            D5 = self.fun_Dxx(5, xvalues, yvalues, voiname)
            D95 = self.fun_Dxx(95, xvalues, yvalues, voiname)
            HI = D5 - D95

            if dIrrVolcc != 0:
                CI = self.ExtV95 / float(dIrrVolcc)
                # CN = float(V95 * V95 / (CI * 100.))
            else:
                print("Calculation of CI/CN not possible. CI=0 cases!")
                CI = 'NA'
                writeinfo = "Voi " + voiname + " CI = NA"
                related_funs.writelog(self.path2log, writeinfo)
        for temp in voiVxx:
            try:
                Vxxlist.append(self.fun_Vxx(temp, xvalues, yvalues, voiname))
            except:
                Vxxlist.append(-1)
        for temp in voiDxx:
            DVal = temp
            try:
                Dxxlist.append(self.fun_Dxx(temp, xvalues, yvalues, voiname))
            except:
                Dxxlist.append(-1)
        for temp in voiDcc:
            Dval = float(temp) / float(dIrrVolcc) * 100.00
            # try:
            #Dcclist.append(self.fun_Dxx(Dval, xvalues, yvalues, voiname) * float(voidose) / 100.00)
            Dcclist.append(self.fun_Dxx(Dval, xvalues, yvalues, voiname))
            # except:
            #    Dcclist.append(-1)
        # print(Dcclist)
        return Dmin, Dmax, Dmean, CI, HI, Vxxlist, Dxxlist, Dcclist

    def fun_determinTargetDoseTech(self,targetinfo):
        # from given target info number and self.targetdoselist, determin target belongs to SIBH, SIBL or Normal
        NoofdosesinPD=len(set(self.targetdoselist))
        neckintarget=False
        for targetname in self.targetnamelist:
            if 'neck' in targetname:
                neckintarget=True
                break
        if NoofdosesinPD==1: # ctv or gtv case
            TargetDoseTech='Normal'
        elif NoofdosesinPD==2: # ctv and gtv case
            if self.targetdoselist[targetinfo] == max(self.targetdoselist):
                TargetDoseTech = 'SIBH'
            else:
                TargetDoseTech = 'SIBL'
        elif NoofdosesinPD==3 and not neckintarget: # gboost, gtv, ctv case. ~!neck case not included.
            if self.targetdoselist[targetinfo] == min(self.targetdoselist):
                TargetDoseTech = 'SIBL'
            else:
                TargetDoseTech = 'SIBH'
        elif NoofdosesinPD==3 and neckintarget: # gtv, ctv, ctvneck case.
            if self.targetdoselist[targetinfo] == min(self.targetdoselist):
                TargetDoseTech = 'Normal'
            elif self.targetdoselist[targetinfo] == max(self.targetdoselist):
                TargetDoseTech = 'SIBH'
            else:
                TargetDoseTech = 'SIBL'
        elif NoofdosesinPD==4:
            non_duplicate_sorted_list = sorted(set(self.targetdoselist))
            if self.targetdoselist[targetinfo] == max(self.targetdoselist):
                TargetDoseTech = 'SIBH'
            elif self.targetdoselist[targetinfo] == min(self.targetdoselist):
                TargetDoseTech = 'SIBL'
            elif self.targetdoselist[targetinfo]==non_duplicate_sorted_list[1]:
                TargetDoseTech = 'SIBL'
            elif self.targetdoselist[targetinfo]==non_duplicate_sorted_list[2]:
                TargetDoseTech = 'SIBH'
        else:
            TargetDoseTech='9999'
            print('error in counting the target doses')
        return TargetDoseTech
    def fun_determinOarDoseTech(self):
        # from given target info number and self.targetdoselist, determin target belongs to SIBH, SIBL or Normal
        NoofdosesinPD=len(set(self.targetdoselist))
        if NoofdosesinPD==1:
            OARDoseTech = 'Normal'
        else:
            OARDoseTech = 'SIB'
        return OARDoseTech
    def getExternalV95(self, filename, voiname, voitype, voidose):
        # print voiname
        if voidose != self.PlanDose:
            dTarget_Fraction_Dose = float(voidose) / float(self.fractions)
            coefficient = float(self.PlanDose) / dTarget_Fraction_Dose
        else:
            coefficient = 1
        xvalues = []
        yvalues = []
        ExtV95 = -1
        dIrrVolcc = self.getIrrVolByVOI(filename, voiname) / 1000  # mm^3 --> cc
        if dIrrVolcc > 0.0:
            xvalues, yvalues = self.getDVHValuesFromFileByVOI(filename, voiname)
        xvalues = np.array(xvalues)
        yvalues = np.array(yvalues)
        xvalues = xvalues * coefficient

        if (voitype != "EXT"):
            print('wrong input, this is for External structure only')
        else:
            percen = self.fun_Vxx(95, xvalues, yvalues, voiname)
            ExtV95 = percen / 100 * dIrrVolcc  # abs v95%
        return ExtV95
    def getDVHValuesFromFileByVOI(self, filename, VOIstr, **kwargs):
        bReadDifferential = kwargs.get("bReadDifferential", False)
        """Read x and y values from GD file for DVH by VOIstr"""
        with open(filename, 'r') as fin:
            startline = 0
            stopline = 0
            stoplines = []
            j = 0
            for line in fin:
                j += 1
                if re.match("c:[ ]*{}".format(VOIstr + ' '), line):
                    startline = j + 2
                if re.match("c:[ ]*VOI", line):
                    stoplines.append(j - 1)
            for elem in stoplines:
                if elem > startline:
                    stopline = elem
                    break
            if stopline == 0:
                stopline = j
            if startline == 0:
                print("VOI %s not found!" % VOIstr)
                exit(1)
        with open(filename, 'r') as fin:
            xvalues = []
            yvalues = []
            jj = 0
            for line in fin:
                jj += 1
                if jj < startline or jj > stopline: continue
                values = line.split()
                if len(values) == 3:
                    xvalues.append(float(values[0]))
                    if bReadDifferential:
                        yvalues.append(float(values[2]))
                    else:
                        yvalues.append(float(values[1]))

        return (xvalues, yvalues)
    def Dmin(self, VOIstr, filename):
        mindose=9999
        with open(filename, 'r') as fin:
            for line in fin:
                b = line.split()
                if line.startswith("c: " + VOIstr) and b[1] == VOIstr:
                    mindose = float(line.split()[2])
            result = mindose
        return result

    def Dmax(self, VOIstr, filename):
        maxdose=9999
        with open(filename, 'r') as fin:
            for line in fin:
                b = line.split()
                if line.startswith("c: " + VOIstr) and b[1] == VOIstr:
                    maxdose = float(line.split()[3])
            result = maxdose
        return result

    def Dmean(self, VOIstr, path):
        meandose=9999
        with open(path, 'r') as fin:
            for line in fin:
                b = line.split()
                if line.startswith("c: " + VOIstr) and b[1] == VOIstr:
                    meandose = float(line.split()[4])
            result = meandose
        return result

    def getIrrVol(self, path):
        with open(path, 'r') as fin:
            for line in fin:
                values = line.split()
                if len(values) != 3:
                    continue
                else:
                    irrvol = values[1]
                    break
            result = irrvol
        return result

    def getIrrVolByVOI(self, path, VOIstr):
        irrvol = -1
        try:
            with open(path, 'r') as fin:
                for line in fin:
                    b=line.split()
                    if line.startswith("c: " + VOIstr) and b[1]==VOIstr:
                        irrvol = float(line.split()[7])
                        break
        except:
            pass
        return irrvol

    def V(self, n, path):
        result = []
        with open(path, 'r') as fin:
            for line in fin:
                values = line.split()
                if len(values) != 3:
                    continue
                else:
                    values[0] = int(values[0])
                    values[1] = float(values[1])
                if values[0] == n:
                    relevantDose = values[0]
                    relevantVolume = values[1]
                    result = relevantVolume
        return result

    def D(self, n, path):
        result = []
        with open(path, 'r') as fin:
            for line in fin:
                values = line.split()
                if len(values) != 3:
                    continue
                else:
                    values[0] = float(values[0])
                    values[1] = float(values[1])
                if values[1] < n:
                    smallervaluevolume = values[1]
                    smallervaluedose = values[0]
                if values[1] > n:
                    biggervaluevolume = values[1]
                    biggervaluedose = values[0]
                    continue
                finalvalue = (smallervaluedose - biggervaluedose) / (biggervaluevolume - smallervaluevolume) * (
                            n - smallervaluevolume) + smallervaluedose
                result = finalvalue  # print str(biggervaluevolume) + " " + str(smallervaluevolume)
                # print str(biggervaluedose) + " " + str(smallervaluedose)
                # print str(finalvalue)
                break
        return result

    def D_n(self, nn, yDataList, xData):
        nn = float(nn)
        D_left = []
        D_right = []
        for elem in yDataList:
            bDataHIprocessable = False
            for j, e in list(enumerate(elem)):
                if e >= nn:
                    bDataHIprocessable = True
            if not bDataHIprocessable:
                print("D%s not possible to calculate! Is it an OAR?!?!?!" % nn)
                continue
            for j, e in list(enumerate(elem)):
                if e <= nn:
                    D_right.append([xData[j], e])
                    break
            for j, e in reversed(list(enumerate(elem))):
                if e >= nn:
                    D_left.append([xData[j], e])
                    break

        # interpolation between D_left and D_right
        if len(D_right) != len(D_left):
            print("List length from D_left and D_right don't match.")
            exit(1)
        else:
            Dnn = []
            for i in range(len(D_left)):
                if D_right[i][1] == D_left[i][1]:
                    Dnn.append(D_left[i][0])
                else:
                    m = (D_right[i][1] - D_left[i][1]) / (D_right[i][0] - D_left[i][0])
                    n = D_left[i][1] - m * D_left[i][0]
                    Dnn.append((nn - n) / m)
            return Dnn

    def V_n(self, nn, yDataList, xData):
        nn = float(nn)
        V_left = []
        V_right = []
        for elem in yDataList:
            bDataHIprocessable = False
            for j, e in list(enumerate(elem)):
                if e >= nn:
                    bDataHIprocessable = True
            if not bDataHIprocessable:
                print("D%s not possible to calculate! Is it an OAR?!?!?!" % nn)
                continue
            for j, e in list(enumerate(elem)):
                if e <= nn:
                    V_right.append([xData[j], e])
                    break
            for j, e in reversed(list(enumerate(elem))):
                if e >= nn:
                    V_left.append([xData[j], e])
                    break

        # interpolation between D_left and D_right
        if len(V_right) != len(V_left):
            print("List length from D_left and D_right don't match.")
            exit(1)
        else:
            Vnn = []
            for i in range(len(V_left)):
                if V_right[i][1] == V_left[i][1]:
                    Vnn.append(V_left[i][0])
                else:
                    m = (V_right[i][1] - V_left[i][1]) / (V_right[i][0] - V_left[i][0])
                    n = V_left[i][1] - m * V_left[i][0]
                    Vnn.append((nn - n) / m)
            return Vnn

    def fun_Vxx(self, xx_per, Dose_values, Volume_values, voiname):
        Vxx = -1
        for i in range(0, len(Dose_values)):
            if Dose_values[i] == xx_per:
                Vxx = Volume_values[i]
                break
            elif Dose_values[i] < xx_per:
                continue
            elif Dose_values[i] > xx_per:
                Vxx = self.fun_trend(float(xx_per), Dose_values[i - 1], Dose_values[i], Volume_values[i - 1],
                                     Volume_values[i])
                break
        if Vxx == -1:
            # print('V',xx_per,' of ' + voiname + ' is not possible to calculate')
            pass
        return Vxx

    def fun_Dxx(self, xx_per, Dose_values, Volume_values, voiname):
        Dxx = -1
        for i in range(0, len(Volume_values)):
            if Volume_values[i] == xx_per:
                Dxx = Dose_values[i]
                break
            elif Volume_values[i] > xx_per:
                continue
            elif Volume_values[i] <= xx_per:
                Dxx = self.fun_trend(float(xx_per), Volume_values[i - 1], Volume_values[i], Dose_values[i - 1],
                                     Dose_values[i])
                break
        if Dxx == -1:
            # print('D',xx_per,' for '+voiname+' is not possible to calculate')
            pass
        return Dxx

    def fun_trend(self, newX, X1, X2, Y1, Y2):
        newY = Y1 + (Y2 - Y1) * (newX - X1) / (X2 - X1)
        return newY

    def fun_ana_gamma(self, path2dose1, path2dose2, gammacri):
        gammaana = gamma.class_gammaanalysis()
        criterialist, gammalist = gammaana.fun_gamma_analysis(additionalinfo='', dose1=path2dose1, dose2=path2dose2,
                                                              dosediscrit=gammacri, cuoff=10, maxdose='global',
                                                              interfra=5, maxgamma=1.1, fraction=1, saveresults='False',
                                                              pronecase=False, moreinfo=False)
        return gammalist

    def fun_calculateDevPerc(self, ref, com):  #
        if float(ref) != 0:
            return (float(com) - float(ref)) / float(ref)
        else:
            return 9999
    def fun_calculateDev_abs(self, ref, com):  # consider reference as 100%. return 100*A/ref
        if float(ref) != 0:
            return 100*float(com)/float(ref)
        else:
            return 9999

    # according to target and oar and Vxx, Dxx, Dcc. return 1 and 0 matrix,
    # 1 means max is worst, 0 means min is the worst.
    def fun_parameterstobeanalysised(self):
        whoistheworst=[]
        for targetinfo in range(0, len(self.targetnamelist)):
            # one target: Min, Max, Mean, CI, HI,
            whoistheworst.append(0)
            whoistheworst.append(1)
            whoistheworst.append(0)
            whoistheworst.append(1)
            whoistheworst.append(1)
            for Vxxinfo in self.Vxx:
                if float(Vxxinfo)<=100:
                    whoistheworst.append(0)
                else:
                    whoistheworst.append(1)
            for Dxxinfo in self.Dxx:
                if float(Dxxinfo)>=50:
                    whoistheworst.append(0)
                else:
                    whoistheworst.append(1)
            for Dccinfo in self.Dcc: # D1cc max is the worst
                whoistheworst.append(1)
        for oarinfo in range(0, len(self.oarnamelist)):
            # for oar mean
            whoistheworst.append(1)
            for Dccinfo in self.Dcc:  # D1cc max is the worst
                whoistheworst.append(1)
        return whoistheworst


