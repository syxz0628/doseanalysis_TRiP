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
    def __init__(self, patientID, planname, targetnamelist, targetdoselist, oarnamelist, externalname, fractions,
                 savepath, gammaEva, robustevaluation, path2gdlist, nameofgdlist):
        self.patientID = patientID
        self.planname = planname
        self.targetnamelist = targetnamelist
        self.targetdoselist = targetdoselist
        self.oarnamelist = oarnamelist
        self.externalname = externalname
        self.lowerdoseforext = 999  # take lower prescribed dose for calculate of V95 for Ext
        self.ExtV95 = 0
        self.fractions = fractions
        self.savepath = savepath
        self.gammaEva = gammaEva
        self.robustevaluation = robustevaluation
        self.robusteva = True
        if self.robustevaluation == '21':
            self.robust_suffix = ['_nom', '_nx', '_ny', '_nz', '_px', '_py', '_pz',
                                  '_hd', '_hd_nx', '_hd_ny', '_hd_nz', '_hd_px', '_hd_py', '_hd_pz',
                                  '_ld', '_ld_nx', '_ld_ny', '_ld_nz', '_ld_px', '_ld_py', '_ld_pz']
        elif self.robustevaluation == '9':
            self.robust_suffix = ['_nom', '_hd', '_ld', '_nx', '_ny', '_nz', '_px', '_py', '_pz']
        else:
            self.robusteva=False
            self.robust_suffix = ['']

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

    def fun_analysis_gd(self, dose_shown_in_gd):
        # dose_shown_in_gd is the dose written in the exec file, for SPHIC momi cases this value was set to 3 for all plans.
        self.PlanDose = dose_shown_in_gd
        if len(self.nameofgdlist) != len(self.FileList):
            errormess = 'Detects wrong input:"gdlist path and name is not the same length, please check"'
            related_funs.writelog(self.path2log, errormess)
            sys.exit()
        related_funs.writelog(self.path2log, 'Start a new analysis')
        savedata_fildname = self.savepath + self.patientID + '_' + self.planname + '_' + \
                            "_".join(m for m in self.nameofgdlist) + '.txt'
        # write log
        writeloginfo = 'running patient: ' + self.patientID + ' plan: ' + self.planname
        related_funs.writelog(self.path2log, writeloginfo)
        referencedata = True
        for fileNo in range(0, len(self.FileList)):
            fileToanalysis = self.FileList[fileNo]
            Definednameofdata = self.nameofgdlist[fileNo]

            VOI_data = self.AnalyzeDVHvoidata(fileToanalysis, Definednameofdata, referencedata)
            if referencedata:
                # patientIDToW,plannameToW,VOI_names,VOI_volumes,VOI_pres_Dose,VOI_Parameter,VOI_data = \
                #     self.AnalyzeDVHReference(fileToanalysis, Definednameofdata)
                patientIDToW, plannameToW, VOI_names, VOI_volumes, VOI_pres_Dose, VOI_Parameter = self.AnalyzeDVHPre(
                    fileToanalysis)
                self.writelinesinfo.append(patientIDToW)
                self.writelinesinfo.append(plannameToW)
                self.writelinesinfo.append(VOI_names)
                self.writelinesinfo.append(VOI_volumes)
                self.writelinesinfo.append(VOI_pres_Dose)
                self.writelinesinfo.append(VOI_Parameter)
                referencedata = False
            self.fun_append_listonebyone(self.writelinesinfo, VOI_data)

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

    def fun_analysis_refonly(self,dose_shown_in_gd):
        # dose_shown_in_gd is the dose written in the exec file, for SPHIC momi cases this value was set to 3 for all plans.
        self.PlanDose = dose_shown_in_gd
        fileNo=0
        fileToanalysis = self.FileList[fileNo]
        Definednameofdata = self.nameofgdlist[fileNo]
        referencedata=True
        VOI_data = self.AnalyzeDVHvoidata(fileToanalysis, Definednameofdata, referencedata)

        # patientIDToW,plannameToW,VOI_names,VOI_volumes,VOI_pres_Dose,VOI_Parameter,VOI_data = \
        #     self.AnalyzeDVHReference(fileToanalysis, Definednameofdata)
        patientIDToW, plannameToW, VOI_names, VOI_volumes, VOI_pres_Dose, VOI_Parameter = self.AnalyzeDVHPre(
            fileToanalysis)
        self.writelinesinfo.append(patientIDToW)
        self.writelinesinfo.append(plannameToW)
        self.writelinesinfo.append(VOI_names)
        self.writelinesinfo.append(VOI_volumes)
        self.writelinesinfo.append(VOI_pres_Dose)
        self.writelinesinfo.append(VOI_Parameter)
        self.fun_append_listonebyone(self.writelinesinfo, VOI_data)
        savedata_fildname = self.savepath +'reference_analysis' + '.tx'
        with open(savedata_fildname, 'a+') as savefileinfo:
            # savefileinfo.writelines('patientID plan VOI volume pre_dose parameter 3D 4D1 4D2 4D3 ...')
            reverseinfo = list(zip(*self.writelinesinfo))
            # reverseinfo = self.writelinesinfo
            for oneline in reverseinfo:
                for onedata in oneline:
                    savefileinfo.writelines(str(onedata) + ' ')
                savefileinfo.write('\n')

    def AnalyzeDVHPre(self, fileToanalysis):  # return write data: vol, pdose, parameter.
        # get the voiname, voivolume, voiprescirbeddose, voiparameter info
        patientIDToW = ['ID']
        plannameToW = ['Planname']
        VOI_names = ['OAR/Tname']
        VOI_volumes = ['Volume']
        VOI_pres_Dose = ['Dose/Fractions']
        VOI_Parameter = ['parameter']
        for targetinfo in range(0, len(self.targetnamelist)):
            targetName = self.targetnamelist[targetinfo]

            dIrrVolcc = -1
            dIrrVolcc = self.getIrrVolByVOI(fileToanalysis + self.robust_suffix[0] + '.dvh.gd',
                                            targetName) / 1000  # mm^3 --> cc
            if dIrrVolcc <= 0.0:
                writeinfo = "Voi " + targetName + " did not receive any dose/doesn't exist."
                related_funs.writelog(self.path2log, writeinfo)

            for i in range(0, len(self.Vxx) + len(self.Dxx) + len(self.Dcc) + 5):  # vxx dxx dcc
                patientIDToW.append(self.patientID)
                plannameToW.append(self.planname)
                VOI_names.append(self.targetnamelist[targetinfo])
                VOI_volumes.append(dIrrVolcc)
                VOI_pres_Dose.append(self.targetdoselist[targetinfo])

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
            dIrrVolcc = self.getIrrVolByVOI(fileToanalysis + self.robust_suffix[0] + '.dvh.gd',
                                            oarinfo) / 1000  # mm^3 --> cc
            if dIrrVolcc <= 0.0:
                writeinfo = "Voi " + oarinfo + " did not receive any dose/doesn't exit."
                related_funs.writelog(self.path2log, writeinfo)

            for i in range(0, len(self.Dcc) + 1):  # dcc mean
                patientIDToW.append(self.patientID)
                plannameToW.append(self.planname)
                VOI_names.append(oarinfo)
                VOI_volumes.append(str(dIrrVolcc))
                VOI_pres_Dose.append('/' + str(self.fractions) + 'Fxs/')
            VOI_Parameter.append('mean')
            for n in self.Dcc:
                VOI_Parameter.append('D' + str(n) + 'cc')
        # start get dose DVH info.
        return patientIDToW, plannameToW, VOI_names, VOI_volumes, VOI_pres_Dose, VOI_Parameter

    def AnalyzeDVHvoidata(self, fileToanalysis, Definednameofdata,
                          referencedata):  # return write data: vol, pdose, parameter.
        # start get dose DVH info.
        VOI_data = []
        voidata_worst = ['Worst_' + Definednameofdata]
        voidata_mean = ['Mean_' + Definednameofdata]
        voidata_median = ['Median_' + Definednameofdata]
        voidata_SD = ['SD_' + Definednameofdata]
        ContainsReference = False
        startaveragepoint = 0
        if referencedata:
            ContainsReference = True
            startaveragepoint = 1

        for filesuffex in self.robust_suffix:
            gdfilestoanalysis = fileToanalysis + filesuffex + '.dvh.gd'
            VOI_data.append([Definednameofdata + filesuffex + '_dev'])  # add first row(name). for reference colume, it will change acorrodingly.
            countRefereindex = 0
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
                    VOI_data[-1][0] = Definednameofdata + filesuffex  # modify the first row(name)
                    pass
                else:
                    Dmin = self.fun_calculateDevPerc(self.referenceDATAforCompare[countRefereindex], Dmin)
                    countRefereindex += 1
                    Dmax = self.fun_calculateDevPerc(self.referenceDATAforCompare[countRefereindex], Dmax)
                    countRefereindex += 1
                    Dmean = self.fun_calculateDevPerc(self.referenceDATAforCompare[countRefereindex], Dmean)
                    countRefereindex += 1
                    CI = self.fun_calculateDevPerc(self.referenceDATAforCompare[countRefereindex], CI)
                    countRefereindex += 1
                    HI = self.fun_calculateDevPerc(self.referenceDATAforCompare[countRefereindex], HI)
                    countRefereindex += 1
                    for vxxno in range(0, len(Vxxlist)):
                        if Vxxlist[vxxno] == -1:
                            print('V', self.Vxx[vxxno], ' of ' + targetName + ' is not possible to calculate')
                        Vxxlist[vxxno] = self.fun_calculateDevPerc(self.referenceDATAforCompare[countRefereindex],
                                                                   Vxxlist[vxxno])
                        countRefereindex += 1
                    for Dxxno in range(0, len(Dxxlist)):
                        if Dxxlist[Dxxno] == -1:
                            print('V', self.Dxx[Dxxno], ' of ' + targetName + ' is not possible to calculate')
                        Dxxlist[Dxxno] = self.fun_calculateDevPerc(self.referenceDATAforCompare[countRefereindex],
                                                                   Dxxlist[Dxxno])
                        countRefereindex += 1
                    for Dccno in range(0, len(Dcclist)):
                        if Dcclist[Dccno] == -1:
                            print('V', self.Dcc[Dccno], ' of ' + targetName + ' is not possible to calculate')
                        Dcclist[Dccno] = self.fun_calculateDevPerc(self.referenceDATAforCompare[countRefereindex],
                                                                  Dcclist[Dccno])
                        countRefereindex += 1

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

                if ContainsReference:
                    self.referenceDATAforCompare.append(Dmin)
                    self.referenceDATAforCompare.append(Dmax)
                    self.referenceDATAforCompare.append(Dmean)
                    self.referenceDATAforCompare.append(CI)
                    self.referenceDATAforCompare.append(HI)
                    for Vxxinfo in Vxxlist:
                        self.referenceDATAforCompare.append(Vxxinfo)
                    for Dxxinfo in Dxxlist:
                        self.referenceDATAforCompare.append(Dxxinfo)
                    for Dccinfo in Dcclist:
                        self.referenceDATAforCompare.append(Dccinfo)

            for oarinfo in self.oarnamelist:
                targetDose = max(self.targetdoselist)
                Dmin, Dmax, Dmean, CI, HI, Vxxlist, Dxxlist, Dcclist = self.getDVHMetricsFromFileByVOI(
                    gdfilestoanalysis,
                    oarinfo, 'OAR',
                    float(targetDose),
                    self.Vxx, self.Dxx,
                    self.Dcc)
                if ContainsReference:
                    pass
                else:
                    Dmean = self.fun_calculateDevPerc(self.referenceDATAforCompare[countRefereindex], Dmean)
                    countRefereindex += 1
                    for Dccno in range(0, len(Dcclist)):
                        Dcclist[Dccno] = self.fun_calculateDevPerc(self.referenceDATAforCompare[countRefereindex],
                                                                   Dcclist[Dccno])
                        countRefereindex += 1

                VOI_data[-1].append(str('%.4f' % Dmean))
                for Dccinfo in Dcclist:
                    VOI_data[-1].append(str('%.4f' % Dccinfo))

                if ContainsReference:
                    self.referenceDATAforCompare.append(Dmean)
                    for Dccinfo in Dcclist:
                        self.referenceDATAforCompare.append(Dccinfo)
                    # reference data start average from voidata 1. skip the reference.

            ContainsReference = False
        if self.robusteva:
            for i in range(1, len(VOI_data[0])):  # calculate worst, mean, median, sd
                floatvoidata = []
                for j in range(startaveragepoint, len(VOI_data)):
                    floatvoidata.append(float(VOI_data[j][i]))
                voidata_np = np.array(floatvoidata)
                voidata_worst.append(related_funs.lambda_abs_max(voidata_np, 0, np.abs))
                voidata_mean.append(np.mean(voidata_np))
                voidata_median.append(np.median(voidata_np))
                voidata_SD.append(np.std(voidata_np))
            VOI_data.append(voidata_worst)
            VOI_data.append(voidata_mean)
            VOI_data.append(voidata_median)
            VOI_data.append(voidata_SD)
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
        with open(filename, 'r') as fin:
            for line in fin:
                b = line.split()
                if line.startswith("c: " + VOIstr) and b[1] == VOIstr:
                    maxdose = float(line.split()[2])
            result = maxdose
        return result

    def Dmax(self, VOIstr, filename):
        with open(filename, 'r') as fin:
            for line in fin:
                b = line.split()
                if line.startswith("c: " + VOIstr) and b[1] == VOIstr:
                    maxdose = float(line.split()[3])
            result = maxdose
        return result

    def Dmean(self, VOIstr, path):
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
            writeinfo = "Voi " + VOIstr + ' of ' + path + ' did not receive any dose/doesn\'t exist.'
            related_funs.writelog(self.path2log, writeinfo)
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

    def fun_append_listonebyone(self, write2list, datalist):
        for list2write in datalist:
            write2list.append(list2write)
        return write2list

    def fun_calculateDevPerc(self, ref, com):  #
        if float(ref) != 0:
            return (float(com) - float(ref)) / float(ref)
        else:
            return 9999
