#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 13:58:34 2022

@author: ysheng
"""
import re
# import sys,glob ## command line stuff
import numpy as np

import related_funs
import gamma

class class_analysis_gd:
    def __init__(self, patientID,planname,targetnamelist,targetdoselist,oarnamelist,externalname,fractions,path2gdlist,savename,savepath):
        self.patientID=patientID
        self.planname=planname
        self.targetnamelist=targetnamelist
        self.targetdoselist=targetdoselist
        self.oarnamelist=oarnamelist
        self.externalname=externalname
        self.lowerdoseforext=999 #take lower prescribed dose for calculate of V95 for Ext
        self.ExtV95=0
        self.fractions=fractions
        self.FileList=path2gdlist
        self.savename=savename
        self.savepath=savepath

        if self.savename==None:
            self.savename=''
        if self.savepath == None:
            self.savepath='./'
        self.Vxx=[90,95,100,105]
        self.Dxx=[5,95,98]
        self.Dcc=[1]
        self.path2log=self.savepath+'00_Doseana_processing.log'
        #self.path2log = '/home/yurii/Sheng/patient_data/00_Doseana_processing.log'

        self.voilist = []
        self.VOI_names=[]
        self.VOI_volumes=[]
        self.VOI_pres_Dose=[]
        self.VOI_Parameter=[]
        self.VOI_data=[] # 2dimension, [1st all voi-parameter data][2nd all voi-parameter data][...]
        self.writelinesinfo=[]
        
        
    def fun_analysis_gd(self,dose_shown_in_gd,gamma):#dose_shown_in_gd is the dose written in the exec file, for SPHIC momi cases this value was set to 3 for all plans.
        self.PlanDose = dose_shown_in_gd
        related_funs.writelog(self.path2log, 'Start a new analysis')
        # add target and voi all info to voilist
        for targetinfo in range(0,len( self.targetnamelist)):
            self.SetVOI(voiname=self.targetnamelist[targetinfo],voidose=self.targetdoselist[targetinfo],\
                        voitype='Target',Vxx=self.Vxx,Dxx=self.Dxx,Dcc=self.Dcc)
        for oarinfo in self.oarnamelist:
            self.SetVOI(voiname=oarinfo,voitype='OAR',voidose=float(self.PlanDose) * float(self.fractions),Dcc=self.Dcc)
        self.SetVOI(voiname=self.externalname, voitype='EXT', voidose=float(self.PlanDose) * float(self.fractions), Vxx=self.Vxx)
        # path of save data .txt
        # save in specific patient folder
        # savedata_fildname='/u/ysheng/MyAIXd/projects/patients/'+self.patientID+'/dose_ana_'+self.savename+'_'+self.patientID+'_'+self.planname+'.txt'
        # save in command/dose-compare_logs folder
        savedata_fildname=self.savepath+self.patientID+'_'+self.planname+'_' + self.savename+'.txt'
        # write log
        writeloginfo='running patient: '+self.patientID+' plan: '+self.planname
        related_funs.writelog(self.path2log, writeloginfo)
        # write analysis data to line

        writegddata=self.AnalyzeDVHs()

        # write gamma data to line
        gammacri = []
        gammaresult = []
        writegammadata=''
        gammacri=['3/3','2/2']
        if gamma:
            print('detects show gamma')
            if 'phy' in self.FileList[0][self.FileList[0].rfind('/'):]:
                path23dnrrd = self.FileList[0][:self.FileList[0].rfind('/')] + '/totalnom_phys.nrrd'
            else:
                path23dnrrd = self.FileList[0][:self.FileList[0].rfind('/')] + '/totalnom_bio.nrrd'
            for path2gd in self.FileList[1:]:
                nominname=''
                if '3D' in path2gd:
                    nominname='nom_'
                if 'phy' in path2gd[path2gd.rfind('/'):]:
                    path24Dnrrd = path2gd[:path2gd.rfind('/')] + '/total'+nominname + 'phys.nrrd'
                elif 'bio' in path2gd[path2gd.rfind('/'):]:
                    path24Dnrrd = path2gd[:path2gd.rfind('/')] + '/total'+nominname + 'bio.nrrd'

                for cri in gammacri:
                    try:
                        gammaresult.append(self.fun_ana_gamma(path23dnrrd, path24Dnrrd,cri))
                    except:
                        writeloginfo = 'check if reference/compare dose nrrd file exist for patient ' + self.patientID + ' plan ' + self.planname
                        related_funs.writelog(self.path2log, writeloginfo)

            for i in range(0,len(gammacri)):
                writegammadata += self.patientID + ' ' + self.planname + ' - - - 0gamma' + gammacri[i] + ' - '
                for j in range(0,int(len(gammaresult)/len(gammacri))):
                    writegammadata+=str(gammaresult[i+2*j][0])+'% '
                writegammadata+='\n'
                #
        # save info and analysis data
        with open (savedata_fildname,'w+') as savefileinfo:
            #savefileinfo.writelines('patientID plan VOI volume pre_dose parameter 3D 4D1 4D2 4D3 ...')
            for oneline in writegddata:
                savefileinfo.writelines(self.patientID+' '+self.planname+' ')
                for onedata in oneline:
                    savefileinfo.writelines(str(onedata)+' ')
                savefileinfo.write('\n')
            savefileinfo.writelines(writegammadata)
            
    def SetVOI(self,voiname, voitype, **kwargs):
        voiVxx = kwargs.get('Vxx', self.Vxx)
        voiDxx = kwargs.get('Dxx', self.Dxx)
        voiDcc = kwargs.get('Dcc', self.Dcc)
        voidose = kwargs.get('voidose', self.PlanDose)
        self.voilist.append([voiname,voitype,voidose,voiVxx,voiDxx,voiDcc])
        
    def AnalyzeDVHs(self): # return write data : vol,pdose,parameter,3D,4D1,4D2...
        #get the voiname, voivolume, voiprescirbeddose, voiparameter info
        for voiname,voitype,voidose,voiVxx,voiDxx,voiDcc in self.voilist:
            dIrrVolcc = -1
            dIrrVolcc = self.getIrrVolByVOI(self.FileList[0],voiname)/1000 # mm^3 --> cc
            if dIrrVolcc <= 0.0:
                writeinfo="Voi "+voiname+" did not receive any dose/doesn't exit."
                related_funs.writelog(self.path2log, writeinfo)
            if voitype=='Target':
                for i in range (0, len(voiVxx)+len(voiDxx)+len(voiDcc)+3+2): # vxx dxx dcc min max mean hi ci
                    self.VOI_names.append(voiname)
                    self.VOI_pres_Dose.append(voidose)
                    if self.lowerdoseforext>float(voidose):
                        self.lowerdoseforext=float(voidose)
                    self.VOI_volumes.append(str(dIrrVolcc))
                self.VOI_Parameter.append('min')
                self.VOI_Parameter.append('max')
                self.VOI_Parameter.append('mean')
                self.VOI_Parameter.append('CI')
                self.VOI_Parameter.append('HI')
                for j in voiVxx:
                    self.VOI_Parameter.append('V'+str(j)+'%')
                for m in voiDxx:
                    self.VOI_Parameter.append('D'+str(m)+'%')
                for n in voiDcc:
                    self.VOI_Parameter.append('D'+str(n)+'cc')    
            elif voitype=='OAR':
                for i in range (0, len(voiDcc)+1): # dcc mean
                    self.VOI_names.append(voiname)
                    #print('voidose',voidose)
                    self.VOI_pres_Dose.append('/'+str(self.fractions)+'Fxs/')
                    self.VOI_volumes.append(str(dIrrVolcc))
                self.VOI_Parameter.append('mean')
                for n in voiDcc:
                    self.VOI_Parameter.append('D'+str(n)+'cc')
            # for external
            elif voitype=='EXT':
                for i in range (0, len(voiVxx)): # vxx dxx dcc min max mean hi ci
                    self.VOI_names.append(voiname)
                    self.VOI_pres_Dose.append(self.lowerdoseforext)
                    self.VOI_volumes.append(str(dIrrVolcc))
                for j in voiVxx:
                    self.VOI_Parameter.append('V'+str(j)+'%')
            else:
                wirteloginfo='a wrong type found, check voitype'
                related_funs.writelog(self.path2log, wirteloginfo)

            
        self.writelinesinfo.append(self.VOI_names)
        self.writelinesinfo.append(self.VOI_volumes)
        self.writelinesinfo.append(self.VOI_pres_Dose)
        self.writelinesinfo.append(self.VOI_Parameter)
        
        # get the voi corresponding data
        for filename in self.FileList:
            if 'NO_CT' in filename:
                break
            self.VOI_data.append([]) 
            print("File: "+filename)
            # get External V95 info for CI calculation for target.
            for voiname, voitype, voidose, voiVxx, voiDxx, voiDcc in self.voilist:
                if voitype=='EXT':
                    self.ExtV95=self.getExternalV95(filename, voiname, voitype, self.lowerdoseforext)
            for voiname,voitype,voidose,voiVxx,voiDxx,voiDcc in self.voilist:
                Dmin,Dmax,Dmean,CI,HI,Vxxlist,Dxxlist,Dcclist=self.getDVHMetricsFromFileByVOI(filename,voiname,voitype,voidose,voiVxx,voiDxx,voiDcc)
                if voitype=='Target':
                    self.VOI_data[-1].append(str('%.2f' % Dmin))
                    self.VOI_data[-1].append(str('%.2f' % Dmax))
                    self.VOI_data[-1].append(str('%.2f' % Dmean))
                    self.VOI_data[-1].append(str('%.2f' % CI))
                    self.VOI_data[-1].append(str('%.2f' % HI))
                    for Vxxinfo in Vxxlist:
                        self.VOI_data[-1].append(str('%.2f' % Vxxinfo + '%'))
                    for Dxxinfo in Dxxlist:
                        self.VOI_data[-1].append(str('%.2f' % Dxxinfo + '%'))
                    for Dccinfo in Dcclist:
                        self.VOI_data[-1].append(str('%.2f' % Dccinfo))
                elif voitype=='OAR':
                    self.VOI_data[-1].append(str('%.2f' % Dmean))
                    for Dccinfo in Dcclist:
                        self.VOI_data[-1].append(str('%.2f' % Dccinfo))
                elif voitype=='EXT':
                    for Vxxinfo in Vxxlist:
                        self.VOI_data[-1].append(str('%.2f' % Vxxinfo + '%'))
        for voidata in self.VOI_data:
            self.writelinesinfo.append(voidata)
        reverseinfo=list(zip(*self.writelinesinfo))
        return(reverseinfo)

    def getDVHMetricsFromFileByVOI(self,filename,voiname,voitype,voidose,voiVxx,voiDxx,voiDcc):
        
        #print voiname
        print("VOI: "+voiname)
        
        if voidose!=self.PlanDose:
            dTarget_Fraction_Dose=float(voidose)/float(self.fractions)
            coefficient=float(self.PlanDose)/dTarget_Fraction_Dose
        else:
            coefficient=1
        
        bGotIrrVol = False
        bVOIinData = False
        dIrrVolcc = -1.0
        
        
        Dmin=-1
        Dmax=-1
        Dmean=-1
        CI=-1
        HI=-1
        Vxxlist=[]
        Dxxlist=[]
        Dcclist=[]
        
        xvalues = []
        yvalues = []

        
        bVOIinData = True
        bGotIrrVol = False
        dIrrVolcc = self.getIrrVolByVOI(filename,voiname)/1000 # mm^3 --> cc
        
        if dIrrVolcc>0.0:
            bGotIrrVol = True
            #print("Voi: "+voiname+" has volume "+str(dIrrVolcc)+" cc")
            xvalues, yvalues = self.getDVHValuesFromFileByVOI(filename,voiname)

        xvalues = np.array(xvalues)
        yvalues = np.array(yvalues)
        xvalues = xvalues*coefficient
        
        Dmin = self.Dmin(voiname,filename)*float(self.PlanDose)*float(self.fractions)
        Dmean = self.Dmean(voiname,filename)*float(self.PlanDose)*float(self.fractions)
        Dmax = self.Dmax(voiname,filename)*float(self.PlanDose)*float(self.fractions)
        
        DVal = 0
        if(voitype=="Target"):
            #finding D95 & D5
            #D5  = self.D_n(5.0,[yvalues],xvalues)[0] 
            #D95 = self.D_n(95,[yvalues],xvalues)[0]
            D5=self.fun_Dxx(5, xvalues, yvalues)
            D95=self.fun_Dxx(95, xvalues, yvalues)
            HI = D5 - D95
            #CI = self.getCIByVOI(filename,voiname)

            V95 = self.fun_Vxx(95, xvalues, yvalues)
            if dIrrVolcc != 0:
                CI = self.ExtV95/float(dIrrVolcc)
                #CN = float(V95 * V95 / (CI * 100.))
            else:
                print("Calculation of CI/CN not possible. CI=0 cases!")
                CI = 'NA'
                writeinfo = "Voi " + voiname + " CI = NA"
                related_funs.writelog(self.path2log, writeinfo)
        for temp in voiVxx:
            try:
                Vxxlist.append(self.fun_Vxx(temp,xvalues,yvalues))
            except:
                Vxxlist.append(-1)
        for temp in voiDxx:
            DVal=temp
            try:
                Dxxlist.append(self.fun_Dxx(temp,xvalues,yvalues))
            except:
                Dxxlist.append(-1)
        for temp in voiDcc:
            Dval=float(temp)/float(dIrrVolcc)*100.00
            print("Dval=",Dval)
            #try:
            Dcclist.append(self.fun_Dxx(Dval,xvalues,yvalues)*float(voidose)/100.00)
            #except:
            #    Dcclist.append(-1)
        #print(Dcclist)
        return Dmin,Dmax,Dmean,CI,HI,Vxxlist,Dxxlist,Dcclist

    def getExternalV95(self, filename, voiname, voitype, voidose):
        # print voiname
        print("External: " + voiname)

        if voidose != self.PlanDose:
            dTarget_Fraction_Dose = float(voidose) / float(self.fractions)
            coefficient = float(self.PlanDose) / dTarget_Fraction_Dose
        else:
            coefficient = 1
        xvalues = []
        yvalues = []
        dIrrVolcc = self.getIrrVolByVOI(filename, voiname) / 1000  # mm^3 --> cc
        if dIrrVolcc > 0.0:
            xvalues, yvalues = self.getDVHValuesFromFileByVOI(filename, voiname)
        xvalues = np.array(xvalues)
        yvalues = np.array(yvalues)
        xvalues = xvalues * coefficient
        ExtV95=0
        if (voitype != "EXT"):
            print('wrong input, this is for External structure only')
        else:
            percen=self.fun_Vxx(95, xvalues, yvalues)
            ExtV95 = percen/100*dIrrVolcc # abs v95%
        return ExtV95

    def getDVHValuesFromFileByVOI(self, filename, VOIstr,**kwargs):
        bReadDifferential = kwargs.get("bReadDifferential",False)
        """Read x and y values from GD file for DVH by VOIstr"""
        with open(filename,'r') as fin:
            startline=0
            stopline=0
            stoplines=[]
            j=0
            for line in fin:
                j+=1
                if re.match("c:[ ]*{}".format(VOIstr+' '), line):
                    startline=j+2
                if re.match("c:[ ]*VOI",line):
                    stoplines.append(j-1)
            for elem in stoplines:
                if elem > startline:
                    stopline = elem
                    break
            if stopline==0:
                stopline=j
            if startline==0:
                print("VOI %s not found!" % VOIstr)
                exit(1)
        with open(filename,'r') as fin:
            xvalues=[]
            yvalues=[]
            jj=0
            for line in fin:
                jj+=1
                if jj < startline or jj > stopline: continue
                values=line.split()
                if len(values)==3:
                    xvalues.append(float(values[0]))
                    if bReadDifferential:
                        yvalues.append(float(values[2]))
                    else:
                        yvalues.append(float(values[1]))
        			
        return (xvalues,yvalues) 

    def Dmin(self,voi,filename):
        with open(filename,'r') as fin:
            for line in fin:
                if line.startswith("c: "+ voi):
                    maxdose=float(line.split()[2])
            result=maxdose
        return result
    
    def Dmax(self,voi,filename):
        with open(filename,'r') as fin:
            for line in fin:
                if line.startswith("c: "+ voi):
                    maxdose=float(line.split()[3])
            result=maxdose
        return result
	
    def Dmean(self,voi, path):
        with open(path,'r') as fin:
            for line in fin:
                if line.startswith("c: "+ voi):
                    meandose=float(line.split()[4])
            result=meandose
        return result
	
    def getIrrVol(self,path):
        with open(path,'r') as fin:
            for line in fin:
                values=line.split()
                if len(values)!=3:
                    continue
                else:
                    irrvol=values[1]
                    break
            result=irrvol
        return result
	
    def getIrrVolByVOI(self,path, VOIstr):
        irrvol=-1
        with open(path,'r') as fin:
            for line in fin:
                if line.startswith("c: "+VOIstr):
                    irrvol=float(line.split()[7])
                    break
        return irrvol
	
    def getCIByVOI(self,path, VOIstr):
        CI=-1
        with open(path,'r') as fin:
            for line in fin:
                if line.startswith("c: "+VOIstr):
                    CI=float(line.split()[10])
                    break
        return CI
							
    def V(self,n, path):
        result=[]
        with open(path,'r') as fin:
            for line in fin:
                values=line.split()
                if len(values)!=3:
                    continue
                else: 
                    values[0]=int(values[0])
                    values[1]=float(values[1])
                if values[0]==n:
                    relevantDose=values[0]
                    relevantVolume=values[1]
                    result=relevantVolume
        return result

    def D(self,n, path):
        result=[]
        with open(path,'r') as fin:
            for line in fin:
                values=line.split()
                if len(values)!=3:
                    continue
                else: 
                    values[0]=float(values[0])
                    values[1]=float(values[1])
                if values[1]<n:
                    smallervaluevolume=values[1]
                    smallervaluedose=values[0]
                if values[1]>n:
                    biggervaluevolume=values[1]
                    biggervaluedose=values[0]
                    continue
                finalvalue=(smallervaluedose - biggervaluedose)/(biggervaluevolume - smallervaluevolume)*(n - smallervaluevolume)+smallervaluedose
                result=finalvalue#print str(biggervaluevolume) + " " + str(smallervaluevolume)
                #print str(biggervaluedose) + " " + str(smallervaluedose)
                #print str(finalvalue)
                break
        return result
    
    def D_n(self,nn,yDataList,xData):
        nn=float(nn)
        D_left=[]
        D_right=[]
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
                    D_right.append( [xData[j], e] )
                    break
            for j, e in reversed(list(enumerate(elem))):
                if e >= nn:
                    D_left.append( [xData[j], e] )
                    break

        #interpolation between D_left and D_right
        if len(D_right) != len(D_left):
            print("List length from D_left and D_right don't match.")
            exit(1)
        else:
            Dnn=[]
            for i in range(len(D_left)):
                if D_right[i][1] == D_left[i][1]:
                    Dnn.append(D_left[i][0])
                else:
                    m = ( D_right[i][1] - D_left[i][1] ) / ( D_right[i][0] - D_left[i][0] )
                    n = D_left[i][1] - m * D_left[i][0]
                    Dnn.append( (nn - n)/m )
            return Dnn    
        
    def V_n(self,nn,yDataList,xData):
        nn=float(nn)
        V_left=[]
        V_right=[]
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
                    V_right.append( [xData[j], e] )
                    break
            for j, e in reversed(list(enumerate(elem))):
                if e >= nn:
                    V_left.append( [xData[j], e] )
                    break

        #interpolation between D_left and D_right
        if len(V_right) != len(V_left):
            print("List length from D_left and D_right don't match.")
            exit(1)
        else:
            Vnn=[]
            for i in range(len(V_left)):
                if V_right[i][1] == V_left[i][1]:
                    Vnn.append(V_left[i][0])
                else:
                    m = ( V_right[i][1] - V_left[i][1] ) / ( V_right[i][0] - V_left[i][0] )
                    n = V_left[i][1] - m * V_left[i][0]
                    Vnn.append( (nn - n)/m )
            return Vnn    
    def fun_Vxx(self,xx_per,Dose_values,Volume_values):
        Vxx=-1
        for i in range(0, len(Dose_values)):
            if Dose_values[i]==xx_per:
                Vxx=Volume_values[i]
                break
            elif Dose_values[i]<xx_per:
                continue
            elif Dose_values[i]>xx_per:
                Vxx=self.fun_trend(float(xx_per),Dose_values[i-1],Dose_values[i],Volume_values[i-1],Volume_values[i])
                break
        if Vxx==-1:
            print('V',xx_per,'is not possible to calculate')
        return Vxx
    
    def fun_Dxx(self,xx_per,Dose_values,Volume_values):
        Dxx=-1
        for i in range(0, len(Volume_values)):
            if Volume_values[i]==xx_per:
                Dxx=Dose_values[i]
                break
            elif Volume_values[i]>xx_per:
                continue
            elif Volume_values[i]<=xx_per:
                Dxx=self.fun_trend(float(xx_per),Volume_values[i-1],Volume_values[i],Dose_values[i-1],Dose_values[i])
                break
        if Dxx==-1:
            print('D',xx_per,'is not possible to calculate')
        return Dxx
    
    def fun_trend(self,newX,X1,X2,Y1,Y2):
        newY=Y1+(Y2-Y1)*(newX-X1)/(X2-X1)
        return newY

    def fun_ana_gamma(self,path2dose1,path2dose2,gammacri):
        gammaana = gamma.class_gammaanalysis()
        criterialist,gammalist=gammaana.fun_gamma_analysis(additionalinfo='', dose1=path2dose1, dose2=path2dose2, dosediscrit=gammacri, cuoff=10,maxdose='global', interfra=5, maxgamma=1.1, fraction=1, saveresults='False',pronecase=False,moreinfo=False)
        return gammalist
