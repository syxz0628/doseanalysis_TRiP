"""
This file contains several function to analyse GD files.
Anna Eichhorn (a.eichhorn@gsi.de)

==> Copy of Moritz Wolf (mo.wolf@gsi.de)
changes in 	- getDVHValuesFromFile
"""
#import HatchHack
import re
import sys
import numpy as np
import pylab
import os
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc,rcParams
from matplotlib.patches import Rectangle
from datetime import datetime
from scipy.interpolate import interp1d

class DVH:
    Colors = [
            "black","lightgray",				# 0
            "green","lightgreen",				# 1
            "blue","lightskyblue",			# 2
            "red","pink", 					# 3
            "orange","navajowhite", 			# 4
            "darkmagenta","plum",				# 5
            "dimgray","dimgray",				# 6
            "darkblue","darkblue",			# 7
            "darkgreen","darkgreen",			# 8
            "midnightblue","midnightblue",	# 9
            "cyan","cyan",	# 10
            "deepskyblue","deepskyblue",	# 11
    ]

    def __init__(self, FileList, PlanDose=1, SaveName="DVHAnalysis"):
        print("Note: a log file named "+SaveName+".log was automatically created and must be closed at the end of the program using DVH.Close")
        self.FileList = FileList
        self.SaveName = SaveName
        self.PlanDose = PlanDose
        self.logFile = open(SaveName+".log", "w")
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S\n")
        self.logFile.write("Current date and time: \n")
        self.logFile.write(dt_string)
        self.logFile.write("List of DVH files to analyze:\n")
        self.voilist = []
        for name in FileList: self.logFile.write(name+"\n")

    def SetVOI(self,voiname, voilabel, voicolor, voitype):
        self.voilist.append([voiname,voilabel,voicolor,voitype])

        
    def Plot(self, labels, **kwargs):
        if len(self.voilist)==0: 
            print("Please first set the voi's to plot")
            exit(1)
        if not len(labels)==len(self.FileList): 
            print("List of labels must be equal to number of files to analyize")
            exit(1)

        lsty = kwargs.get('lsty',["-"]*len(self.FileList))
        lw  = kwargs.get('lw',[1]*len(self.FileList))
        alpha = kwargs.get('alpha', [1]*len(self.FileList))
        xlims = kwargs.get('xlims', [0,130])
        ylims = kwargs.get('ylims', [0,105])
        ylabel = kwargs.get('ylabel', r"Volume [%]")
        xlabel = kwargs.get('xlabel',r"Dose [%]")
        fontScale = kwargs.get('fontScale', 1)
        tickScale = kwargs.get('tickScale', 1)
        title = kwargs.get('title', "")
        showPlot = kwargs.get('showPlot', True)
        zoomBox = kwargs.get('zoomBox',[])
        legendLoc = kwargs.get('legendLoc',0)
        useTex = kwargs.get('useTex',True)
        axisGrid = kwargs.get('axisGrid',True)
        bPMB = kwargs.get('bPMB', False)
        bRBE = kwargs.get('bRBE', False)
        plotHline = kwargs.get('plotHline', 0)
        plotVline =  kwargs.get('plotVline', 0)
        bReadDifferential = kwargs.get('bReadDifferential', False)
        plotformat = kwargs.get('format','pdf')
        fancyFont = kwargs.get('fancyFont',False)
        normalize = kwargs.get('normalize',1)
        ## steering general plot behavior 
        if fancyFont: rc('font',**{'family':'sans-serif','sans-serif': ['DejaVu Sans']}) 
        if useTex: rc({'text.usetex': True})

	### set some constant values
        fontSize=11*fontScale
        fontSizeLabel=12*fontScale
        majorTick=20
        minorTick=5

        ### build a nicer grid    
        fig = plt.figure(figsize=(7.4, 4.57), dpi=300)
        nCol = 1
        if len(zoomBox): nCol = 2
        ax = fig.add_subplot(1,nCol,1)

        #ax.minorticks_on()

        major_xticks = np.arange(xlims[0], xlims[1], majorTick*tickScale)
        minor_xticks = np.arange(xlims[0], xlims[1], minorTick*tickScale)
        ### Y ticks!!!
        major_yticks = np.arange(ylims[0], ylims[1], majorTick)
        minor_yticks = np.arange(ylims[0], ylims[1], minorTick)

        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.set_yticks(major_yticks)
        ax.set_yticks(minor_yticks, minor=True)

        if axisGrid:
            ax.set_axisbelow(True)
            ax.set_facecolor('#E6E6E6')
            ax.grid(color='w',linestyle='solid')

        #plt.rc('axes',lw=3)

        LegendColors = []
        LegendLabels = []
        ##Makes nice legend
        for label,sty,w,a in zip(labels,lsty,lw,alpha):
            l = ax.plot([998,999],[98,99],ls=sty,color='k', alpha=a,lw=w) 
            LegendColors.append(l[0])
            LegendLabels.append(label)

        ax.set_xlim(xlims[0],xlims[1]) # X-Werte max
        ax.set_ylim(ylims[0],ylims[1])
		
        if bPMB: #ATTENTION for PMB only
            ax.set_xlabel("RBE-weighted dose [Gy]",size=fontSize)
        elif bRBE:
            ax.set_xlabel("RBE factor [a.u.]",size=fontSize)			
        else:
            ax.set_xlabel(xlabel,size=fontSize)

        ax.set_ylabel(ylabel,size=fontSize)
        if len(zoomBox): axins = fig.add_subplot(1,nCol,2,sharey=ax) #ax.inset_axes([0.5, 0.5, 0.47, 0.47])

        for voiname, voilabel,voicolor,voitype in self.voilist: 
            addedVOItoLegend = False
            for en,filename in enumerate(self.FileList):
                xvalues, yvalues = self.getDVHValuesFromFileByVOI(filename, voiname, bReadDifferential=bReadDifferential)
                if normalize!=1: xvalues = [x*normalize for x in xvalues]		
                if bPMB: #ATTENTION for PMB only
                    xvalues = [x*self.PlanDose/100.0 for x in xvalues]
                elif bRBE:
                    xvalues = [x/10.0 for x in xvalues]
						
                p1 = ax.plot( xvalues, yvalues,color=voicolor,ls=lsty[en], alpha=alpha[en], lw =lw[en])
                if len(zoomBox): axins.plot(xvalues,yvalues,color=voicolor,ls=lsty[en],alpha=alpha[en],lw=lw[en]) #zoomBox
                if not addedVOItoLegend:
                    aux = plt.plot([999,1000],[999,1000],color=voicolor,ls='-',linewidth=3.0)
                    LegendLabels.append(voilabel)
                    LegendColors.append(aux[0])
                    addedVOItoLegend = True

        if len(zoomBox):
            if axisGrid:
                axins.set_axisbelow(True)
                axins.set_facecolor('#E6E6E6')
                axins.grid(color='w',linestyle='solid')

            axins.set_xlim([zoomBox[0],zoomBox[1]])
            axins.set_ylim([zoomBox[2],zoomBox[3]])
            axins.set_xlabel(xlabel,size=fontSize)
            #ax.indicate_inset_zoom(axins, edgecolor="black")
            ax.add_patch(Rectangle((zoomBox[0], zoomBox[2]), zoomBox[1]-zoomBox[0], zoomBox[3]-zoomBox[2], facecolor="none", ec='red', lw=1))

        #plt.tight_layout()
        if title: plt.title(title, y=1.01) # size=fontSizeLabel,

        if plotHline:
            ax.hlines(y=95,xmin=0,xmax=130,color='grey',linestyle="-",linewidth=1.0)
        if plotVline:	
            if bPMB: #ATTENTION for PMB only
                ax.vlines(x=0.95*24.0,ymin=0,ymax=130,color='red',linestyle="-",linewidth=1.0)
            else:
                ax.vlines(x=95,ymin=0,ymax=130,color='grey',linestyle="-",linewidth=1.0)
        if legendLoc == 'out':
            if len(zoomBox):
                box = axins.get_position()
                axins.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                # Put a legend to the right of the current axis
                axins.legend(LegendColors, LegendLabels, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':fontSizeLabel},framealpha=0.8,fancybox=True)
            else:
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                # Put a legend to the right of the current axis
                ax.legend(LegendColors, LegendLabels, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':fontSizeLabel},framealpha=0.8,fancybox=True)
 
        else: ax.legend(LegendColors,LegendLabels,loc=legendLoc, prop={'size':fontSizeLabel},framealpha=0.8,fancybox=True) ### standard fuer FC plots

        fig.savefig(self.SaveName+"."+plotformat, bbox_inches='tight',dpi=300, transparent=False, format=plotformat)
        if showPlot: 
            plt.show()

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
	###debug
	#print startline,stopline
	
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


    
    def Close(self): 
        self.logFile.close()
#    
    def AnalyzeDVHsIndividual(self, **kwargs):
        for n in self.FileList: 
            for voiname,voilabel,voicolor,voitype in self.voilist:
                self.getDVHMetricsFromFileByVOI(n,voiname,voitype,**kwargs)

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
				
#def checkVOIStringInDVHfile(path, VOIstr):
#	with open(path,'r') as fin:
#		for line in fin:
#			### for debug only
#			if line.find(VOIstr+' ') >= 0:
#				#print line
#				return True
#			#if VOIstr in line:
#				#return True
#		return False
#				
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
#
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

#===============================================================================================#
##
##                               get average dvh values from file list
##                               add.: V(N), Homogeneityindex, D50, spread on D50
##
##===============================================================================================#
    def getDVHMetricsFromFileByVOI(self,name,voiname,voitype,**kwargs): #, bHI, dvDxx, bSpreadD50, bMaxDose, dMaxDoseLimit, sType, dTargetDose, dNormDose, f , bPrintToFile=False, bReadDifferential=False):

        Dxx = kwargs.get('Dxx',95)
        Vxx = kwargs.get('Vxx',95)
        bMaxDose = kwargs.get('bMAxDose',True)
        dMaxDoseLimit = kwargs.get('dMAxDoseLimit',self.PlanDose)
        normalize = kwargs.get('normalize',1)
        dTargetDose = self.PlanDose

        xvaluesMaxLen=0
        xvaluesMax = []
        yvaluesMaxLen=0
        bGotIrrVol = False
        bVOIinData = False
        dIrrVolcc = -1.0
        CI_data= []
        xvalues = []
        yvalues = []
        #print filenameList
        if not self.logFile.closed: self.logFile.write("File: "+name+"\n")
        print("File :"+name)
        bVOIinData = True
        if not bGotIrrVol:
            dIrrVolcc = self.getIrrVolByVOI(name,voiname)/1000 # mm^3 --> cc
        if dIrrVolcc>0.0:
            bGotIrrVol = True
            if not self.logFile.closed: self.logFile.write("Voi: "+voiname+" has volume "+str(dIrrVolcc)+" cc\n")
            print("Voi: "+voiname+" has volume "+str(dIrrVolcc)+" cc")
            xvalues, yvalues = self.getDVHValuesFromFileByVOI(name,voiname)
        else: 
            print("Voi "+voiname+" did not receive any dose.")
            if not self.logFile.closed: self.logFile.write("Voi "+voiname+" did not receive any dose.\n")

        xvalues = np.array(xvalues)*normalize
        yvalues = np.array(yvalues)
        

        if bVOIinData ==  False:
            print("Structure %s not found in given file list!" % VOIstr)
            return 0
        print("%s (Vol.: %.1f cc):\n" % (voiname,dIrrVolcc) )
        if not self.logFile.closed: self.logFile.write("%s (Vol.: %.1f cc):\n" % (voiname,dIrrVolcc) )
        Dmean = self.Dmean(voiname,name)*self.PlanDose
        if not self.logFile.closed: self.logFile.write("Dmean for voi "+voiname+" = "+str(Dmean))
        print("Dmean for voi "+voiname+" = "+str(Dmean))

        DVal = 0
        if(voitype=="target"):
    	    #finding D95 & D5
            D5  = self.D_n(5.0,[yvalues],xvalues)[0] 
            D95 = self.D_n(95,[yvalues],xvalues)[0]
            HI = D5 - D95		
            print(voiname+" Homogeneity index = %.1f" % HI)
            if not self.logFile.closed: self.logFile.write("Homogeneity index = "+str(HI)+"\n")
            CI = self.getCIByVOI(name,voiname)
            if not self.logFile.closed: self.logFile.write("CI of voi "+voiname+" of type "+voitype+" is: "+str(CI)+"\n")
            print("CI of voi "+voiname+" of type "+voitype+" is: "+str(CI))
            
            x_fine = np.linspace(xvalues.min(),xvalues.max(),num=10000, endpoint=True)
            y_interp1d = interp1d(xvalues, yvalues, kind='quadratic')
            y_fine = y_interp1d(x_fine)
            V95 = y_fine[np.argmin(np.abs(x_fine-95))]

            if CI: 
                CN=float(V95*V95/(CI*100.))
                print("CN of voi "+voiname+ " = " +str(CN))
                if not self.logFile.closed: self.logFile.write("CN of voi "+voiname+" = "+str(CN)+"\n")
            else:
                print("Calculation of CN not possible. CI=0 cases!")
            DVal = Dxx
        else:
            DVal = Dxx/dIrrVolcc*100

        DX = self.D_n(DVal,[yvalues],xvalues)[0]
        if voitype == "target":
            print(voiname+" D%.1f = %.2f" % (DVal,DX))
            if not self.logFile.closed: self.logFile.write(voiname+" D%.1f = %.2f" % (DVal,DX)+"\n")
        else:
            DxxNormalizedToConstraint = DX*self.PlanDose/Dxx
        if Vxx > 0.0:
            VX=yvalues[Vxx]
            print(voiname+" V"+str(Vxx)+" = %.2f\n" %VX)
            if not self.logFile.closed: self.logFile.write(voiname+" V"+str(Vxx)+" = %.2f\n" % VX)
		
        if bMaxDose:
            maxDose = xvalues[np.where(yvalues>0)][-1]

            print("target dose: %.1f" % dTargetDose)
            print("max. dose limit: %.1f" % dMaxDoseLimit)
            normfac=float(dTargetDose)/dMaxDoseLimit
            print("local normfac: %.2f" % normfac)
            #print maxDose
            print(voiname+" max Dose = %.1f" % (maxDose*normfac))
            if not self.logFile.closed: self.logFile.write(voiname+" max Dose = %.1f" % (maxDose*normfac)+"\n")
            print("\n-------------------------------------")
            if not self.logFile.closed: self.logFile.write("\n-------------------------------------\n\n")
        return 1
