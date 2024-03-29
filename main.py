# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 15:27:59 2022

@author: ysheng
"""
import argparse
import analysis_gd
import acc_dos_nrrd
import sys

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--patientID", required=False, help="IDcd . of the patient")
    parser.add_argument("-p","--planname", required=False, help="name of the specific plan")
    parser.add_argument("-t","--targetnamelist", required=False, help="Target name list")
    parser.add_argument("-d","--targetdoselist", required=False, help="Target dose list")
    parser.add_argument("-o","--oarnamelist", required=False, help="OAR name list")
    parser.add_argument("-e", "--external", required=False, help="External name")
    parser.add_argument("-f","--fractions", required=False, help="fractions for this plan", default=1.0)
    parser.add_argument("-ga", "--gamma", required=False, action='store_true',help="add gamma analysis to the result")
    # parser.add_argument("-s","--savename", required=False, help="txt file save to name")
    parser.add_argument("-sp", "--savepath", required=False, help="file path to save logs")
    parser.add_argument("-r", "--robustevaluation", required=False, help="flag if DVH files are robust evaluation files, 21 or 9, Default None")
    parser.add_argument("-n", "--nameofgdlist", required=False, help="name of the data, such as 3DRecon, 4DPerRSC,4DMBR,,3DReview_CTxx, 4DRePerRSC_CTxx")
    parser.add_argument("-g", "--path2gdlist", required=False, help="path to gd file")
    parser.add_argument("-rs", "--referenceSpecial", required=False, action='store_true', help="analysis Reference data only")
    parser.add_argument("-fd", "--fourDgdlist", required=False, help="path to 4D gd files for accumulation dose analysis")
    parser.add_argument("-op", "--OptMethod", required=False, help="IMPT or SFO")
    parser.add_argument("-dp", "--doseshowninplansgd", required=False, help="dose_shown_in_gd is the dose written in "
                                                                            "the exec file, for SPHIC momi cases this "
                                                                            "value was set to 3 for all plans.", default=3)
    parser.add_argument("-sa", "--showallresult", required=False, action='store_true',help="select to show results of only worst cases")
    parser.add_argument("-ac", "--accumulatedose", required=False, help="input path of doses nrrd files and path to new dose file, accumulate doses")
    parser.add_argument("-rr", "--randomrobustana", required=False,help="analysis radnom senerio based DVH, 1,30 or more, default None")
    parser.add_argument("-fa", "--fxacc", required=False,
                        help="analysis radnom senerio based DVH, accumulate each 1-N fractions, Default None, give fx value" )
    parser.add_argument("-re", "--referencegd", required=False, help="path to referece gd file")
    parser.add_argument("-rf", "--referencefraction", required=False, help="referece gd file fractions")
    parser.add_argument("-rp", "--referencepd", required=False, help="referece gd file fractions")
    parser.add_argument("-ex", "--exportpath", required=False, help="export gds to.txt files with  x ref y1 y2 ...")
    #parser.add_argument("-t", "--timeoffset", required=False, type=int, nargs='+',
    #                    help="Time offset in msec,to adjust results in ~250ms level that was added to system determined timeoffset value;multiple values are acceptable, e.g. -t 250 -250 100",
    #                    default=250)
    #parser.add_argument("-g", "--log", required=False, nargs='?',
    #                    help="write error/successed information to .log file")
    args = parser.parse_args()

# define mandatory parameters.
# call analysis_gd function
    if args.path2gdlist!=None:
        patientID = args.patientID
        planname = args.planname
        targetnamelist = args.targetnamelist.split(',')
        targetdoselist = args.targetdoselist.split(',')
        oarnamelist = args.oarnamelist.split(',')
        fractions = args.fractions
        savepath = args.savepath
        gammaEva = args.gamma
        externalname = args.external
        robustevaluation = args.robustevaluation
        path2gdlist = args.path2gdlist.split(',')
        nameofgdlist = args.nameofgdlist.split(',')
        referenceSpecial = args.referenceSpecial
        OptMethod = args.OptMethod
        Planneddose = args.doseshowninplansgd
        Showallresult = args.showallresult
        randomrobust=args.randomrobustana # 1,30
        fractionsacc =args.fxacc # 5,10,20,23
        referecegd=args.referencegd
        referencefraction=args.referencefraction
        referenceplanneddose=args.referencepd
        exportpath=args.exportpath

        analysis_gd_data=analysis_gd.class_analysis_gd(patientID,planname,OptMethod,targetnamelist,targetdoselist,
                                                   oarnamelist, externalname,fractions,savepath,gammaEva,
                                                   robustevaluation,path2gdlist,nameofgdlist,Planneddose,Showallresult,
                                                       randomrobust,fractionsacc,referecegd,referencefraction,referenceplanneddose)
        if exportpath!=None:
            analysis_gd_data.fun_analysis_export(exportpath)
        else:
            analysis_gd_data.fun_analysis_gd()


    if args.path2gdlist==None and args.accumulatedose!=None:
        path2doslist = args.accumulatedose.split(',')
        analysis_dos_data=acc_dos_nrrd.class_analysis_dos_nrrd(path2doslist)
        analysis_dos_data.fun_accumulatedose()