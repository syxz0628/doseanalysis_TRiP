# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 15:27:59 2022

@author: ysheng
"""
import argparse
import analysis_gd
import sys

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--patientID", required=True, help="IDcd . of the patient") 
    parser.add_argument("-p","--planname", required=True, help="name of the specific plan") 
    parser.add_argument("-t","--targetnamelist", required=True, help="Target name list") 
    parser.add_argument("-d","--targetdoselist", required=True, help="Target dose list") 
    parser.add_argument("-o","--oarnamelist", required=True, help="OAR name list")
    parser.add_argument("-e", "--external", required=True, help="External name")
    parser.add_argument("-f","--fractions", required=True, help="fractions for this plan")
    parser.add_argument("-ga", "--gamma", required=False, action='store_true',help="add gamma analysis to the result",default=False)
    # parser.add_argument("-s","--savename", required=False, help="txt file save to name")
    parser.add_argument("-sp", "--savepath", required=False, help="file path to save logs")
    parser.add_argument("-r", "--robustevaluation", required=False, help="flag if DVH files are robust evaluation files, 21 or 9, Default None")
    parser.add_argument("-n", "--nameofgdlist", required=True, help="name of the data, such as 3DRecon, 4DPerRSC,4DMBR,,3DReview_CTxx, 4DRePerRSC_CTxx")
    parser.add_argument("-g", "--path2gdlist", required=True, help="path to gd file")
    parser.add_argument("-rs", "--referenceSpecial", required=False, action='store_true', help="analysis Reference data only")
    parser.add_argument("-fd", "--fourDgdlist", required=False, help="path to 4D gd files for accumulation dose analysis")
    #parser.add_argument("-t", "--timeoffset", required=False, type=int, nargs='+',
    #                    help="Time offset in msec,to adjust results in ~250ms level that was added to system determined timeoffset value;multiple values are acceptable, e.g. -t 250 -250 100",
    #                    default=250)
    #parser.add_argument("-g", "--log", required=False, nargs='?',
    #                    help="write error/successed information to .log file")
    args = parser.parse_args()

# define mandatory parameters.
    patientID=args.patientID
    planname=args.planname
    targetnamelist=args.targetnamelist.split(',')
    targetdoselist=args.targetdoselist.split(',')
    oarnamelist=args.oarnamelist.split(',')
    fractions=args.fractions
    savepath=args.savepath
    gammaEva=args.gamma
    externalname=args.external
    robustevaluation=args.robustevaluation
    path2gdlist = args.path2gdlist.split(',')
    nameofgdlist = args.nameofgdlist.split(',')
    referenceSpecial=args.referenceSpecial

# call analysis_gd function

    analysis_gd_data=analysis_gd.class_analysis_gd(patientID,planname,targetnamelist,targetdoselist,oarnamelist,externalname,fractions,savepath,gammaEva,robustevaluation,path2gdlist,nameofgdlist)
    if referenceSpecial:
        analysis_gd_data.fun_analysis_refonly(dose_shown_in_gd=3)
    else:
        analysis_gd_data.fun_analysis_gd(dose_shown_in_gd=3)

