import os
import sys
import time
import numpy as np
def max_index(lst_int):
    index = []
    max_n = max(lst_int)
    for i in range(len(lst_int)):
        if lst_int[i] == max_n:
            index.append(i)
    return index  #返回一个列表

def writelog(path2log,writeinfo):
    if path2log != None:
        print(path2log,' --> ',writeinfo)
        with open(path2log, 'a+') as logfile:
            logfile.writelines(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())+os.linesep)
            logfile.writelines(writeinfo+os.linesep)
            
def fun_trend(newX,X1,X2,Y1,Y2):
    newY=0
#    newY=Y1+(Y2-Y1)*(newX-X1)/(X2-X1)
    print('newY=',newY)
    return newY

def lambda_abs_max(arr, axis=None, key=None, keepdims=False):
    if callable(key):
        idxs = np.argmax(key(arr), axis)
        if axis is not None:
            idxs = np.expand_dims(idxs, axis)
            result = np.take_along_axis(arr, idxs, axis)
            if not keepdims:
                result = np.squeeze(result, axis=axis)
            return result
        else:
            return arr.flatten()[idxs]
    else:
        return np.amax(arr, axis)