import numpy as np
import nrrd

# 读取第一个NRRD文件
filename1 = '/home/yurii/Sheng/patient_data/dosetest/1_Lung_P_10_10phys.nrrd'
datafile1,headfile1 = nrrd.read(filename1)

# 读取第二个NRRD文件
filename2 = '/home/yurii/Sheng/patient_data/dosetest/1_Lung_P_10_11phys.nrrd'
datafile2,headfile2 = nrrd.read(filename2)


# 将两个数据文件相加
result_data = datafile1 + datafile2

# 保存结果为一个新的NRRD文件
result_header = nrrd.read_header(filename1)
result_header['data file']='senario1.dos.gz'
nrrd.write('/home/yurii/Sheng/patient_data/dosetest/senario1.nhdr', result_data,result_header)

