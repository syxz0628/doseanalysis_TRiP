import numpy as np
import nrrd

# 读取第一个NRRD文件
data1, header1 = nrrd.read('/home/yurii/Sheng/patient_data/02/dose/TRiP_ana/3D-dose/1_Lung_P_10_10phys.nrrd')
# 初始化累加结果
sum_data = data1.astype(np.float64)
data, header = nrrd.read('/home/yurii/Sheng/patient_data/02/dose/TRiP_ana/3D-dose/1_Lung_P_10_11phys.nrrd')
sum_data += data.astype(np.float64)
#
# # 循环读取剩余的NRRD文件
# for i in range(2, 6):
#     filename = f'file{i}.nrrd'
#     data, header = nrrd.read(filename)
#     # 将当前文件的矩阵数据加到累加结果中
#     sum_data += data.astype(np.float64)
#
# # 将累加结果保存为新的NRRD文件
nrrd.write('/home/yurii/Sheng/patient_data/02/dose/TRiP_ana/3D-dose/sum.nrrd', sum_data, header1)
