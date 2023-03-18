import numpy as np

# 生成100个随机数
data=[95.9533,99.4857,94.812,97.977,99.7943,97.7815,99.9121,100.086,99.7453,97.2225,98.2949,100.0683,98.7291,96.4699,95.9228,97.8746,99.2943,96.4591,99.9688,100.0019,99.9627,95.2902,99.5243,99.4742,97.3819,99.3235,100.1076,96.0596,97.89,98.2235,100.1297,98.444,98.6283,98.9658,97.15,96.1378,99.4401,99.1133,96.0341,94.9394,99.4482,98.9821,97.621,96.0185,99.1293,99.9326,96.5422,100.0025,96.9471,98.2422
]
data2 = np.array(data)

# 依次计算第1到第2，第1到第3，第1到第4，...，第1到第100个数的中位数
for i in range(2, 50):
    median = np.median(data[:i])
    per5th = np.percentile(data[:i], 5)
    per95th = np.percentile(data[:i], 95)
    print(per95th)