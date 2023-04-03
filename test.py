import pandas as pd
import numpy as np

# 生成示例数据集
dfA = pd.DataFrame(np.random.normal(size=(5, 6)), columns=['A', 'B', 'C', 'D', 'E', 'F'])

# 计算每列的5th profile和95%上下置信度区间
per5th = dfA.iloc[:, 1:].quantile(0.05, axis=0)
ci95_low = dfA.iloc[:, 1:].apply(lambda x: np.percentile(x, 2.5), axis=0)
ci95_high = dfA.iloc[:, 1:].apply(lambda x: np.percentile(x, 97.5), axis=0)

# 将结果转换为DataFrame对象，并将它们与原始数据集按列连接
df_per5th = pd.DataFrame(per5th).T.add_prefix('Per5th_')
df_ci95_low = pd.DataFrame(ci95_low).T.add_prefix('Ci95%low_')
df_ci95_high = pd.DataFrame(ci95_high).T.add_prefix('Ci95%high_')

dfA = pd.concat([dfA.iloc[:, 0], df_per5th, df_ci95_low, df_ci95_high, dfA.iloc[:, 1:-1]], axis=1)
# 按照指定顺序对列进行排序
dfA = dfA.sort_index(key=lambda x: x.str.split('_').str[-1])
print(dfA)
