import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


def df_deal(file,num=10):
    df = pd.read_table(file,index_col=2)
    df = -np.log10(df.loc[:,['p.adjust']])
    df = df.sort_values('p.adjust',ascending=False)
    return df.iloc[:num,:]


up_df = df_deal(sys.argv[1])
up_df = up_df.sort_values('p.adjust')
down_df = -df_deal(sys.argv[2])

fig, ax = plt.subplots()

# 分别绘制正值和负值
ax.barh(down_df.index.to_list(), down_df['p.adjust'].values, color='blue',label='down')
ax.barh(up_df.index.to_list(), up_df['p.adjust'].values, color='red',label='up')

# 设置y轴字体大小
ax.tick_params(axis='y', labelsize=15)  # 你可以改变这里的数字，以达到你想要的字体大小
ax.tick_params(axis='x', labelsize=15)  # 你可以改变这里的数字，以达到你想要的字体大小

ax.set_xticks(list(ax.get_xticks()))
ax.set_xticklabels([str(abs(x)) for x in ax.get_xticks()])

# ax.legend(bbox_to_anchor=(1.05, 1),fontsize=20)
handles, labels = ax.get_legend_handles_labels()
# 反转句柄和标签的顺序
handles = handles[::-1]
labels = labels[::-1]

# 添加图例
ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), fontsize=20)

#取消显示框线
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)



ax.set_xlabel('-Log10(p.adjust)',fontsize=20)
# ax.set_ylabel('Category')
plt.title(sys.argv[3],fontsize=30)

plt.savefig('Treg_%s.GO.png'%sys.argv[3],bbox_inches='tight')