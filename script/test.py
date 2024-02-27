import scanpy as sc
import matplotlib.pyplot as plt
import sys,os
import warnings
import numpy as np
from glob import glob
import matplotlib.colors as mcolors
import pandas as pd
from collections import OrderedDict
warnings.filterwarnings("ignore")

##每个样本的预处理
def data_predeal(data,sample):
    data.var_names_make_unique() #去除重复的基因
    sc.pl.highest_expr_genes(data, n_top=20,save='exp.png') # 每一个基因在所有细胞中的平均表达量（这里计算了百分比含量）
    sc.pp.filter_cells(data, min_genes=200) # 每一个细胞至少表达200个基因，过滤表达的基因数目少于200个的细胞
    sc.pp.filter_cells(data, max_genes=8000) # 每一个细胞最多表达8000个基因，过滤表达的基因数目多于8000个的细胞
    data.var['mt'] = data.var_names.str.startswith('MT')
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    ##绘制线粒体占比情况图
    sc.pl.scatter(data,x='total_counts', y='pct_counts_mt',save='%s.png'%sample)
    sc.pl.scatter(data,x='total_counts', y='n_genes_by_counts',save='%s.png'%sample)
    #过滤线粒体占比在20%以上的
    data = data[data.obs.pct_counts_mt < 20,:]
    #提取基因数量小于8000的细胞大于200
    # data = data[(200 < data.obs.n_genes_by_counts)&(data.obs.n_genes_by_counts < 8000),:]
    return data

#数据标准化
def data_deal(data,list_gene):
    list_gene = [gene for gene in list_gene if gene in data.var.index]
    #正则化数据集
    sc.pp.normalize_total(data, target_sum=1e4) #这个函数对每个单细胞的原始基因表达计数进行规范化，使得每个细胞的所有基因表达计数加起来等于一个固定的目标总和（这里是 1e4 或 10000）。这样做是为了消除不同细胞之间的测序深度差异，使得不同细胞之间的基因表达量可以进行比较。
    sc.pp.log1p(data) #这个函数对数据应用 log1p 转换，即对每个数值应用 log(x+1) 函数。这是为了把长尾分布的基因表达数据转换为更接近正态分布的数据，便于后续的统计分析
    #标准化数据集
    # 中心化
    sc.pp.scale(data, zero_center=True, max_value=None)
    data.X = data.X.clip(-10,10) #将所有的数据剪裁在-10和10之间
    # sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.highly_variable_genes(data)
    #筛选出高突变基因
    data.var.highly_variable[list_gene]=True
    data = data[:, data.var.highly_variable]
    # 过滤掉没用的东西
    # sc.pp.regress_out(data, ['total_counts', 'pct_counts_mt'])
    sc.pp.pca(data)
    data.write('pbmc3k_norm.h5ad')
    return data

##提取所有含有Mtx的文件夹
def dir_get(base_dir):
    list_dirs = []
    for dirs,dir_ch, file in os.walk(base_dir):
        if dirs == '.':continue
        if 'matrix.mtx.gz' in file:list_dirs.append(dirs)
    return list_dirs

def read_mtx(base_dir):
    dir_matrix = dir_get(base_dir)
    list_total, list_key = [], []
    for i in dir_matrix:
        prefix = i.strip().split('/')[-1]
        data = data_predeal(sc.read_10x_mtx(i),prefix)
        list_total.append(data)
        list_key.append(prefix)
    total = sc.concat(list_total,label = 'batch', keys=list_key)
    total.write('pbmc3k.h5ad')
    return total

#批量读取h5格式的文件
def read_h5(base_dir):
    list_file = glob(base_dir+"*/filtered_feature_bc_matrix.h5")
    list_total, list_key = [], []
    for i in list_file:
        prefix = i.strip().split('/')[-2]
        data = sc.read_10x_h5(i)
        list_total.append(data)
        list_key.append(prefix)
    total = sc.concat(list_total,label = 'batch', keys=list_key)
    total.write('pbmc3k.h5ad')
    return total

#绘制富集
def mark_cluster_new(data,check='louvain'):
    sc.pl.umap(data, color=check,show=False)
    # sc.pl.tsne(data, color="leiden",show=False)
    leiden = data.obs[check]
    site_info = data.obsm['X_umap']
    cluster_list = leiden.cat.categories
    for i in cluster_list:
        x,y = site_info[leiden==i].mean(axis=0)
        plt.text(x,y,i,fontsize=12)
    ax = plt.gca()
    ax.spines['top'].set_visible(False) #取消显示上边框
    ax.spines['right'].set_visible(False)#取消显示右边框
    ax.set_title('')
    plt.savefig("figures/umap_%s_cluster.png"%check,bbox_inches = 'tight')

#通过wilcoxon秩和检验比较差异基因
def diff_gene(data,thread=10):
    sc.tl.rank_genes_groups(data, 'louvain', method='wilcoxon',n_jobs=thread)
    result = data.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    dict_total = OrderedDict()
    dict_total['cluster'] = []
    for i in range(len(result['names'])):
        for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']:
            if key not in dict_total:dict_total[key] = []
            for group in groups:
                dict_total[key].append(result[key][i][group])
                if key=='names':dict_total['cluster'].append(group)
    return pd.DataFrame(dict_total),groups

def diff_gene_pic(data,groups,cmap,out_dir,total_info):
    for group in groups:
        # list_gene = total_info[(total_info['cluster']==group)&(total_info['logfoldchanges']>3)&(total_info['pvals']<0.01)&(total_info['pvals_adj']<0.01)].sort_values('logfoldchanges',ascending=False)['names'][:20].to_list() #提取每个cluster中logfc>3以及p_valur<0.01的显著差异基因用于后续绘制小提琴图以及umap图
        list_gene = total_info[(total_info['cluster']==group)&(total_info['pvals']<0.01)&(total_info['pvals_adj']<0.01)].sort_values('logfoldchanges',ascending=False)['names'][:20].to_list() #提取每个cluster中p_valur<0.01的显著差异基因用于后续绘制小提琴图以及umap图
        if len(list_gene)==1:
            sc.pl.violin(data,keys = list_gene, groupby='louvain',rotation=45,show=False)
            fig = plt.gcf()
            ax = plt.gca()
            ylabel = ax.get_ylabel()
            ax.set_title(ylabel)
            ax.set_ylabel('')
            plt.tight_layout()
            plt.savefig("%s/violin%s_cluster.png"%(out_dir,group),bbox_inches = 'tight',dpi=300)
        elif len(list_gene)>1:
            draw_violin(data,list_gene,out_dir,group)
        else:
            return 
        plt.clf()
        sc.pl.umap(data,color=list_gene,cmap=cmap,show=False)
        plt.savefig("%s/umap%s_cluster.png"%(out_dir,group),bbox_inches = 'tight')

####绘制小提琴图用于展示差异基因表达情况
def draw_violin(data,list_gene,out_dir,cluster):
    n = len(list_gene)  
    num_rows = (n + 3) // 4  # 计算所需行数
    num_cols = min(n, 4)  # 每行最多显示 4 张图片
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, 14))  # 创建子图
    # 循环绘制图片
    for i in range(n):
        row = i // num_cols
        col = i % num_cols
        ax = axs[row, col] if num_rows > 1 else axs[col]  # 处理只有一行的情况
        # 绘制图片
        sc.pl.violin(data,ax=ax,keys = list_gene[i], groupby='louvain',rotation=45,ylabel='')  # 添加其他参数以绘制图片
        ax.set_title(list_gene[i])
    for ax in axs.flat[n:]:
        ax.axis('off')
    plt.tight_layout()  # 调整子图布局
    plt.savefig("%s/violin%s_cluster.png"%(out_dir,cluster),bbox_inches = 'tight',dpi=300)  # 显示图像

###绘制marker基因的图
def draw_gene(data,gene,out_dir,cmap):
    sc.pl.umap(data, color=gene,show=False,cmap=cmap)
    leiden = data.obs['louvain']
    site_info = data.obsm['X_umap']
    cluster_list = leiden.cat.categories
    for i in cluster_list:
        x,y = site_info[leiden==i].mean(axis=0)
        plt.text(x,y,i,fontsize=12)
    ax = plt.gca()
    ax.spines['top'].set_visible(False) #取消显示上边框
    ax.spines['right'].set_visible(False)#取消显示右边框
    ax.set_title('')
    plt.savefig("%s/umap%s_cluster.png"%(out_dir,gene),bbox_inches = 'tight')

def gene_cluster(data,cell_gene,groups,total_info):
    colors = ["lightgray", "darkblue"]
    cmap = mcolors.LinearSegmentedColormap.from_list("mycmap", colors)
    index_gene = data.var.index
    for key,value in cell_gene.items():
        list_gene = [gene for gene in value if gene in index_gene]
        if not list_gene:continue
        dirs = "figures/%s"%key
        if not os.path.exists(dirs):
            os.mkdir(dirs)
        sc.pl.umap(data, color=list_gene ,save='%s_all.png'%key,cmap=cmap)
        for i in list_gene:
            draw_gene(data,i,dirs,cmap)
    dirs = 'figures/diff_gene'
    if not os.path.exists(dirs):
        os.mkdir(dirs)
    diff_gene_pic(data,groups,cmap,dirs,total_info)

def gene_get(file):
    dict_gene = {}
    with open(file,'r') as f:
        for line in f:
            if line.startswith('type'):continue
            line = line.strip().split("\t")
            if line[0] not in dict_gene:dict_gene[line[0]] = []
            dict_gene[line[0]].append(line[1])
    return dict_gene

#合并数据集
# data = read_mtx(sys.argv[1])
data = read_h5(sys.argv[1])
#获取marker基因
dict_gene = gene_get(sys.argv[2])
#合并数据集处理
data = data_deal(data,sum(dict_gene.values(),[]))
# data = sc.read_h5ad(sys.argv[1])
# sc.external.pp.bbknn(data,batch_key='batch') ##去除批次效应
sc.pp.neighbors(data)
sc.tl.umap(data)
# sc.tl.leiden(data)
sc.tl.louvain(data)
total_info, groups = diff_gene(data)  #通过wilcoxon检验获取每个基因的p值、logfc值
total_info.to_csv('Allgenepvalue.csv',sep='\t')
gene_cluster(data,dict_gene,groups,total_info)
mark_cluster_new(data,check='louvain')
sc.pl.umap(data,color='batch',save='sample.png') ##绘制样本富集图
# data.write('result.h5ad')