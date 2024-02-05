binf=script
volc=$binf/volcanoplot.R
umap=$binf/RDS_TypePlot.R
barh=$binf/scprop.r
go=$binf/draw_bar.py
go_buble=$binf/pathway_Bubble_plot_y6.R
umap_gene=$binf/umap_plot_240110.R
diff_gene=$binf/scRNA_diffgene.R
gene_func=$binf/clusterProfiler.test.R


###单纯用来了解各个文件的格式
##umap绘制 需要输入RDS文件
Rscript $umap -r rds -t test/type.txt -c test/type_col.txt -o test/ -s test
##umap gene 绘制
Rscript $umap_gene -r rds -g Ighg1 -o test/umap

##barh图绘制
Rscript $barh -s test/sample.txt -t test/type.txt -o test/barh -h 6 -w 10 -T test/type_col.txt -S test/sample_col.txt

#差异基因分析
Rscript $diff_gene -r rds -c test/compare.txt -o test/diff_gene -g test/sample.txt -t test/type.txt

##火山图绘制
Rscript $volc -d test/Treg.mBrain_vs_tLung.txt -f mBrain_vs_tLung -o test/vocalo/

#功能富集
Rscript $gene_func -f test/Treg.mBrain_vs_tLung.txt -n Treg -o test/gene_func -a true -s org.Hs.eg.db,hsa,human -g 1 -t SYMBOL -d  KEGG -C 0.05 