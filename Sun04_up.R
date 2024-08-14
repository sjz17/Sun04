

###############【pqtl】###############
setwd("D:\\BaiduNetdiskDownload\\Sun04")
source("D:/BaiduNetdiskDownload/Sun00/function/code_function.R", encoding = 'GB18030')
get_library();library(QTLMR);library(data.table)
load("D:/BaiduNetdiskDownload/Sun04/fenland_3.Rdata")
ukb_ppp <- data_cispQTL_ukb_ppp_exp(clump = 2)  #decode cis-pQTL的暴露数据
decode <- data_cispQT_decode_exp(clump = 3)      #ukb-ppp cis-pQTL的暴露数据
fenland <- data_cispQT_fenland_exp()             #fenland cis-pQTL的暴露数据
fenland_3=rbindlist(fenland_3)
pqtl_data=list(ukb_ppp=ukb_ppp,decode=decode,fenland=fenland_3)
#i=2

out.mr.pqtl=lapply(c(1:3),function(i){
  data=pqtl_data[[i]]
  data=calculation_Fvalue_R2(data, samplesize = NULL)
  clump_data <- remove_MHC_data(dat = data,chr_col = "chr.exposure",pos_col = "pos.exposure",MHC_start = 28477797,MHC_end = 33448354)
  outcome_data <- extract_outcome_data(snps = clump_data$SNP,outcomes = "ebi-a-GCST005523")
  harm_data <- harmonise_data(clump_data,outcome_data)
  res_hete <- mr_heterogeneity(harm_data)
  res_pleioTab=mr_pleiotropy_test(harm_data)
  out1=xQTL_mr(harm_data, FDR_method = "fdr", PVE = TRUE)
  out1=generate_odds_ratios(out1)
  out2=subset(out1,out1$pval<0.05)
  out3=subset(out1,out1$fdr<0.05)
  out=list(clump_data=clump_data,harm_data=harm_data,
           res_hete=res_hete,res_pleioTab=res_pleioTab
           ,out1=out1,out2=out2,out3=out3)
})
names(out.mr.pqtl)=names(pqtl_data)
save(out.mr.pqtl,file = "D:\\BaiduNetdiskDownload\\Sun04\\out.mr.pqtl.Rdata")

out.mr.pqtl_dup=lapply(c(1:3),function(i){
  data=pqtl_data[[i]]
  data=calculation_Fvalue_R2(data, samplesize = NULL)
  clump_data <- remove_MHC_data(dat = data,chr_col = "chr.exposure",pos_col = "pos.exposure",MHC_start = 28477797,MHC_end = 33448354)
  outcome_data <- extract_outcome_data(snps = clump_data$SNP,outcomes = "ieu-a-1059")
  harm_data <- harmonise_data(clump_data,outcome_data)
  res_hete <- mr_heterogeneity(harm_data)
  res_pleioTab=mr_pleiotropy_test(harm_data)
  out1=xQTL_mr(harm_data, FDR_method = "fdr", PVE = TRUE)
  out1=generate_odds_ratios(out1)
  out2=subset(out1,out1$pval<0.05)
  out3=subset(out1,out1$fdr<0.05)
  out=list(clump_data=clump_data,res_hete=res_hete,res_pleioTab=res_pleioTab
           ,out1=out1,out2=out2,out3=out3)
})
names(out.mr.pqtl_dup)=names(pqtl_data)
save(out.mr.pqtl_dup,file = "D:\\BaiduNetdiskDownload\\Sun04\\out.mr.pqtl_dup.Rdata")



out.mr_pleioTab=lapply(out.mr.pqtl,function(data){out=data$res_pleioTab})
out.mr_pleioTab=rbindlist(out.mr_pleioTab,idcol =T )

out.mr_hete=lapply(out.mr.pqtl,function(data){out=data$res_hete})
out.mr_hete=rbindlist(out.mr_hete,idcol =T )

out.mr_F=lapply(out.mr.pqtl,function(data){out=data$clump_data;out=out[,c("SNP","exposure","Fvalue_Scheme4")]})
out.mr_F=rbindlist(out.mr_F,idcol =T )

write.csv(out.mr_pleioTab,"D:\\BaiduNetdiskDownload\\Sun04\\out.mr_pleioTab.csv", row.names = T, quote = T)
write.csv(out.mr_hete,"D:\\BaiduNetdiskDownload\\Sun04\\out.mr_hete.csv", row.names = T, quote = T)
write.csv(out.mr_F,"D:\\BaiduNetdiskDownload\\Sun04\\out.mr_F.csv", row.names = T, quote = T)


###########=
data.gene.all=lapply(out.mr.pqtl_dup,function(data){out=data$out1})
data.gene.all=do.call(rbind,data.gene.all)
data.gene.fdr=lapply(out.mr.pqtl_dup,function(data){out=data$out2})
data.gene.fdr=rbindlist(data.gene.fdr,idcol = T)
out.dup=subset(data.gene.all,data.gene.all$exposure %in% uni)


uni=data.gene.fdr[!duplicated(data.gene.fdr$exposure),];uni=uni$exposure

write.csv(out.mr.pqtl_dup.list,"D:\\BaiduNetdiskDownload\\Sun04\\out.mr.pqtl_dup.csv", row.names = T, quote = T)


######plot

load("D:/BaiduNetdiskDownload/Sun04/out.mr.pqtl.Rdata")
xQTL_volcano_plot(
  dat=out.mr.pqtl$decode$out1,
  breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
  scale_color = c("#686789","grey","#D89C7A"),
  FDR_method = "fdr",save_plot = TRUE,
  pdf_name = "volcano_plot_decode",
  width = 6,height = 5,save_path = "./"
)

#####filter

load("D:/BaiduNetdiskDownload/Sun04/out.mr.pqtl.Rdata")

gene=out.mr.pqtl$fenland$out3$exposure;gene
data.gene.all=lapply(out.mr.pqtl,function(data){out=data$out2})
data.gene.all=do.call(rbind,data.gene.all)
data.gene.fdr=lapply(out.mr.pqtl,function(data){out=data$out3})
data.gene.fdr=rbindlist(data.gene.fdr,idcol = T)
uni=data.gene.fdr[!duplicated(data.gene.fdr$exposure),];uni=uni$exposure

data.gene.fdr.neg=subset(data.gene.fdr,data.gene.fdr$b<0)
dput(data.gene.fdr.neg$exposure)
plot_forest=mr_plot_forest(data.gene.fdr,'CeD');plot_forest

write.csv(data.gene.fdr,"D:\\BaiduNetdiskDownload\\Sun04\\data.gene.fdr.csv", row.names = T, quote = T)


pdf(file="D:\\BaiduNetdiskDownload\\Sun04\\plot_forest.pdf",  width=15, height=15)
plot_forest
dev.off();dev.off()



###############【SMR】###############
#options(scipen=999)
setwd("D:\\BaiduNetdiskDownload\\Sun04")
source("D:/BaiduNetdiskDownload/Sun00/function/code_function.R", encoding = 'GB18030')
library(QTLMR);library(tidyr)


mr.SMR.blood=function(uni){
  
  info=lapply(dput(uni),function(gene){ tem=Gene_start_end_info(gene)})
  info=rbindlist(info);
  
  out.SMR.eqtlgen=SMR_multi_HEIDI(
    GWAS_file = "./out.data/outcome/GCST005523_SMR.ma",
    GWAS_name = "",
    xQTL_file = "./out.data/eqtl/westra_eqtl_hg19",
    bfile_1000G = "C:/1kg.v3/EUR",
    Gene_name = dput(uni),
    save_path = "./out.data/SMR/",
    diff_freq = 0.2,diff_freq_prop = 0.05,MAF = 0.01,
    cis_wind = 1000,pval = 5e-08,smr_multi = TRUE,
    smr_multi_set_wind = NULL,smr_multi_ld_snp = 0.1,
    heidi_mtd = 1,HEIDI_test_pval = 0.00157,ld_upper = 0.9,ld_lower = 0.05,
    nSNPs_min = 3,nSNPs_max = 20, thread_num = 4)
  
  out.SMR.eqtlgen1=merge(out.SMR.eqtlgen,info,by.x='probeID',by.y='gene_ensembl')
  out.SMR.eqtlgen1$Gene=out.SMR.eqtlgen1$gene_symbol
  out.SMR.eqtlgen1=out.SMR.eqtlgen1[,-c(23:29)]
  out.SMR.eqtlgen=out.SMR.eqtlgen1
}

mr.SMR.eqtlgen=function(uni){
  
  info=lapply(dput(uni),function(gene){ tem=Gene_start_end_info(gene)})
  info=rbindlist(info);
  
  out.SMR.eqtlgen=SMR_multi_HEIDI(
    GWAS_file = "./out.data/outcome/GCST005523_SMR.ma",
    GWAS_name = "",
    xQTL_file = "./out.data/eqtl/CAGE.sparse.lite",
    bfile_1000G = "C:/1kg.v3/EUR",
    Gene_name = dput(uni),
    save_path = "./out.data/SMR/",
    diff_freq = 0.2,diff_freq_prop = 0.05,MAF = 0.01,
    cis_wind = 1000,pval = 5e-08,smr_multi = TRUE,
    smr_multi_set_wind = NULL,smr_multi_ld_snp = 0.1,
    heidi_mtd = 1,HEIDI_test_pval = 0.00157,ld_upper = 0.9,ld_lower = 0.05,
    nSNPs_min = 3,nSNPs_max = 20, thread_num = 4)
  
  out.SMR.eqtlgen1=merge(out.SMR.eqtlgen,info,by.x='probeID',by.y='gene_ensembl')
  out.SMR.eqtlgen1$Gene=out.SMR.eqtlgen1$gene_symbol
  out.SMR.eqtlgen1=out.SMR.eqtlgen1[,-c(23:29)]
  out.SMR.eqtlgen=out.SMR.eqtlgen1
}
out.SMR.eqtlgen=mr.SMR.eqtlgen(uni)
save(out.SMR.eqtlgen,file = "D:\\BaiduNetdiskDownload\\Sun04\\out.SMR.eqtlgen.Rdata") 

mr.SMR.GTEx=function(uni){
  # file=list.files('D:\\BaiduNetdiskDownload\\Sun04\\out.data\\eqtl')
  # write.csv(file, file = 'D:\\BaiduNetdiskDownload\\file.csv', row.names = T) 
  # file=sapply(strsplit(file,"\\."),"[",1) 
  # save(file.list,file='D:\\BaiduNetdiskDownload\\Sun04\\file.list.Rdata')
  # file.list=paste0('./out.data/eqtl/',file.list)
  
  file=sapply(strsplit(file.list,"\\/"),"[",4)
  file=sapply(strsplit(file,"\\."),"[",1)
  
  file.list1=file.list[1]
  out.SMR.GTEx=lapply(file.list,function(file.name){
    SMR_multi_HEIDI(
      GWAS_file = "./out.data/outcome/GCST005523_SMR.ma",
      GWAS_name = "",
      xQTL_file = file.name,
      bfile_1000G = "C:/1kg.v3/EUR",
      Gene_name = dput(uni),
      save_path = "./out.data/SMR/",
      diff_freq = 0.2,diff_freq_prop = 0.05,MAF = 0.01,
      cis_wind = 1000,pval = 5e-08,smr_multi = TRUE,
      smr_multi_set_wind = NULL,smr_multi_ld_snp = 0.1,
      heidi_mtd = 1,HEIDI_test_pval = 0.00157,ld_upper = 0.9,ld_lower = 0.05,
      nSNPs_min = 3,nSNPs_max = 20, thread_num = 4
    )
  })
  names(out.SMR.GTEx)=file
}
out.SMR.GTEx=mr.SMR.GTEx(uni)
save(out.SMR.GTEx,file = "D:\\BaiduNetdiskDownload\\Sun04\\out.SMR.GTEx.Rdata")  


###########plot
load("D:/BaiduNetdiskDownload/Sun04/out.SMR.GTEx.Rdata")
load("D:/BaiduNetdiskDownload/Sun04/out.SMR.eqtlgen.Rdata")

#####=
out.SMR.GTEx$eqtlgen=out.SMR.eqtlgen
gene=data.frame(Gene=uni);tissue=data.frame('.id'=names(out.SMR.GTEx))
out.SMR.GTEx=rbindlist(out.SMR.GTEx,idcol=T)


#####=
out1=as.data.frame(out.SMR.GTEx[,c(1,4,18)]) 
merge1=full_join(gene,out1,by='Gene')
merge2=full_join(tissue,merge1,by='.id')
r1 <- pivot_wider(merge2, names_from = ".id", values_from = "b_SMR")

#####=
out2=as.data.frame(out.SMR.GTEx[,c(1,4,20,22)])
out2[is.na(out2)]=1
out2$p_SMR=ifelse(out2$p_HEIDI>0.05,out2$p_SMR,NA)
out2=out2[,-4]
out2[is.na(out2)]=1
merge1=full_join(gene,out2,by='Gene')
merge2=full_join(tissue,merge1,by='.id')
p1 <- pivot_wider(merge2, names_from = ".id", values_from = "p_SMR")

r1=bio_first(r1);p1=bio_first(p1)
r1=r1[-14,]
p1=p1[-14,]

r1=bio_data_process(r1);p1=bio_data_process(p1);

#####data
names=rownames(r1)
data1=ifelse(r1<0,'-','+')
out.data1=lapply(c(1:22),function(i){out=as.data.frame(table(data1[i,]))})
names(out.data1)=names
out.data1=rbindlist(out.data1,idcol = T,fill = T)
out.SMR=pivot_wider(out.data1,names_from =Var1,values_from = Freq)
write.csv(out.SMR,"D:\\BaiduNetdiskDownload\\Sun04\\out.SMR.csv", row.names = T, quote = T)

#####plot
#提取pvalue值矩阵；
p1=round(p1,3)

#预览转置后的相关性系数矩阵和pvalue矩阵；
r2 <- t(r1)
p2 <- t(p1)


#使用显著性星号标记进行替换；

#p2[p2>=0 & p2 < 0.001] <- "***"
#p2[p2>=0.001 & p2 < 0.01] <- "**"
p2[p2 < 0.05/22] <- "**"
p2[p2>=0.05/22 & p2 < 0.05] <- "*"
p2[p2>=0.05 & p2 <= 1] <- "-"


p2[is.na(p2)]=""
r2[is.na(r2)]=0
r2=as.matrix(t(r2))
p2=as.matrix(t(p2))

r2=bio_data_process(r2);p2=bio_data_process(p2);
#载入pheatmap包；
library(pheatmap)
#自定义颜色；
mycol = c(colorRampPalette(c("#686789", "white"))(100),
          colorRampPalette(c("white", "#D89C7A"))(100))

#绘制热图；
plot=pheatmap::pheatmap(r2,scale = "none",
                        border_color ="grey",
                        number_color='black',
                        legend = T,   #标签
                        fontsize_number=14,
                        fontsize_row=8,
                        fontsize_col=9,
                        cellwidth=15,
                        cellheight=15,
                        cluster_rows=F,
                        cluster_cols=F,
                        #gaps_row =c(length(colnames(GSE146853.mult$data.tryptophan))),
                        color = mycol,
                        breaks=c(seq(-2,2,length.out=200)),
                        #breaks = c(-0.2,0.2),
                        display_numbers = p2,
                        show_rownames=T,
                        shoe_colnames=T,
                        angle_col = 315,  #角度
                        annotation_legend = FALSE)

plot
pdf(file="D:\\BaiduNetdiskDownload\\Sun04\\SMR.pdf",  width=15, height=15)
plot
dev.off();dev.off()

#ieu-a-983


###############【coloc】###############
setwd("D:\\BaiduNetdiskDownload\\Sun04")
source("D:/BaiduNetdiskDownload/Sun00/function/code_function.R", encoding = 'GB18030')
get_library();library(QTLMR)
path='D:\\BaiduNetdiskDownload\\Sun04'


mr_extrect_decode=function(uni){
  info_decode=data_info_decode_pQTL()
  int=intersect(info_decode$gene_name_Ensembl,uni)
  tem=subset(info_decode,info_decode$gene_name_Ensembl %in% int );
  gene_name=tem$gene_name_Ensembl;gene_ID=tem$ID
  
  info_pQTL=lapply(c(1:length(gene_name)),function(i){
    tem1=get_cispQTL_decode_ukb_Online(
      gene =gene_name[[i]],gene_ID=gene_ID[[i]] ,resource = "decode", #ukb_ppp
      cis_wind_kb = 1000,build = 38,save_name = "tem",save_path = path,SMR_data = F)
  })
}
mr_extrect_ukb_ppp=function(uni){
  info_ukb=data_info_ukb_ppp_pQTL()
  int=intersect(info_ukb$HGNC.symbol,uni)
  
  tem=subset(info_ukb,info_ukb$HGNC.symbol %in% int );
  gene_name=tem$HGNC.symbol;gene_ID=tem$ID
  
  tem1=lapply(c(1:length(gene_name)),function(i){
    
    tem2=get_cispQTL_decode_ukb_Online(
      gene =gene_name[[i]],gene_ID  =gene_name[[i]]
      ,resource = "ukb_ppp", #ukb_ppp
      cis_wind_kb = 1000,build = 38,save_name = "tem",save_path = path,SMR_data = F)
  })
  
  tem1
}
mr_extrect_fenland=function(uni){
  info_fenland=data_info_fenland_pQTL()
  int=intersect(info_fenland$gene_name,uni)
  tem=subset(info_fenland,info_fenland$gene_name %in% int );
  gene_name=tem$gene_name;gene_ID=tem$SomaScan_ID
  
  tem1=lapply(1:length(gene_name),function(i){
    tem2=get_cispQTL_fenland_Online(
      gene = gene_name[[i]],
      ID = gene_ID[[i]],
      cis_wind_kb = 1000,
      build = 37,
      save_name = " ",
      save_path = "./",
      SMR_data = F
    )
  })
}

data.gene.ukb_ppp=mr_extrect_ukb_ppp(out.mr.pqtl$ukb_ppp$out3$exposure)
data.gene.fenland=mr_extrect_fenland(out.mr.pqtl$fenland$out3$exposure)
data.gene.decode=mr_extrect_decode(out.mr.pqtl$decode$out3$exposure)

names(data.gene.decode)=out.mr.pqtl$decode$out3$exposure
names(data.gene.fenland)=out.mr.pqtl$fenland$out3$exposure
names(data.gene.ukb_ppp)=out.mr.pqtl$ukb_ppp$out3$exposure
save(data.gene.decode,data.gene.ukb_ppp,data.gene.fenland,file = ".\\out.data\\Coloc\\data.gene.Rdata")  

outcome_dat <- read.delim("D:/BaiduNetdiskDownload/Sun04/out.data/outcome/GCST005523_TwosampleMR.txt")

out.coloc.decode=lapply(data.gene.decode,function(tem){
  out.coloc=coloc_GWAS(tem,outcome_dat ,
                       type_exposure = "quant",
                       col_pvalues_exposure = "pval",
                       col_N_exposure = "N",
                       col_MAF_exposure = "MAF",
                       col_beta_exposure = "beta",
                       col_se_exposure = "se",
                       col_snp_exposure = "SNP",
                       col_chr_exposure = "CHR",
                       col_pos_exposure = "BP",
                       sd_exposure = NA,
                       type_outcome = "cc",
                       col_pvalues_outcome = "pval.outcome",
                       col_N_outcome = "samplesize.outcome",
                       col_MAF_outcome = NA,
                       col_beta_outcome = "beta.outcome",
                       col_se_outcome = "se.outcome",
                       col_snp_outcome = "SNP",
                       prevalence_outcome = NA,
                       title1 = "xQTL",
                       title2 = "GWAS",
                       build = 37,
                       save_locus = FALSE,
                       save_stacked_dat = FALSE,
                       plot_pdf = "locuscompare",
                       width = 8,
                       height = 5,
                       save_path = path
  )
})
names(out.coloc.decode)=names(data.gene.decode)

out.coloc.list=list(out.coloc.decode=out.coloc.decode,
                    out.coloc.fenland=out.coloc.fenland,
                    out.coloc.ukb_ppp=out.coloc.ukb_ppp)
save(out.coloc.list,file = ".\\out.coloc.list.Rdata")


out.coloc.summary=lapply(out.coloc.list,function(list){
  tem=lapply(names(list),function(gene){
    tem=as.data.frame(t(list[[gene]]$summary))
  })
  names(tem)=names(list)
  tem=rbindlist(tem,idcol = T)
})

out.coloc.summary=rbindlist(out.coloc.summary,idcol = T)

out.coloc.summary$direct=ifelse(out.coloc.summary$PP.H4.abf>=0.8,'Yes','No')
save(out.coloc.summary,file = ".\\out.coloc.summary.Rdata")



###############【MAGMA】###############

setwd('D:/BaiduNetdiskDownload/Sun04')
source("D:/BaiduNetdiskDownload/Sun00/function/code_function.R", encoding = 'GB18030')
library(QTLMR)
path='D:\\BaiduNetdiskDownload\\Sun04'
#5.1 使用MAGMA软件进行基因层面和基因集合层面的关联分析
gene_based_dat <- MAGMA_gene_based(
  GWAS_file = "./out.data/outcome/GCST005523_METAL.txt",  # 输入的GWAS文件路径
  bfile_1000G = "C:/MR/g1000_eur/g1000_eur",  # 1000G参考数据路径
  gene_loc = "C:/MR/ENSGv110.coding.genes.txt",  # 基因位置文件路径
  set_annot = "C:/MR/MSigDB_20231Hs_MAGMA.txt",  # 基因集合注释文件路径
  SNP_P_col = c(3, 10),  # SNP p值列的位置
  samplesize_col = "N",  # 样本量列的名称
  save_name = "CeD",  # 保存的名称
  save_path = "./out.data/MAGMA/"  # 保存路径
)
# #5.2 使用MAGMA进行组织特异性分析
# Tissue_specific_dat <- MAGMA_Tissue_specific(
#   genes_raw = "./6.MAGMA/Migraine.genes.raw",  # 原始基因数据路径
#   gene_covar = "./6.MAGMA/gtex_v8_ts_avg_log2TPM.txt",  # 组织特异性表达数据路径
#   save_name = "Migraine",  # 保存的名称
#   save_path = "./6.MAGMA/"  # 保存路径
# )
#5.3 对MAGMA基因层面分析结果进行基因注释及曼哈顿图绘制
MAGMA_genes <- MAGMA_genes_Manhattanplot(
  genes_out = "./out.data/MAGMA/CeD.genes.out.txt",  # 基因分析结果文件路径
  gene_loc = "C:/MR/ENSGv110.coding.genes.txt",  # 基因位置文件路径
  Manhtn_plot = TRUE,  # 是否绘制曼哈顿图
  threshold_sig = "bonferroni",  # 显著性阈值调整方法
  Manhtn_gene_sig = 10,  # 曼哈顿图的基因标记标签的数目
  signal_cex = 1,  # 显著性基因图形圆点的大小
  width = 9,  # 图形宽度
  height = 7,  # 图形高度
  save_name = "MAGMA_gsa"  # 保存的名称
)

# 筛选FDR值小于0.05的基因
out.MAGMA=subset(MAGMA_genes,MAGMA_genes$GENE %in% uni)
out.MAGMA$FDR=p.adjust(out.MAGMA$P,method="fdr")
save(MAGMA_genes,file = ".\\out.MAGMA_genes.Rdata")


write.csv(out.MAGMA,"D:\\BaiduNetdiskDownload\\Sun04\\tab.MAGMA.csv", row.names = T, quote = T)






###############【enrich】###############

library(BioEnricher);library(DESeq2);library(tidyverse);library(clusterProfiler);library(org.Hs.eg.db);library(BioEnricher)

out.enrich <- lzq_ORA.integrated(
  genes =uni ,
  background.genes = NULL,
  GO.ont = 'ALL',
  perform.WikiPathways = F,
  perform.Reactome = T,
  perform.MsigDB = F,
  MsigDB.category = 'ALL',
  perform.Cancer.Gene.Network = T,
  perform.disease.ontoloty = T,
  perform.DisGeNET = T,
  perform.CellMarker = T,
  perform.CMAP = T,
  min.Geneset.Size = 3
)
save(out.enrich,file="D:\\BaiduNetdiskDownload\\Sun04\\out.enrich.Rdata")
#********
load("D:/BaiduNetdiskDownload/Sun04/out.enrich.Rdata")
plot=lzq_ORA.barplot1(enrich.obj = out.enrich$KEGG,show.term.num = 20,colors = colorRampPalette(c('#D89C7A', "#F3EEEA",'#686789'))(100));plot
pdf(file="D:\\BaiduNetdiskDownload\\Sun04\\gene.kegg.pdf",  width=11, height=8)
plot
dev.off();dev.off()

enrichGO=out.enrich$GO@result
plot=bio_plot_go(enrichGO);plot
pdf(file="D:\\BaiduNetdiskDownload\\Sun04\\gene.go.pdf",  width=10, height=6)
plot$my_plot
dev.off();dev.off()



###############【bulk-seq】###############


library(data.table)
load("D:/BaiduNetdiskDownload/Sun00/GEO/CD.GSE134900.Rdata")
load("D:/BaiduNetdiskDownload/Sun00/GEO/CD.GSE131705.Rdata")

table(GSE131705$cli$group)
gene.list=intersect(uni,rownames(GSE134900$exp_unique));gene.list
data=GSE134900

out=lapply(gene.list,function(gene){
  merge=merge(as.data.frame(t(data$exp_unique[gene,])),data$cli,by=0)
  merge=bio_data_process(merge)
  #out2=bio_plot_box(merge,'group',gene,comparisons = list(c('healthy','celiac disease')));out2
  out1=wilcox.test(merge[,gene] ~ group, data = merge);out1
  out1=as.data.frame(out1$p.value)
  data1=subset(merge,merge$group=='celiac disease');res1=mean(data1[,gene])
  data2=subset(merge,merge$group=='healthy');res2=mean(data2[,gene])
  
  out1$direct=ifelse(res1>res2,'+','-')
  out1
})

names(out)=gene.list
out=rbindlist(out,idcol = T)
out$direct=ifelse(out$`out1$p.value`<0.05,out$direct,'×')

colnames(data.gene.fdr)
out1=data.gene.fdr[,c(4,7)]
out2=merge(out,out1,by.x='.id',by.y="exposure")
out3=out2[,c(1,3)]
out.GSE131705=out3[!duplicated(out3)]
colnames(out.GSE131705)=c('gene','GSE131705')
save(out.GSE131705,out.GSE134900,file = ".\\out.GEO.Rdata")


#########plot

library(data.table)
load("D:/BaiduNetdiskDownload/Sun00/GEO/CD.GSE134900.Rdata")
load("D:/BaiduNetdiskDownload/Sun00/GEO/CD.GSE131705.Rdata")
#########=
table(GSE131705$cli$group)
gene.list=intersect(uni,rownames(GSE131705$exp_unique));gene.list
data=GSE131705
plot.list1=lapply(gene.list,function(gene){
  merge=merge(as.data.frame(t(data$exp_unique[gene,])),data$cli,by=0)
  merge=bio_data_process(merge)
  p1=bio_plot_box(merge,'group',gene,comparisons = list(c('healthy','celiac disease')));p1
})

#########=
table(GSE134900$cli$group)
gene.list=intersect(uni,rownames(GSE134900$exp_unique));gene.list
data=GSE134900
plot.list2=lapply(gene.list,function(gene){
  merge=merge(as.data.frame(t(data$exp_unique[gene,])),data$cli,by=0)
  merge=bio_data_process(merge)
  p1=bio_plot_box(merge,'group',gene,comparisons = list(c('healthy','celiac disease')));p1
})


plot=ggarrange(plot.list1[[1]],plot.list1[[2]],plot.list1[[3]],plot.list1[[4]],plot.list1[[5]],
               plot.list1[[6]],plot.list1[[7]],plot.list1[[8]],plot.list1[[9]],plot.list1[[10]],
               plot.list1[[11]],plot.list1[[12]],plot.list1[[13]],plot.list1[[14]],plot.list1[[15]],
               plot.list1[[16]],plot.list1[[17]],plot.list1[[18]],plot.list1[[19]],plot.list1[[20]],
               plot.list2[[1]],plot.list2[[2]],plot.list2[[3]],plot.list2[[4]],plot.list2[[5]],
               plot.list2[[6]],plot.list2[[7]],plot.list2[[8]],plot.list2[[9]],plot.list2[[10]],
               plot.list2[[11]],plot.list2[[12]],plot.list2[[13]],plot.list2[[14]],plot.list2[[15]],
               plot.list2[[16]],plot.list2[[17]],plot.list2[[18]],plot.list2[[19]]
               ,nrow = 8,ncol = 5,heights = c(1,1,1))
plot
pdf(file="D:\\BaiduNetdiskDownload\\Sun04\\gene2.pdf",  width=18, height=35)
plot
dev.off();dev.off()

###############【MR-Phe】###############
phewas_GeneATLAS <- read_csv("phewas  GeneATLAS.csv")
data_exposure <- read_excel("D:/BaiduNetdiskDownload/Sun00/Sun_exposure_ao_2023.4.2.xlsx")

data=subset(data_exposure,data_exposure$author %in% 'Pan-UKB team')
#############=
setwd("D:\\BaiduNetdiskDownload\\Sun04")
source("D:/BaiduNetdiskDownload/Sun00/function/code_function.R", encoding = 'GB18030')
get_library();library(QTLMR);library(ggplot2);library(data.table)
path='D:\\BaiduNetdiskDownload\\Sun04'

out_phe_BTN2A1 <- read_csv("D:/BaiduNetdiskDownload/Sun04/out.phe.BTN2A1.csv")
out_phe_CTSH <- read_csv("D:/BaiduNetdiskDownload/Sun04/out.phe.CTSH.csv")
out_phe_IL18R1 <- read_csv("D:/BaiduNetdiskDownload/Sun04/out.phe.IL18R1.csv")
out_phe_PTPRC <- read_csv("D:/BaiduNetdiskDownload/Sun04/out.phe.PTPRC.csv")

out_phe=list(BTN2A1=out_phe_BTN2A1,CTSH=out_phe_CTSH,
             IL18R1=out_phe_IL18R1,PTPRC=out_phe_PTPRC)

out_phe=rbindlist(out_phe,idcol=T)
out_phe=out_phe_IL18R1

out_phe$NegativeLogP=-log10(out_phe$`P-value`)
out_phe=subset(out_phe,out_phe$NegativeLogP>0)
colnames(out_phe)[1]='gene';colnames(out_phe)
#-log10(0.05/4756)

out_phe=bio_data_process(out_phe)
plot<-ggplot(out_phe, aes(x=Domain, y=NegativeLogP)) +
  scale_fill_brewer(palette="Set1")+
  scale_y_continuous(expand = c(0, 0))+
  geom_point(aes(fill=gene, color=gene),position = position_dodge(1),show.legend = T)+
  theme_bw()+
  geom_hline(yintercept = -log10(0.05/4756),linetype=4,lwd=0.3,color='orange')+
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 45,hjust = 1),
        legend.title = element_text(size = 12),
        legend.position = 'top',
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        plot.title = element_text(size = 18,face = "bold",hjust = 0.5))
plot
pdf(file="D:\\BaiduNetdiskDownload\\Sun04\\MR_phe.pdf", width=15, height=7)
plot
dev.off();dev.off()


load("D:/BaiduNetdiskDownload/Sun04/out.mr.pqtl.Rdata")

data.snp=lapply(out.mr.pqtl,function(data){out=data$harm_data})
data.snp=rbindlist(data.snp,idcol = T,fill = T)
out.snp=subset(data.snp,data.snp$exposure %in% c('CTSH','IL18R1','PTPRC','BTN2A1'))
out.snp=out.snp[,c(1,2,39)]


p<-ggplot(df, aes(x=OTU, y=NegativeLogP, color=Group, size=RA, shape=Level)) +
  geom_point(alpha=.7)+
  scale_shape_manual(values=c(16,1))+ #选择形状，nosig对应空心圆，enrich对应实心圆
  scale_size_continuous(range = c(0.5,5),
                        breaks = c(0.05, 0.10,0.50,1.00,5.00),
                        labels = c('0.05', '0.10','0.50','1.00','5.00'))+ #设置散点尺寸
  guides(color="none",shape="none")#隐藏图例




###############【table】###############

library(dplyr);library(readr);library(data.table);library(tidyr)
load("D:/BaiduNetdiskDownload/Sun04/out.mr.pqtl.Rdata")

uni=union(union(out.mr.pqtl$fenland$out3$exposure,out.mr.pqtl$decode$out3$exposure),
          out.mr.pqtl$ukb_ppp$out3$exposure)

gene.list=list(fenland=as.data.frame(out.mr.pqtl$fenland$out3$exposure),
               decode=as.data.frame(out.mr.pqtl$decode$out3$exposure),
               ukb_ppp=as.data.frame(out.mr.pqtl$ukb_ppp$out3$exposure))
gene=rbindlist(gene.list,idcol = T)
colnames(gene)[2]='gene'
gene$name='positive'
table <- pivot_wider(gene, names_from = ".id", values_from = "name")
data_cispQT=data_cispQT_fenland_exp()
##########=
info_ukb=data_info_ukb_ppp_pQTL()
info_ukb=info_ukb[,c(3,10)];colnames(info_fenland)[1]='gene'
merge_ukb=merge(info_ukb,out.mr.pqtl$ukb_ppp$out3,by.x='HGNC.symbol',by.y='exposure')
colnames(merge_ukb)[2]='protein'

info_decode=data_info_decode_pQTL()

info_fenland=data_info_fenland_pQTL()
info_fenland=info_fenland[,c(2,13)];colnames(info_fenland)[1]='gene'

Gene_start_end_info_ensembl('OLFM2')
info_fenland

merge_fenland=left_join(table,info_fenland,by=c('gene'))
colnames(merge_ukb)[2]='protein'
tab_protein=merge_ukb[,c(1,2)]
colnames(tab_protein)[1]='gene'
##########=
save(table,file = "D:/BaiduNetdiskDownload/Sun04/table.Rdata") 
load("D:/BaiduNetdiskDownload/Sun04/table.Rdata")

##########=
tab_MR <- read_excel("D:/BaiduNetdiskDownload/Sun04/tab.MR.xlsx")
tab_SMR <- read_csv("D:/BaiduNetdiskDownload/Sun04/tab.SMR.csv")
tab_coloc <- read_csv("D:/BaiduNetdiskDownload/Sun04/tab.coloc.csv")
tab_MAGMA <- read_csv("D:/BaiduNetdiskDownload/Sun04/tab.MAGMA.csv")
load("D:/BaiduNetdiskDownload/Sun04/tab.GEO.Rdata")

table=full_join(table,tab_MR,by=c('gene'))
table=full_join(table,tab_SMR,by=c('gene'))
table=full_join(table,tab_coloc,by=c('gene'))
table=full_join(table,tab_MAGMA,by=c('gene'))
table=full_join(table,out.GSE131705,by=c('gene'))
table=full_join(table,out.GSE134900,by=c('gene'))
table=full_join(table,tab_protein,by=c('gene'))
write.csv(table,"D:\\BaiduNetdiskDownload\\Sun04\\table.csv", row.names = T, quote = T)









