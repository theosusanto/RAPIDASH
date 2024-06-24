if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

expt_file_raw = read.table(file='data/combined_file.txt',sep='\t', header=TRUE, row.names=1)
expt_cpm = cpm(expt_file_raw)
row_keep = filterByExpr(expt_cpm)
expt_file = expt_file_raw[row_keep,]

expt_names = colnames(expt_file)
#expt_file = expt_file[,!grepl('rep2',expt_names)] #uncomment to exclude rep2
#expt_names = colnames(expt_file)
expt_parsed = strsplit(expt_names,'_')
repid = sapply(expt_parsed,function(t) str_replace(t[4], "rep",))
conditionid = sapply(expt_parsed,function(t) t[1])
riboid = sapply(expt_parsed,function(t) t[2])

#BiocManager::install("tRanslatome")
##loading the tRanslatome package
library(biomaRt)
library(tRanslatome)
##sets the comparions - C vs A/E
condition_id = "Nterm"
expt_input = as.matrix(expt_file)
translatome.analysis <- newTranslatomeDataset(expt_input,
                                              expt_names[grepl("wt_RNA_.*",expt_names)], # control, RNA
                                              #expt_names[grepl(paste("_",condition_id,".*_mRNA",sep=""),expt_names)], # condition, RNA
                                              expt_names[grepl("Nterm_RNA_.*",expt_names)], # condition, RNA
                                              expt_names[grepl("wt_RFP_.*",expt_names)], # control, RFP
                                              expt_names[grepl("Nterm_RFP_.*",expt_names)], # condition, RFP
                                              label.level= c("RNA", "RFP"),
                                              data.type='ngs',
                                              #data.type='array',
                                              label.condition=c("control", paste("condition",condition_id)))

##identification of DE in RFP/RNA conditions. edgeR and DEseq handle counts directly.
# for other methods such as limma, there needs to be an additional normalization step beforehand and data.type needs to be set to array
edgeR.DEGs <- computeDEGs(translatome.analysis,
                          method= "edgeR", mult.cor=TRUE, significance.threshold=0.1)
dev.new(width=10, height=4)
Histogram(edgeR.DEGs)

Scatterplot(edgeR.DEGs)

MAplot(edgeR.DEGs)

# below here wont work with ensembl gene ids

##enrichment analysis of the selected DEGs
#CCEnrichment <- GOEnrichment(edgeR.DEGs,ontology="all", classOfDEGs="both",
#                             test.method="classic", test.threshold = 0.05)
##performing a comparison of the biological themes enriched
##in the two levels of gene expression
#CCComparison <- GOComparison(CCEnrichment)

#Radar(CCEnrichment)

#find genes that are significantly changing at the RFP level.
edgeR.DEGs@DEGs.table[(edgeR.DEGs@DEGs.table[,"edgeR.pval.mtc(RFP)"] < 0.1),]
write.csv(edgeR.DEGs@DEGs.table[(edgeR.DEGs@DEGs.table[,"edgeR.pval.mtc(RFP)"] < 0.1),], file = "S146A_vs_control_significant0.1_RFP.csv")

newtab = edgeR.DEGs@DEGs.table[(edgeR.DEGs@DEGs.table[,"edgeR.pval.mtc(RFP)"] < 0.1),]
#Check individual genes and how they change across replicates. Based on Adele's code.
#We want to check whether the genes that are pulled up as significant are mainly driven by one replicate or not.
rnactrl = which(grepl("wt_RNA_.",expt_names))
rfpctrl = which(grepl("wt_RFP_.",expt_names))
rnacond = which(grepl("Nterm_RNA", expt_names))
rfpcond = which(grepl("Nterm_RFP",expt_names))

interleave <- function(a, b) {
  mlab <- min(length(a), length(b)) 
  seqmlab <- seq_len(mlab) 
  c(rbind(a[seqmlab], b[seqmlab]), a[-seqmlab], b[-seqmlab]) 
}

anyexpt = c(rnactrl, rnacond, rfpctrl, rfpcond)
expt_cond = sapply(conditionid[anyexpt],function(x){substr(x,1,1)})
expt_type = riboid[anyexpt]
expt_rep = repid[anyexpt]

require(ggplot2)

gene_name = rownames(newtab)[1]

for(gene_name in rownames(newtab)){
  df = data.frame(log_expr_value=as.double(log(expt_file[gene_name,anyexpt]+1)), expt_cond=expt_cond, expt_type=expt_type, expt_rep=expt_rep)
  qplot(data=df, x=log_expr_value, y=expt_rep, 
        group=expt_rep, geom=c('point','line'), color=expt_cond, main=gene_name)+
    facet_wrap(~expt_type)
  ggsave(paste(gene_name,"_diff.pdf",sep=""),device="pdf")
}

