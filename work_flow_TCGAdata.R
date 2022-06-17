###work flow of TCGA data###
###Qiuyi Tang###############
###2022-6-17################
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

  library(data.table)
  library(stringr)
  library(TCGAbiolinks)
  library(RNAseqStat2)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(EnvStats)
  library(ggpubr)
  library(dplyr)
  library(survival)
  library(survminer)
  library(patchwork)
##expression datasets of a cancer
projects <- getGDCprojects()
project <-  "TCGA-PRAD"
query <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts")
GDCdownload(query = query)
data <- GDCprepare(query = query)

countsdata <- SummarizedExperiment::assay(data)

saveRDS(countsdata,file = paste0("TCGA_RNA_data/",project,".rds"))

##RNAseq_analysis
group_list <- ifelse(as.numeric(str_sub(colnames(countsdata),14,15)) < 10,'tumor','normal')
table(group_list)


rownames(countsdata) <- as.data.frame(str_split(rownames(countsdata),"\\.",simplify = T))[,1]

data_i <- Create_DEGContainer(species = "Human",
                              dataType = "Counts",
                              idType = "ENSEMBL",
                              expMatrix = countsdata,
                              groupInfo = group_list,
                              caseGroup = "tumor")
data_o <- runALL(object = data_i,dir = "RNAseq_Result")

##clinical informarion datasets

clin_data <- GDCquery_clinic(project = project, type = "clinical")


##expression comparison of a specific gene__boxplot
gene_select_1 <- "GCNT4"
gene_select_2 <- "MCF2"
TPM_data <- SummarizedExperiment::assay(data,4)
rownames(TPM_data) <- as.data.frame(str_split(rownames(TPM_data),"\\.",simplify = T))[,1]
TPM_data <- toSYMBOL(TPM_data,"Human")


gene_data_1 <- TPM_data[which(rownames(TPM_data) == gene_select_1),]
gene_data_1 <- gather(gene_data_1,submitter_id,)
gene_data_1$sampletype <- ifelse(as.numeric(str_sub(gene_data_1$submitter_id,14,15)) < 10,'tumor','normal')
colnames(gene_data_1)[2] <- gene_select_1
gene_data_1$submitter_id <- str_sub(gene_data_1$submitter_id,1,12)

gene_data_2 <- TPM_data[which(rownames(TPM_data) == gene_select_2),]
gene_data_2 <- gather(gene_data_2,submitter_id,)
gene_data_2$sampletype <- ifelse(as.numeric(str_sub(gene_data_2$submitter_id,14,15)) < 10,'tumor','normal')
colnames(gene_data_2)[2] <- gene_select_2
gene_data_2$submitter_id <- str_sub(gene_data_2$submitter_id,1,12)

gene_data <- cbind(gene_data_1,gene_data_2)[,c(1,2,3,5)]

boxplot <- ggplot(gene_data_1,aes(sampletype,gene_data_1[,2]))+
  geom_boxplot(aes(color=sampletype))+
  theme(axis.text.x = element_text(angle = 65,vjust = 0.6))+
  labs(
       y=paste0("Expression of ",gene_select," (TPM) ")
    )+
  theme_classic2()+
  stat_compare_means(label = "p.format")+
  stat_n_text()
boxplot


scatterplot <- ggplot(data = gene_data,aes(gene_data[,2],gene_data[,4]))+
  geom_point(size = 3.0, shape = 16)+
  geom_smooth(method = "lm")+
  theme_classic2()+
  labs(y=paste0("Expression of ",gene_select_2," (TPM) "),
       x=paste0("Expression of ",gene_select_1," (TPM) "))
scatterplot


gene_grade_data <-inner_join(gene_data,clin_data,by="submitter_id")[,c("submitter_id",gene_select_1,gene_select_2,"primary_gleason_grade")]
gene_grade_data <- gather(gene_grade_data,key = "gene",value = "expr",all_of(gene_select_1),all_of(gene_select_2))
gene_grade_data$primary_gleason_grade <- as.factor(gene_grade_data$primary_gleason_grade)

g <- ggplot(gene_grade_data,aes(primary_gleason_grade,expr))+
  geom_boxplot(aes(color=primary_gleason_grade))+
  theme(axis.text.x = element_text(angle = 65,vjust = 0.6))+
  labs(       y="Expression Level")+
  theme_classic2()+
  facet_wrap(~gene)+
  stat_compare_means(label = "p.format")+
  stat_n_text()
g

##survival analysis

gene_grade_data <-inner_join(gene_data,clin_data,by="submitter_id")
gene_grade_data <- gene_grade_data[which(gene_grade_data$sampletype=="tumor"),]


df<-subset(gene_grade_data,select =c(submitter_id,vital_status,days_to_death,days_to_last_follow_up,GCNT4,MCF2))

df[df$vital_status=='Dead',]$vital_status <- 2
df[df$vital_status=='Alive',]$vital_status <- 1
df$vital_status <- as.numeric(df$vital_status)
df$time <- df$days_to_death
df$time[which(is.na(df$time))] <- df$days_to_last_follow_up[which(is.na(df$time))]

with(df,Surv(time,vital_status))

fit <- survfit(Surv(time,vital_status)~exp,data = df)

p1 <- ggsurvplot(fit,pval = TRUE,
           conf.int = F, risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           legend.title="Expression of GCNT4")

p1
df$exp <- ''
df[df$MCF2 >= mean(df$MCF2),]$exp <- "H"
df[df$MCF2 < mean(df$MCF2),]$exp <- "L"

with(df,Surv(time,vital_status))

fit <- survfit(Surv(time,vital_status)~exp,data = df)

p2 <- ggsurvplot(fit,pval = TRUE,
                 conf.int = F, risk.table = F, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "hv", # Specify median survival
                 ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("#E7B800", "#2E9FDF"),
                 legend.title="Expression of MCF2")
p2

splots <- list()
splots[[1]] <- p1
splots[[2]] <- p2
q <- arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 2, nrow = 1)

(boxplot|scatterplot)/g
