# Load required library
library(dplyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(stats)
library(ggpubr)
library(viridis)
library(scales)
library(gprofiler2)

######### scMM lmmfit functions
source("~/FLASH-MM/R/lmmfit.nt.R")
source("~/FLASH-MM/R/lmmfitSS.R")
source("~/FLASH-MM/R/lmmtest.R")
source("~/FLASH-MM/R/qqpvalue.R")
library(FLASHMM)

get_gprofiler_enrich <- function(markers, model_animal_name){
  gostres <- gost(query = markers,
                  ordered_query = TRUE, exclude_iea =TRUE, 
                  sources=c('GO:BP' ,'REAC', 'GO:MF', 'GO:CC', 'KEGG', 'CORUM', 'HP', 'WP'), #'TF',
                  organism = model_animal_name)
  return(gostres)
}

model_animal_name ='hsapiens'
num_genes = 200

########################################################
############## data exploration and visualization
########################################################
umap_coord = read.table('~/scLMM/LMM-scRNAseq/Data/kidney_atlas/UMAP.coords.tsv.gz')
cell_meta = read.csv('~/scLMM/LMM-scRNAseq/Data/kidney_atlas/meta.tsv', sep = '\t')

sum(cell_meta$Cell != umap_coord$V1)
cell_meta = cbind(cell_meta, umap_coord)

ggplot(cell_meta, aes(x=V2, y=V3, color=sex))+geom_point(alpha=0.4, size=2)+
  theme_classic()+scale_color_manual(values = c("orchid1", 'lightskyblue'))+xlab('UMAP 1')+ylab('UMAP 2')+
  theme(text = element_text(size=16),legend.title = element_blank())#
a_cell_type = 'CNT'
cell_meta$a_cell_type = ifelse(cell_meta$Cell_Types_Broad==a_cell_type, a_cell_type, '')

ggplot(cell_meta, aes(x=V2, y=V3, color=a_cell_type))+geom_point(alpha=0.4, size=1.5)+
  theme_classic()+scale_color_manual(values = c('grey85', "green4"))+xlab('UMAP 1')+ylab('UMAP 2')+
  theme(text = element_text(size=16),legend.title = element_blank())#


unique_cell_types <- unique(cell_meta$Cell_Types_Broad)  # Get unique cell types
new_palette <- setNames(
  muted(rainbow(length(unique_cell_types)), l = 94, c = 30),  # More muted colors
  unique_cell_types
)
new_palette["CNT"] <- "#6A3D9A"  # Assign bold purple to cTAL

# Plot using the new palette
ggplot(cell_meta, aes(x=V2, y=V3, color=Cell_Types_Broad)) +
  geom_point(alpha=1, size=2) +
  theme_classic() +
  scale_color_manual(values = new_palette) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  theme(text = element_text(size=17), legend.title = element_blank())


df = data.frame(table(cell_meta$Cell_Types_Broad))
df[order(df$Freq, decreasing = T),]
#### based on sample-type ##### 
cell_meta_counts <- ddply(cell_meta, .(cell_meta$sex, cell_meta$Cell_Types_Broad), nrow)
names(cell_meta_counts) <- c("Sex", "CellType", "Freq")
#names(c25) = cell_meta_counts$CellType
#cell_meta_counts$CellType= factor(cell_meta_counts$CellType, levels = as.character(c25) ) 

ggplot(data=cell_meta_counts, aes(x=CellType, y=Freq, fill=Sex)) +
  geom_bar(stat="identity",color='black')+theme_classic()+#+scale_fill_brewer(palette = "Blues")+
  ylab('Counts')+xlab('Clusters')+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'),
        legend.title = element_blank())+
  xlab('')

cell_meta_counts_split <- split( cell_meta_counts , f = cell_meta_counts$CellType )
cell_meta_counts_split_norm <- lapply(cell_meta_counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,cell_meta_counts_split_norm )
head(counts_norm)

ggplot(data=counts_norm, aes(x=CellType, y=Freq, fill=Sex)) +
  geom_bar(stat="identity",color='black',alpha=0.9)+theme_classic()+
  ylab('Fraction per cell type (%)')+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=15,angle=90,color='black'),
        legend.title = element_blank()) +  
  xlab('')


######### importing kidney data  
#load('LMM-scRNAseq-jan2024/Kidney_reanalysis_CC/data/kidney-counts-lmmfit.RData')
load("~/FLASH-MM/Results_data//kidney-counts-lmm.beta.RData")

##running time
rtlmm
##t-values
tvlmm <- t(fit$t)
##p-values
pvlmm <- t(fit$p)
dim(pvlmm)
sum(apply(is.na(pvlmm), 1, any))
felmm <- t(fit$coef)
slmm <- fit$theta

##LMM tests
#test <- lmmtest(fit) ##t-values
#tvlmm <- test[, grep("_t", colnames(test)), drop = F]


##p-values
#plmm <- test[, grep("_pvalue", colnames(test)), drop = F]
hist(pvlmm)
pvalues_df=data.frame(pvalues=as.vector(pvlmm))
# basic histogram
ggplot(pvalues_df, aes(x=pvalues))+
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9)+theme_classic()+
  xlab('p-values')+ylab('Counts')+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=15,color='black'),
        legend.title = element_blank()) 


pvlmm_adj = sapply(1:ncol(pvlmm), function(i) p.adjust(pvlmm[,i], method = "fdr"))
colnames(pvlmm_adj) = colnames(pvlmm)
colnames(felmm) == colnames(pvlmm)
sum(apply(is.na(pvlmm), 1, any))
dim(pvlmm)
dim(felmm)

#### counts extracted from data object shared via paper 
dirData = '~/scLMM/LMM-scRNAseq/Data/kidney_atlas/'
datafile <- paste0(dirData, "/Human_Kidney_data.rds")
data = readRDS(file = datafile)
coldata = data@meta.data
Idents(data) = data$Cell_Types_Broad


# Initialize a data frame to store DE gene counts
de_counts <- data.frame(cell_type = colnames(felmm), n_DE_genes = NA)

# Loop over each cell type column (e.g., "CNT:Male", "PC:Male", etc.)
for (col_name in colnames(felmm)) {
  coef_vec <- felmm[, col_name]
  p_adj_vec <- pvlmm_adj[, col_name]
  de_count <- sum(p_adj_vec < 0.05 & abs(coef_vec) > 0.5)
  de_counts[de_counts$cell_type == col_name, "n_DE_genes"] <- de_count
}

# Sort by number of DE genes
de_counts <- de_counts[order(-de_counts$n_DE_genes), ]

# Display top few rows
head(de_counts)

################################################################################
############## cell-type-specific sex-biased gene identification ###############
################################################################################

#col_name = 'Proximal.Tubule:Male'
col_name ='CNT:Male' #'cTAL:Male'
cell_type_name = 'CNT'#'cTAL'
P_VAL_THR = 0.05
# The coefficients of the interaction term are the log-FC between Male and Female within a cell-type. 
PT_male_df = data.frame(pvalue=pvlmm[,col_name],
                        pvalue_adj = pvlmm_adj[,col_name],
                        tvalue = tvlmm[,col_name],
                        #coefficients(interaction)=log-FC(Male/Female per cell-type 
                        coef = felmm[,col_name],
                        gene=rownames(pvlmm))

PT_male_df$pvalue_adj_log = -log(PT_male_df$pvalue_adj+1e-800)
head(PT_male_df)
PT_male_df$pvalue_log = -log(PT_male_df$pvalue+1e-800)

sum(PT_male_df$pvalue_adj < 0.05 & abs(PT_male_df$coef) > 0.5) # 200

# Categorize significance with direction
PT_male_df <- PT_male_df %>%
  mutate(significance_group = case_when(
    pvalue_adj < 0.05 & coef > 0.5 ~ "Male-biased",
    pvalue_adj < 0.05 & coef < -0.5 ~ "Female-biased",
    TRUE ~ "Not Significant"
  ))

# Top 10 male-biased and female-biased genes
top_male <- PT_male_df %>%
  filter(significance_group == "Male-biased") %>%
  arrange(desc(coef)) %>%
  head(10)

top_female <- PT_male_df %>%
  filter(significance_group == "Female-biased") %>%
  arrange(coef) %>%
  head(10)

top_genes <- bind_rows(top_male, top_female)

# Volcano plot
ggplot(PT_male_df, aes(x = coef, y = pvalue_log)) +
  geom_point(aes(color = significance_group), size = 2.5, alpha = 0.85) +
  scale_color_manual(
    values = c(
      "Male-biased" = "#0072B2",    # Blue
      "Female-biased" = "#CC79A7",    # Pink
      "Not Significant" = "grey70"
    )
  ) +
  labs(
    x = "LogFC",
    y = expression(-log[10]("P-value")),
    color = "Significance"
  ) +
  theme_minimal(base_size = 18) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_label_repel(
    data = top_genes,
    aes(label = gene),
    box.padding = 0.5,
    point.padding = 0.4,
    label.size = 0.2,
    size = 6,
    segment.color = "black",
    segment.size = 0.6,
    min.segment.length = 0
  )


hist(PT_male_df$tvalue+1e-10)
hist(PT_male_df$coef)
hist(PT_male_df$pvalue)
hist(PT_male_df$pvalue_adj)

PT_male_df[is.na(PT_male_df$pvalue_adj),]
#data_sub = data[,data$Cell_Types_Broad %in% cell_type_name]


#write.csv(PT_male_df,"~/FLASHMM-analysis/SourceData/Figure3/Figure3C.csv",row.names = FALSE)

############################################################
for(i in 1:10){
  gene_name = PT_male_df_ord$gene[i]
  gene_name
  #hist(counts[gene_name,], main = gene_name)
  df = data.frame(gene=counts[gene_name,], status=coldata$sex)
  head(df)
  df = df[coldata$Cell_Types_Broad==names(table(coldata$Cell_Types_Broad))[15],]
  #ggplot(df, aes(x=status, y=gene))+geom_boxplot()+ggtitle(paste0(gene_name,' ' ,col_name))
  p=ggplot(df, aes(x=gene, color=status))+geom_density()+ggtitle(paste0(gene_name,' ' ,col_name))
  print(p)
}

PT_male_df$score = -log(PT_male_df$pvalue_adj+1e-600)*PT_male_df$coef
PT_male_df_ord = PT_male_df[order(PT_male_df$score, decreasing = TRUE),]
PT_male_df_sort = PT_male_df[order(PT_male_df$score, decreasing = T),]


PT_female_df_sort = PT_male_df[order(PT_male_df$score, decreasing = F),]

################################
#### Pathway analysis 
############################
num_genes=200
enrich_res_male = get_gprofiler_enrich(markers=PT_male_df_sort$gene[1:num_genes], model_animal_name)

enrich_res = enrich_res_male# enrich_res_female

head(enrich_res$result,30)
enrich_res_pos = data.frame(enrich_res$result)
enrich_res_pos = enrich_res_pos[1:20,]
enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'p_value')]
enrich_res_pos$log_p = -log(enrich_res_pos$p_value)
ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
  theme_classic()+ylab('')#+ggtitle(paste0(title))


################################
enrich_res_pos = data.frame(enrich_res$result)
enrich_res_pos = enrich_res_pos[,!colnames(enrich_res_pos)%in%'evidence_codes']
enrich_res_pos$log_p = -log(as.numeric(enrich_res_pos$p_value))
enrich_res_pos = enrich_res_pos[order(enrich_res_pos$log_p, decreasing = T),]

enrich_res_pos_vis = enrich_res_pos
enrich_res_pos_vis$term_name = gsub('metabolic process', 'metabolism',enrich_res_pos_vis$term_name)
rownames(enrich_res_pos_vis) <- 1:nrow(enrich_res_pos_vis)

### filtering the pathway analysis results based on the following criteria: 
#1.	Prioritize low p-values (high log_p) → to keep the most statistically significant pathways.
#2.	Reduce redundancy by removing terms with high parent overlap or similar names.
#3.	Favor specificity (avoid very broad terms with large term_size unless highly enriched).
#4.	Ensure interpretability → prefer pathway names that are well-defined, specific, and not overly technical.
keep_terms <- c(
  "Collecting duct acid secretion",
  "vacuolar proton-transporting V-type ATPase complex",
  "antiporter activity",
  "active monoatomic ion transmembrane transporter activity",
  "basolateral plasma membrane",
  "regulation of systemic arterial blood pressure",
  "Insulin receptor recycling",
  "intracellular monoatomic ion homeostasis",
  "sodium ion import across plasma membrane",
  "Transport of inorganic cations/anions and amino acids/oligopeptides",
  "proton-transporting two-sector ATPase complex",
  "ATPase dependent transmembrane transport complex",
  "Transferrin endocytosis and recycling",
  "circulatory system process",
  "active transmembrane transporter activity"
)

filtered_res <- enrich_res_pos_vis %>%filter(term_name %in% keep_terms)
library(ggplot2)

# Make sure `filtered_res` is already created as shown previously
ggplot(filtered_res, aes(x = log_p, y = reorder(term_name, log_p))) +
  geom_point(color = "steelblue", size = 4) +
  xlab(expression(-log(p~value))) +
  ylab("") +
  ggtitle("") +
  theme_classic(base_size = 17) +
  theme(
    axis.text.x = element_text(color = "grey20", size = 15),
    axis.text.y = element_text(color = "grey20", size = 16),
    axis.title.x = element_text(color = "grey20", size = 16)
  )


filtered_res = data.frame(filtered_res)
filtered_res[] <- lapply(filtered_res, function(x) {
  if (is.list(x)) sapply(x, paste, collapse = ";") else x
})
write.csv(filtered_res,"~/FLASHMM-analysis/SourceData/Figure3/Figure3Dmale.csv",row.names = FALSE)


##############################################################
###################### FEMALE Pathways ######################
PT_female_df_sort = PT_male_df[order(PT_male_df$score, decreasing = F),]
num_genes=c(200,400)
enrich_res_female_1 = get_gprofiler_enrich(markers=PT_female_df_sort$gene[1:num_genes[1]], model_animal_name)
enrich_res_female_2 = get_gprofiler_enrich(markers=PT_female_df_sort$gene[1:num_genes[2]], model_animal_name)

enrich_res_pos = data.frame(rbind(enrich_res_female_1$result, 
                                  enrich_res_female_2$result))

enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'p_value')]
enrich_res_pos$log_p = -log(enrich_res_pos$p_value)
enrich_res_pos = enrich_res_pos[order(enrich_res_pos$log_p, decreasing = T),]

ggplot(enrich_res_pos, aes(x = log_p, y = reorder(term_name, log_p))) +
  geom_point(color = "pink", size = 4) +
  xlab(expression(-log(p~value))) +
  ylab("") +
  ggtitle("") +
  theme_classic(base_size = 17) +
  theme(
    axis.text.x = element_text(color = "grey20", size = 15),
    axis.text.y = element_text(color = "grey20", size = 16),
    axis.title.x = element_text(color = "grey20", size = 16)
  )
################################

enrich_res_pos = data.frame(enrich_res_pos)
enrich_res_pos[] <- lapply(enrich_res_pos, function(x) {
  if (is.list(x)) sapply(x, paste, collapse = ";") else x
})

write.csv(enrich_res_pos,"~/FLASHMM-analysis/SourceData/Figure3/Figure3Dfemale.csv",row.names = FALSE)

    
