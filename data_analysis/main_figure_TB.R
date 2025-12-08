library(ggplot2)
library(dplyr)
library(gprofiler2)
library(RColorBrewer)
library(FLASHMM)
library(ggrepel)

# Define a muted palette for other cell types
muted_palette <- c(
  "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
  "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE"
)

get_gprofiler_enrich <- function(markers, model_animal_name){
  gostres <- gost(query = markers,
                  ordered_query = TRUE, exclude_iea =TRUE, 
                  sources=c('GO:BP' ,'REAC'),
                  organism = model_animal_name)
  return(gostres)
}

model_animal_name ='hsapiens'
num_genes = 300

counts = readRDS('~/sciFA/Data/Nathan_NatImm_2021.rds')
coldata = counts@meta.data

coldata = readRDS('~/scLMM/LMM-scRNAseq/Data/TB_immune_df.rds')
cell_count_df = data.frame(table(coldata$cluster_name))
cell_count_df = cell_count_df[order(cell_count_df$Freq, decreasing = T),]

length(table(coldata$batch))
length(table(coldata$donor))
table(coldata$sex)
table(coldata$age)
table(coldata$season)

a_cell_type = as.character(cell_count_df$Var1)[1]#'CNT' "CD4+ CD27+CD161+"
coldata$a_cell_type = ifelse(coldata$cluster_name==a_cell_type, a_cell_type, '')
cell_types = c("CD4+ activated" , "CD8+ activated")
coldata$selected_cell_type = ifelse(coldata$cluster_name %in% cell_types, coldata$cluster_name, '')
table(coldata$a_cell_type)
table(coldata$selected_cell_type)


# Define bold colors for selected cell types
selected_colors <- c("CD4+ activated" = "deepskyblue3",
                     "CD8+ activated" = "coral3")


# Assign muted colors to non-selected cell types
other_cell_types <- setdiff(unique(coldata$cluster_name), names(selected_colors))
other_colors <- setNames(rep(muted_palette, length.out = length(other_cell_types)), other_cell_types)
all_colors <- c(selected_colors, other_colors)
# Create the UMAP plot
ggplot(coldata, aes(x = UMAP_1, y = UMAP_2)) +
  # Plot non-selected cells with lighter colors and slightly smaller size
  geom_point(data = coldata %>% filter(!cluster_name %in% names(selected_colors)),
             aes(color = cluster_name), size = 2, alpha = 0.9) + #size = 1.6, alpha = 0.4
  # Plot selected cells with bold colors and slightly larger size
  geom_point(data = coldata %>% filter(cluster_name %in% names(selected_colors)),
             aes(color = cluster_name), size = 2, alpha = 1) + #size = 1, alpha = 0.6
  # Apply custom color scale
  scale_color_manual(values = all_colors) +
  theme_classic() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(
    text = element_text(size = 16),
    legend.title = element_blank()
  )

#### based on sample-type ##### 
coldata$cluster_name[is.na(coldata$cluster_name)] = 'unknown'
table(coldata$cluster_name)
coldata_counts <- ddply(coldata, .(coldata$TB_status, coldata$cluster_name), nrow)
names(coldata_counts) <- c("TB_status", "CellType", "Freq") #Sex, season, TB_status

ggplot(data=coldata_counts, aes(x=CellType, y=Freq, fill=TB_status)) +
  geom_bar(stat="identity",color='black')+theme_classic()+
  #+scale_fill_brewer(palette = "Blues")+
  ylab('Counts')+xlab('Clusters')+
  scale_fill_manual(values = c("gray70", '#FFD92F'))+
  #scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'),
        legend.title = element_blank())+
  xlab('')



coldata_counts_split <- split( coldata_counts , f = coldata_counts$CellType )
coldata_counts_split_norm <- lapply(coldata_counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,coldata_counts_split_norm )
head(counts_norm)

ggplot(data=counts_norm, aes(x=CellType, y=Freq, fill=TB_status)) +
  geom_bar(stat="identity",color='black',alpha=0.9)+theme_classic()+
  ylab('Proportion (%)')+
  #scale_fill_manual(values = c("#E69F00", "#009E73"))+
  #scale_fill_brewer(palette = 'Set2')+
  #scale_fill_manual(values = c("pink1", "skyblue1"))+
  scale_fill_manual(values = c("gray70", '#FFD92F'))+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=15,angle=-90,color='black'),
        legend.title = element_blank()) +  coord_flip()+
  xlab('')

ggplot(coldata, aes(x=Cell_Types_Broad, fill=sex))+
  theme_classic()+geom_bar(stat="identity")+
  theme(text = element_text(size=16),legend.title = element_blank())



age_df = data.frame(donor=coldata$donor,age=coldata$age)
age_df_unique = unique.data.frame(age_df)
age_df_unique = age_df_unique[order(age_df_unique$age, decreasing = F),]
dim(age_df_unique)
hist(age_df_unique$age)

ggplot(age_df_unique, aes(x=age)) + theme_classic()+
  geom_histogram(color="black", fill="grey89")+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16,color='black'),
        axis.text.y = element_text(size=16,color='black'),
        legend.title = element_blank())+xlab('Age')+ylab('Donor count')

######### importing TB data  
#load("~/FLASH-MM/Results_data/TB_lmer_default_intercept.RData") 
load("~/FLASH-MM/Results_data/TB_lmm_intercept.RData") 
#fit <- readRDS('~/scLMM/LMM-scRNAseq-jan2024/lmmfitSS_Nathan_NatImm_2021_X_sexAgeSeason_Z_donor.rds') # Time difference of 1.200262 hours
#fit<- readRDS('~/scLMM/LMM-scRNAseq-jan2024/lmmfitSS_Nathan_NatImm_2021.rds') ## Time difference of 6.184403 hours

tvlmm <- t(fit$t)
pvlmm <- t(fit$p)
sum(apply(is.na(pvlmm), 1, any))
felmm <- t(fit$coef)
slmm <- t(fit$theta)

pvlmm_adj = sapply(1:ncol(pvlmm), function(i) p.adjust(pvlmm[,i], method = "fdr"))
colnames(pvlmm_adj) = colnames(pvlmm)
sum(apply(is.na(pvlmm), 1, any))

P_VAL_THR = 0.05
COEF_THR = 0 # 0.2
cov_marker_list = list()
for(i in 31:ncol(pvlmm_adj)){
  df = data.frame(genes=rownames(pvlmm_adj),
                  tvalue=tvlmm[,i],
                  pvalue=pvlmm[,i],
                  coef = felmm[,i],
                  pvalue_adj=pvlmm_adj[,i],
                  score = -log(pvlmm_adj[,i]+1e-20)*felmm[,i])
  df_ord = df[order(df$score, decreasing = TRUE),]
  #df_ord$is_sig = df_ord$pvalue_adj<P_VAL_THR & (df_ord$coef > COEF_THR | df_ord$coef < -COEF_THR)
  #df_ord$is_sig = df_ord$pvalue_adj<P_VAL_THR & (df_ord$coef > COEF_THR | df_ord$coef < -COEF_THR)
  df_ord$is_sig = df_ord$pvalue_adj<P_VAL_THR & df_ord$coef > COEF_THR
  cov_marker_list[[colnames(pvlmm)[i]]] = df_ord
  
}
lapply(cov_marker_list, head)
res= lapply(cov_marker_list, function(x) sum(x$is_sig))
res= data.frame(numDEG=t(data.frame(res)),stringsAsFactors = FALSE)
res$cellType.cov= rownames(res)
res=res[order(res$numDEG, decreasing = T),]
num_de_df = res
num_de_df$cell_type = gsub('.trt', '',num_de_df$cellType.cov)
num_de_df <- num_de_df %>%
  mutate(cell_type = reorder(cell_type, numDEG))

cell_type_labels <- c(
  "CD4p.activated" = "Activated CD4+ T cells",
  "CD8p.activated" = "Activated CD8+ T cells",
  "CD4p.Treg" = "CD4+ Regulatory T cells",
  "Vd2" = "Vdelta2 T cells",
  "CD4p.HLA_DRp" = "CD4+ HLA-DR+ T cells",
  "Vd1" = "Vdelta1 T cells",
  "CD4s8p.PD_1pTIGITp" = "CD4-CD8+ PD-1+ TIGIT+ T cells",
  "CD8p.GZMKp" = "CD8+ GZMK+ T cells",
  "CD4p.CD161p.Th2" = "CD4+ CD161+ Th2 cells",
  "CD4p.Th2" = "CD4+ Th2 cells",
  "CD4p.RORCp.Treg" = "CD4+ RORC+ Tregs",
  "CD4p.cytotoxic" = "CD4+ Cytotoxic T cells",
  "CD4p.CD38pICOSp.central" = "CD4+ Central memory (CD38+ ICOS+)",
  "CD4p.CD161p.cytotoxic" = "CD4+ CD161+ Cytotoxic T cells",
  "CD8p.GZMBp" = "CD8+ GZMB+ T cells",
  "CD4p.Th17s1" = "CD4+ Th17-like subset 1",
  "CD4p.CD161p.Th1" = "CD4+ CD161+ Th1 cells",
  "CD8p.central" = "CD8+ Central memory T cells",
  "CD4p.Th17" = "CD4+ Th17 cells",
  "CD4p.CCR5p.cytotoxic" = "CD4+ CCR5+ Cytotoxic T cells",
  "CD4p.CCR4p" = "CD4+ CCR4+ T cells",
  "CD8p.CXCR3p" = "CD8+ CXCR3+ T cells",
  "CD4p.CCR4pICOSp.central" = "CD4+ Central memory (CCR4+ ICOS+)",
  "CD4p.CCR4p.central" = "CD4+ Central memory (CCR4+)",
  "CD4p.CD27pCD161p" = "CD4+ CD27+ CD161+ T cells",
  "CD4p.CD27p" = "CD4+ CD27+ T cells",
  "CD4p.Th1" = "CD4+ Th1 cells",
  "CD4p.lncRNA" = "CD4+ lncRNA-high T cells",
  "CD4p.central" = "CD4+ Central memory T cells"
)
# Make sure 'cell_type' is a factor with the correct order
num_de_df$cell_type <- factor(num_de_df$cell_type, 
                              levels = num_de_df$cell_type[order(num_de_df$numDEG, decreasing = FALSE)])

write.csv(num_de_df, '~/FLASHMM-analysis/SourceData/Figure4A.csv', row.names = T, col.names = T)
# Plot with simple character labels
ggplot(num_de_df, aes(x = cell_type, y = numDEG)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Number of TB-specific DE genes per cell type",
       x = "Cell type",
       y = "Number of DE genes") +
  geom_text(aes(label = numDEG), hjust = -0.2, color = "black", size = 3.5) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13)
  ) +
  scale_x_discrete(labels = cell_type_labels)



col_name ='CD4p.activated:trt' #'CD8p.activated:trt' 
cell_trt_df = data.frame(pvalue=pvlmm[,col_name],
                        pvalue_adj = pvlmm_adj[,col_name],
                        tvalue = tvlmm[,col_name],
                        coef = felmm[,col_name],
                        gene=rownames(pvlmm))
cell_trt_df$pvalue_adj_log = -log(cell_trt_df$pvalue_adj+1e-800)
hist(cell_trt_df$tvalue+1e-10)
hist(cell_trt_df$coef)
hist(cell_trt_df$pvalue)
hist(cell_trt_df$pvalue_adj)

cell_trt_df[is.na(cell_trt_df$pvalue_adj),]
sum(cell_trt_df$pvalue_adj<P_VAL_THR & (cell_trt_df$coef > COEF_THR | cell_trt_df$coef < -COEF_THR))
cell_trt_df[cell_trt_df$pvalue_adj<P_VAL_THR & (cell_trt_df$coef > COEF_THR | cell_trt_df$coef < -COEF_THR),]
dim(cell_trt_df)
cell_trt_df = cell_trt_df[cell_trt_df$coef>0,]
dim(cell_trt_df)

cell_trt_df$score = -log(cell_trt_df$pvalue_adj+1e-600)*cell_trt_df$coef
cell_trt_df_ord = cell_trt_df[order(cell_trt_df$score, decreasing = TRUE),]
head(cell_trt_df_ord, 20)

############### Volcano Plot for the CD8p.activated population

# Define inputs
col_name = 'CD8p.activated:trt'
cell_type_name = 'CD8+ Activated T'
P_VAL_THR = 0.05
COEF_THR = 0.1

# Build dataframe
cell_trt_df <- data.frame(
  pvalue = pvlmm[, col_name],
  pvalue_adj = pvlmm_adj[, col_name],
  tvalue = tvlmm[, col_name],
  coef = felmm[, col_name],
  gene = rownames(pvlmm)
)

# Log-adjusted p-value
cell_trt_df$pvalue_adj_log <- -log10(cell_trt_df$pvalue_adj + 1e-800)
# Log p-value
cell_trt_df$pvalue_log <- -log10(cell_trt_df$pvalue + 1e-800)
# Mark upregulated genes
cell_trt_df <- cell_trt_df %>%
  mutate(significance_group = ifelse(pvalue_adj < P_VAL_THR & coef > COEF_THR, "TB-upregulated", "Other"))
# Filter to include only genes to plot (keep all for background, or filter here if you want)
cell_trt_df_plot <- cell_trt_df  # optionally: filter(coef > 0)

# Top 10 TB-upregulated genes
top_up <- cell_trt_df %>%
  filter(significance_group == "TB-upregulated") %>%
  arrange(desc(coef)) %>%
  head(10)


write.csv(cell_trt_df_plot, '~/FLASHMM-analysis/SourceData/Figure4C.csv', row.names = T)

# Volcano plot
ggplot(cell_trt_df_plot, aes(x = coef, y = pvalue_log)) +
  geom_point(data = filter(cell_trt_df_plot, significance_group == "Other"),
             color = "grey80", size = 2.5, alpha = 0.6) +
  geom_point(data = filter(cell_trt_df_plot, significance_group == "TB-upregulated"),
             aes(color = significance_group), size = 2.5, alpha = 0.9) +
  scale_color_manual(values = c("TB-upregulated" = "#D55E00")) +
  labs(
    title = paste("Differential expression in", cell_type_name),
    x = "LogFC",
    y = expression(-log[10]("P-value")),
    color = ""
  ) +
  theme_minimal(base_size = 18) +
  geom_vline(xintercept = COEF_THR, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(P_VAL_THR), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_label_repel(
    data = top_up,
    aes(label = gene),
    box.padding = 0.5,
    point.padding = 0.4,
    label.size = 0.2,
    size = 5,
    segment.color = "black",
    segment.size = 0.5,
    min.segment.length = 0
  )


# Define inputs
col_name = 'CD4p.activated:trt'
cell_type_name = 'CD4+ Activated T'
P_VAL_THR = 0.05
COEF_THR = 0.1

# Build dataframe
cell_trt_df <- data.frame(
  pvalue = pvlmm[, col_name],
  pvalue_adj = pvlmm_adj[, col_name],
  tvalue = tvlmm[, col_name],
  coef = felmm[, col_name],
  gene = rownames(pvlmm)
)

# Log-adjusted p-value
cell_trt_df$pvalue_adj_log <- -log10(cell_trt_df$pvalue_adj + 1e-800)
# Log-adjusted p-value
cell_trt_df$pvalue_log <- -log10(cell_trt_df$pvalue + 1e-800)

# Mark TB-upregulated genes
cell_trt_df <- cell_trt_df %>%
  mutate(significance_group = ifelse(pvalue_adj < P_VAL_THR & coef > COEF_THR, "TB-upregulated", "Other"))

# Top 10 TB-upregulated genes
top_up <- cell_trt_df %>%
  filter(significance_group == "TB-upregulated") %>%
  arrange(desc(coef)) %>%
  head(10)

write.csv(cell_trt_df, '~/FLASHMM-analysis/SourceData/Figure4D.csv', row.names = T)


# Volcano plot
ggplot(cell_trt_df, aes(x = coef, y = pvalue_log)) +
  geom_point(data = filter(cell_trt_df, significance_group == "Other"),
             color = "grey80", size = 2.5, alpha = 0.6) +
  geom_point(data = filter(cell_trt_df, significance_group == "TB-upregulated"),
             aes(color = significance_group), size = 2.5, alpha = 0.9) +
  scale_color_manual(values = c("TB-upregulated" = "royalblue4")) +
  labs(
    title = paste("Differential expression in", cell_type_name),
    x = "LogFC",
    y = expression(-log[10]("P-value")),
    color = ""
  ) +
  theme_minimal(base_size = 18) +
  geom_vline(xintercept = COEF_THR, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(P_VAL_THR), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_label_repel(
    data = top_up,
    aes(label = gene),
    box.padding = 0.5,
    point.padding = 0.4,
    label.size = 0.2,
    size = 5,
    segment.color = "black",
    segment.size = 0.5,
    min.segment.length = 0
  )

################################



#####################################################
col_name ='CD8p.activated:trt' #'CD8p.activated:trt' 
cell_trt_df = data.frame(pvalue=pvlmm[,col_name],
                         pvalue_adj = pvlmm_adj[,col_name],
                         tvalue = tvlmm[,col_name],
                         coef = felmm[,col_name],
                         gene=rownames(pvlmm))
cell_trt_df$pvalue_adj_log = -log(cell_trt_df$pvalue_adj+1e-800)
cell_trt_df$score = -log(cell_trt_df$pvalue_adj+1e-600)*cell_trt_df$coef
cell_trt_df_ord = cell_trt_df[order(cell_trt_df$score, decreasing = TRUE),]

enrich_res_cd8 = get_gprofiler_enrich(markers=cell_trt_df_ord$gene[1:num_genes], model_animal_name)
enrich_res = enrich_res_cd8

head(enrich_res$result,30)
enrich_res_pos = data.frame(enrich_res$result)

enrich_res_pos = enrich_res_pos[,!colnames(enrich_res_pos)%in%'evidence_codes']
enrich_res_pos$log_p = -log(as.numeric(enrich_res_pos$p_value))
enrich_res_pos = enrich_res_pos[order(enrich_res_pos$log_p, decreasing = T),]
enrich_res_pos

# Step 1: Compute log_p and apply basic filters
filtered_enrich <- enrich_res_pos %>%
  mutate(
    log_p = -log10(as.numeric(p_value)),
    enrichment_ratio = intersection_size / term_size
  ) %>%
  filter(
    p_value < 0.05,                                   # statistically significant
    term_size < 1000 | enrichment_ratio > 0.2         # specificity: avoid broad unless very enriched
  ) %>%
  arrange(desc(log_p), term_size)                     # prioritize significance, favor smaller terms
# Step 2: Select top ~10 interpretable pathways
top_enrich <- filtered_enrich %>%
  filter(!duplicated(substr(term_name, 1, 30))) %>%   # crude redundancy filter
  slice_head(n = 10)
# Step 3: Rename terms manually for clarity
rename_map <- c(
  "T cell activation" = "T cell activation",
  "regulation of immune effector process" = "Reg. of immune effector response",
  "granzyme-mediated programmed cell death signaling pathway" = "Granzyme-mediated cell death",
  "T cell receptor signaling pathway" = "TCR signaling",
  "lymphocyte mediated immunity" = "Lymphocyte-mediated immunity",
  "positive regulation of immune effector process" = "Positive reg. Immune effector response",
  "Chemokine receptors bind chemokines" = "Chemokine receptor–ligand binding",
  "regulation of T cell activation" = "Reg. of T cell activation",
  "leukocyte mediated immunity" = "Leukocyte-mediated immunity",
  "positive regulation of immune response" = "Positive reg. Immune response"
)
top_enrich$term_label <- ifelse(
  top_enrich$term_name %in% names(rename_map),
  rename_map[top_enrich$term_name],
  top_enrich$term_name
)
ggplot(top_enrich, aes(x = log_p, y = reorder(term_label, log_p))) +
  geom_point(color = "lightsalmon2", size = 4) +
  theme_minimal(base_size = 15) +
  xlab(expression(-log[10](p-value))) +
  ylab("") +
  ggtitle("") +
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = .5, face = "plain"))


top_enrich = data.frame(top_enrich)
top_enrich[] <- lapply(top_enrich, function(x) {
  if (is.list(x)) sapply(x, paste, collapse = ";") else x
})

write.csv(top_enrich,"~/FLASHMM-analysis/SourceData/Figure4E.csv",row.names = FALSE)


#####################################################
col_name ='CD4p.activated:trt' #'CD8p.activated:trt' 
cell_trt_df = data.frame(pvalue=pvlmm[,col_name],
                         pvalue_adj = pvlmm_adj[,col_name],
                         tvalue = tvlmm[,col_name],
                         coef = felmm[,col_name],
                         gene=rownames(pvlmm))
cell_trt_df$pvalue_adj_log = -log(cell_trt_df$pvalue_adj+1e-800)
cell_trt_df$score = -log(cell_trt_df$pvalue_adj+1e-600)*cell_trt_df$coef
cell_trt_df_ord = cell_trt_df[order(cell_trt_df$score, decreasing = TRUE),]

enrich_res_cd4 = get_gprofiler_enrich(markers=cell_trt_df_ord$gene[1:num_genes], model_animal_name)
enrich_res = enrich_res_cd4

head(enrich_res$result,30)
enrich_res_pos = data.frame(enrich_res$result)

enrich_res_pos = enrich_res_pos[,!colnames(enrich_res_pos)%in%'evidence_codes']
enrich_res_pos$log_p = -log(as.numeric(enrich_res_pos$p_value))
enrich_res_pos = enrich_res_pos[order(enrich_res_pos$log_p, decreasing = T),]
enrich_res_pos

# Step 1: Prepare and filter the enrichment results
enrich_res_cd4_df <- enrich_res_cd4$result %>%
  as.data.frame() %>%
  mutate(
    log_p = -log10(as.numeric(p_value)),
    enrichment_ratio = intersection_size / term_size
  ) %>%
  filter(
    p_value < 0.05,                                   # significance
    term_size < 1000 | enrichment_ratio > 0.2         # specificity
  ) %>%
  arrange(desc(log_p), term_size)

# Step 2: Reduce redundancy (remove similar names)
enrich_res_cd4_df <- enrich_res_cd4_df %>%
  filter(!duplicated(substr(term_name, 1, 30))) %>%   # rough redundancy filter
  slice_head(n = 10)

# Step 3: Rename top pathways manually for interpretability
rename_map_cd4 <- c(
  "mitotic cell cycle" = "Mitotic cell cycle",
  "mitotic cell cycle process" = "Mitosis-related process",
  "cell cycle process" = "Cell cycle process",
  "cell cycle" = "Cell cycle",
  "chromosome organization" = "Chromosome organization",
  "DNA metabolic process" = "DNA metabolism",
  "organelle organization" = "Organelle organization",
  "chromosome segregation" = "Chromosome segregation",
  "DNA replication" = "DNA replication",
  "nuclear division" = "Nuclear division"
)

enrich_res_cd4_df$term_label <- ifelse(
  enrich_res_cd4_df$term_name %in% names(rename_map_cd4),
  rename_map_cd4[enrich_res_cd4_df$term_name],
  enrich_res_cd4_df$term_name
)

# Step 4: Plot
ggplot(enrich_res_cd4_df, aes(x = log_p, y = reorder(term_label, log_p))) +
  geom_point(color = "steelblue3", size = 4) +
  theme_minimal(base_size = 15) +
  xlab(expression(-log[10](p-value))) +
  ylab("") +
  ggtitle("") +
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = .5, face = "plain"))


enrich_res_cd4_df = data.frame(enrich_res_cd4_df)
enrich_res_cd4_df[] <- lapply(enrich_res_cd4_df, function(x) {
  if (is.list(x)) sapply(x, paste, collapse = ";") else x
})

write.csv(enrich_res_cd4_df,"~/FLASHMM-analysis/SourceData/Figure4F.csv",row.names = FALSE)


