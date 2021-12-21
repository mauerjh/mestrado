#####################
### GRUPOS NA W1 ###
#####################

#Definir "Control" como ref para comparacao
dds_uni$condition <- relevel(dds_uni$condition, ref = "A")
#results(dds, c("condition", "A", "C"))

# lists the coefficients
dds_uni <- DESeq(dds_uni)
resultsNames(dds_uni)
BA <- results(dds_uni, name="group_B_vs_A")
CA <- results(dds_uni, name="group_C_vs_A")
DA <- results(dds_uni, name="group_D_vs_A")
sex <- results(dds_uni, name="sex_2_vs_1")
batch <- results(dds_uni, name="batch_2_vs_1")
age <- results(dds_uni, name="age")

age
summary(age)
BA
CA
DA

#############################
### COMPARACAO D x A #####
#############################
group_uni <- as.data.frame(DA)
summary(group_uni)

##VOLCANO PLOT##
#Criar coluna com nome de miRNAs
group_uni$miRNA <- rownames(group_uni)
group_uni$miRNA <- sub("[:MIMAT].*", "",group_uni$miRNA)

#Retirada de NAs
sgroup_uni <- group_uni[complete.cases(group_uni), ]

#Criaï¿½ï¿½o de coluna de up ou down regulated
sgroup_uni$diffexpressed <- "No difference"
sgroup_uni$diffexpressed[sgroup_uni$log2FoldChange > 0.6 & sgroup_uni$padj < 0.05] <- "Up regulated"
sgroup_uni$diffexpressed[sgroup_uni$log2FoldChange < -0.6 & sgroup_uni$padj < 0.05] <- "Down regulated"

sgroup_uni$delabel <- NA
sgroup_uni$delabel[sgroup_uni$diffexpressed != "No difference"] <- sgroup_uni$miRNA[sgroup_uni$diffexpressed != "No difference"]

#Grï¿½fico
pdf("volcanogroup.pdf")
ggplot(sgroup_uni, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point(aes(color = diffexpressed)) + theme(legend.position = "right") + 
  theme_minimal() + 
  geom_vline(xintercept=c(-0.6, 0.6), col="gray") +
  geom_hline(yintercept=-log10(0.05), col="gray") + xlim (-5,5) + 
  scale_color_manual(values= c("cyan4", "black", "coral1")) +
  ggtitle("Differential Expressed miRNAs: Group comparison") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color='Dif. Expressed miRNA') + 
  geom_text_repel()
dev.off()

#####################
### GRUPOS NA W2 ###
#####################


#################################
### COMPARACAO W1 x W2 #####
#################################
#Ajuste da tabela de read counts (OPCIONAL tirar possiveis outliers e seus pares
mRCt_uni <- mRC_uni [,-7]
colnames(mRCt_uni) -> ids
as.data.frame(ids) -> ids
colnames(ids) <- "ID"
sampleinfo <- read.csv2("sample_information_140521.csv", h=T)
join(ids, sampleinfo, by="ID", type="inner") -> sampleinfo_order
rownames(sampleinfo_order) <- sampleinfo_order$ID
#pra nao confundir
rm(sampleinfo)
#checar se a ordem estÃ¡ igual
rownames(sampleinfo_order)
colnames(mRCt_uni)

#tirar colunas que nÃ£o pertencem nessa analise. Deixar: wave, batch, sex, age, subjectid (entender se deixa o grupo tb)
sampleinfo_paired <- sampleinfo_order[, c(-1,-2)]

#DESeq2 quick start -> design da analise
ddst_uni <- DESeqDataSetFromMatrix(countData = mRCt_uni,
                                   colData = sampleinfo_paired,
                                   design = ~ subjectid + wave)

#Selecionar apenas miRNAs com pelo menos 10 amostras com no minimo 10 reads 
keept_uni <- rowSums(counts(ddst_uni) >= 3) >= 18
ddst_uni <- ddst_uni[keept_uni,]

#Extrair a matriz de reads normalizada
ddst_uni <- estimateSizeFactors(ddst_uni)
norm_mRCt_uni <- counts(ddst_uni,normalized = T)

#Renomear as linhas pra ficar mais bonito
rownames(norm_mRCt_uni) <- gsub("[:MIMAT].*", "", rownames(norm_mRCt_uni))

#Definir "Control" como ref para comparacao
ddst_uni$wave <- relevel(ddst_uni$wave, ref = "w1")

# lists the coefficients
ddst_uni <- DESeq(ddst_uni)
resultsNames(ddst_uni)

### Analise diferencial #####
wave_uni <- results(ddst_uni, name="wave_w2_vs_w1")
wave_uni <- as.data.frame(wave_uni)
summary(wave_uni)

