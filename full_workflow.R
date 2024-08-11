library(limma)
library(tidyverse)
setwd("E:/ivg-git/phd/phd-thesis/Dissertation!/")
setwd("D:/git/phd-thesis/Dissertation!/")
# setwd("D:/R/FUR/")
# 
# trgts <- readTargets("targets4.csv", sep = ";")
# 
# rough <-  read.maimages(trgts, source="agilent", 
#                         columns = list(R ="rDyeNormSignal", G = "gDyeNormSignal",rIsFeatNonUnifOL = "rIsFeatNonUnifOL", gIsFeatNonUnifOL="gIsFeatNonUnifOL",rIsBGNonUnifOL= "rIsBGNonUnifOL",gIsBGNonUnifOL="gIsBGNonUnifOL",
#                                        rIsFeatPopnOL="rIsFeatPopnOL",gIsFeatPopnOL="gIsFeatPopnOL",rIsBGPopnOL= "rIsBGPopnOL",
#                                        gIsBGPopnOL="gIsBGPopnOL", rIsSaturated="rIsSaturated",gIsSaturated="gIsSaturated"), 
#                         other.columns =  c("rIsFeatNonUnifOL","gIsFeatNonUnifOL", "rIsBGNonUnifOL","gIsBGNonUnifOL",
#                                            "rIsFeatPopnOL","gIsFeatPopnOL", "rIsBGPopnOL",
#                                            "gIsBGPopnOL", "rIsSaturated","gIsSaturated"), 
#                         annotation = c("accessions","chr_coord","Sequence", 
#                                        "ProbeUID", "ControlType", "ProbeName", "GeneName","SystematicName"
#                                        , "Description"))

rough <- readRDS("rough.rds")
roughbet = normalizeBetweenArrays(rough,method="Gquantile")
roughave <- avereps(roughbet,ID=roughbet$genes$ProbeName)
designRC <- modelMatrix(trgts, ref="Col0")
contrast.matrix <- makeContrasts(contrasts = list("M15", "P12", "M20", "P5", "TM3", "TM11", "Tmp", "TP1" ), levels = designRC)

fitRC <- lmFit(roughave, designRC)
fit2 <- contrasts.fit(fitRC, contrast.matrix)
fit2 <- eBayes(fit2)

DEGs <- list()
coefs <- c("M15", "P12", "M20", "P5", "TM3", "TM11", "Tmp", "TP1")

trueCoefs <- c("M15" = "OEM15", "P12" = "OEP12", "M20" = "OEM20", "P5" = "OEP5", "TM3"="Tmp-M3",
               "TM11" = "Tmp-M11", "Tmp" = "rpotmp", "TP1"= "Tmp-P1")
##########################################################
#for 0.05 fdr

for(i in seq_along(coefs)){
  DEGs[[i]] <-  topTable(fit2, adjust.method = "BH", number = Inf, coef = coefs[i], p.value = 0.05, sort.by = "logFC") %>% filter(ControlType == 0) 
  names(DEGs)[i] <- coefs[i]
}

counts <- list()
for(i in seq_along(DEGs)){
  counts[[i]] <- tibble(line = coefs[i], increased = sum(DEGs[[i]]$logFC > 0), suppressed = sum(DEGs[[i]]$logFC < 0))
}
counts <- do.call(rbind, counts)
counts1 <- counts
counts <- pivot_longer(counts, cols = 2:3, names_to = "regulation", values_to = "DEGs")

cpositions <- tibble(line = counts1$line,increased = (counts1$increased+counts1$suppressed)*3 /4, suppressed = counts1$suppressed/2 )
cpositions <- pivot_longer(cpositions, cols = 2:3, names_to = "regulation", values_to = "position")
cpositions <- left_join(cpositions, counts)
cpositions$color <- ifelse(cpositions$regulation == "increased", "black", "white")

counts$line <- trueCoefs[counts$line]
cpositions$line <- trueCoefs[cpositions$line]
p2 <- ggplot(counts)+
  geom_col(aes(x = line, y = DEGs, fill = regulation), col = "black")+
  geom_text(data = cpositions, aes(x = line,y = position,  label = DEGs), angle = 90, fontface = "bold", size = 5, color = cpositions$color)+
  theme_bw()+scale_fill_manual("Регуляция", values = setNames(c("gold", "firebrick1"), c("increased", "suppressed")))+
  theme(axis.text.x = element_text(size = 17, angle = 90, face = "bold", hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))+ylab("Число генов с дифференциальной экспрессией") +xlab("Линии трансгенных растений")


###################################
#for 1|0.2 fdr

for(i in seq_along(coefs)){
DEGs[[i]] <-  topTable(fit2, lfc = 1, adjust.method = "BH", number = Inf, coef = coefs[i], p.value = 0.2, sort.by = "logFC") %>% filter(ControlType == 0) 
names(DEGs)[i] <- coefs[i]
}

counts <- list()
for(i in seq_along(DEGs)){
  counts[[i]] <- tibble(line = coefs[i], increased = sum(DEGs[[i]]$logFC > 0), suppressed = sum(DEGs[[i]]$logFC < 0))
}
counts <- do.call(rbind, counts)
counts1 <- counts
counts <- pivot_longer(counts, cols = 2:3, names_to = "regulation", values_to = "DEGs")

cpositions <- tibble(line = counts1$line,increased = (counts1$increased+counts1$suppressed)*3 /4, suppressed = counts1$suppressed/2 )
cpositions <- pivot_longer(cpositions, cols = 2:3, names_to = "regulation", values_to = "position")
cpositions <- left_join(cpositions, counts)
cpositions$color <- ifelse(cpositions$regulation == "increased", "black", "white")

counts$line <- trueCoefs[counts$line]
cpositions$line <- trueCoefs[cpositions$line]
p1 <- ggplot(counts)+
  geom_col(aes(x = line, y = DEGs, fill = regulation), col = "black")+
  geom_text(data = cpositions, aes(x = line,y = position,  label = DEGs), angle = 90, fontface = "bold", size = 5, color = cpositions$color)+
  theme_bw()+scale_fill_manual("Регуляция", values = setNames(c("gold", "firebrick1"), c("increased", "suppressed")))+
  theme(axis.text.x = element_text(size = 17, angle = 90, face = "bold", hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))+ylab("Число генов с дифференциальной экспрессией") +xlab("Линии трансгенных растений")

for(i in seq_along(coefs)){
  DEGs[[i]] <-  topTable(fit2, lfc = 1, adjust.method = "BH", number = Inf, coef = coefs[i], p.value = 0.2, sort.by = "logFC") %>% filter(ControlType == 0) 
  names(DEGs)[i] <- coefs[i]
}

ggpubr::ggarrange(p1,p2, nrow = 1, ncol = 2, common.legend = T)

# Venn diagrams

upregulated <- list()
downregulated <- list()
for( i in seq_along(DEGs)){
  upregulated[[i]] <- DEGs[[i]]$ProbeName[DEGs[[i]]$logFC > 0]
  downregulated[[i]] <- DEGs[[i]]$ProbeName[DEGs[[i]]$logFC < 0]
}

names(upregulated) <- trueCoefs
names(downregulated) <- trueCoefs



library(venn)
p1 <- venn(upregulated[c(1,2,5,7,8)], zcolor = "style", box = F, ggplot = T, ilabels = "counts", opacity = 0.2, ilcs = 1, sncs = 1)
p2 <- venn(downregulated[c(1,2,5,7,8)], zcolor = "style", box = F, ggplot = T, ilabels = "counts", opacity = 0.2, ilcs = 1, sncs = 1)

ggpubr::ggarrange(p1,p2)



AAE4 <- function(x){
  library(tidyverse)
  raw_data <- x[x$ControlType==0,]
  raw_data <- cbind(raw_data, AGIDB2[match(raw_data$ProbeName, AGIDB2$ID),c(5:16)])
  raw_data <- cbind(raw_data, agiarel[match(raw_data$ProbeName, agiarel$array_element_name), c(5:9)])
  raw_data$locus_unique <- make.unique(as.character(raw_data$locus), sep = " ")
  library(org.At.tair.db)
  mt2 <- separate_rows(raw_data, locus, sep = ";")
  proxy <- aggregate(SYMBOL~TAIR, AnnotationDbi::select(org.At.tair.db, keys = mt2$locus, keytype = "TAIR", columns = "SYMBOL"), paste, collapse = ";")
  mt2$symbol <- proxy$SYMBOL[match(mt2$locus, proxy$TAIR)]
  mt3 <- aggregate(symbol~ProbeName, mt2, paste, collapse = ";")
  raw_data$symbol <- mt3$symbol[match(raw_data$ProbeName, mt3$ProbeName)]
  raw_data$label <- NA
  raw_data$label[raw_data$adj.P.Val < 0.05] <- "*"
  return(raw_data)}

AAE5 <- function(x){
  library(tidyverse)
  raw_data <- x[x$ControlType==0,]
  raw_data <- cbind(raw_data, AGIDB2[match(raw_data$ProbeName, AGIDB2$ID),c(5:16)])
  raw_data <- cbind(raw_data, agiarel[match(raw_data$ProbeName, agiarel$array_element_name), c(5:9)])
  raw_data$locus_unique <- make.unique(as.character(raw_data$locus), sep = " ")
  library(org.At.tair.db)
  mt2 <- separate_rows(raw_data, locus, sep = ";")
  proxy <- aggregate(SYMBOL~TAIR, AnnotationDbi::select(org.At.tair.db, keys = mt2$locus, keytype = "TAIR", columns = "SYMBOL"), paste, collapse = ";")
  mt2$symbol <- proxy$SYMBOL[match(mt2$locus, proxy$TAIR)]
  mt3 <- aggregate(symbol~ProbeName, mt2, paste, collapse = ";")
  raw_data$symbol <- mt3$symbol[match(raw_data$ProbeName, mt3$ProbeName)]
  raw_data$label <- NA
  raw_data$label[raw_data$adj.P.Val < 0.05] <- "*"
  raw_data <- separate_rows(raw_data, "locus",sep= ";")
  raw_data$prot_func <- genefam$Protein_Function[match(raw_data$locus, genefam$Genomic_Locus_Tag)]
  raw_data$gene_family <- genefam$Gene_Family[match(raw_data$locus, genefam$Genomic_Locus_Tag)]
  raw_data$gene_subfamily <- genefam$Sub_Family[match(raw_data$locus, genefam$Genomic_Locus_Tag)]
  raw_data$gene_name2 <- genefam$Gene_Name[match(raw_data$locus, genefam$Genomic_Locus_Tag)]
  raw_data$gene_fam_desc <- genefam$Gene_Family_Criteria[match(raw_data$locus, genefam$Genomic_Locus_Tag)]
  raw_data$description2 <- fdesc$Short_description[match(raw_data$locus, fdesc$Model_name)]
  raw_data$curator_summary <- fdesc$Curator_summary[match(raw_data$locus, fdesc$Model_name)]
  raw_data$tair_comp_desc <- fdesc$Computational_description[match(raw_data$locus, fdesc$Model_name)]
  return(raw_data)}

AAE6 <- function(x){
  library(tidyverse)
  raw_data <- x[x$ControlType==0,]
  raw_data <- cbind(raw_data, AGIDB2[match(raw_data$ProbeName, AGIDB2$ID),c(5:16)])
  raw_data <- cbind(raw_data, agiarel[match(raw_data$ProbeName, agiarel$array_element_name), c(5:9)])
  library(org.At.tair.db)
  mt2 <- separate_rows(raw_data, locus, sep = ";")
  proxy <- aggregate(SYMBOL~TAIR, AnnotationDbi::select(org.At.tair.db, keys = mt2$locus, keytype = "TAIR", columns = "SYMBOL"), paste, collapse = ";")
  raw_data <- separate_rows(raw_data, "locus",sep= ";")
  raw_data$symbol <- proxy$SYMBOL[match(raw_data$locus, proxy$TAIR)]
  raw_data$label <- NA
  raw_data$label[raw_data$adj.P.Val < 0.05] <- "*"
  raw_data$prot_func <- genefam$Protein_Function[match(raw_data$locus, genefam$Genomic_Locus_Tag)]
  raw_data$gene_family <- genefam$Gene_Family[match(raw_data$locus, genefam$Genomic_Locus_Tag)]
  raw_data$gene_subfamily <- genefam$Sub_Family[match(raw_data$locus, genefam$Genomic_Locus_Tag)]
  raw_data$gene_name2 <- genefam$Gene_Name[match(raw_data$locus, genefam$Genomic_Locus_Tag)]
  raw_data$gene_fam_desc <- genefam$Gene_Family_Criteria[match(raw_data$locus, genefam$Genomic_Locus_Tag)]
  raw_data$description2 <- fdesc$Short_description[match(raw_data$locus, fdesc$Model_name)]
  raw_data$curator_summary <- fdesc$Curator_summary[match(raw_data$locus, fdesc$Model_name)]
  raw_data$tair_comp_desc <- fdesc$Computational_description[match(raw_data$locus, fdesc$Model_name)]
  return(raw_data)}

Expressionne <- topTable(fit2, lfc = 1, adjust.method = "BH", number = Inf, p.value = 1, coef = c("M15", "P12", "M20", "P5", "TM3", "TM11", "Tmp", "TP1" ))
Expressionne <- AAE6(Expressionne)
Expressionne <- Expressionne[-grep("no_match", Expressionne$locus),]

#GO Overrepresentation test:

annotated <- list()
for( i in seq_along(DEGs)){
  annotated[[i]] <- AAE6(DEGs[[i]])
  cat(paste("iter",i),"\r")
}
names(annotated) <- trueCoefs

for(i in seq_along(annotated)){
  annotated[[i]] <- annotated[[i]] %>% filter(locus!= "no_match")
}

upregulated <- list()
downregulated <- list()
for(i in seq_along(annotated)){
  upregulated[[i]] <- annotated[[i]] %>% filter(logFC > 0) %>% pull(locus) %>% unique()
  downregulated[[i]] <- annotated[[i]] %>% filter(logFC < 0) %>% pull(locus) %>% unique()
}
names(upregulated) <- trueCoefs
names(downregulated) <- trueCoefs

library(clusterProfiler)
library(org.At.tair.db)

upres <- list()
dores <- list()
for(i in seq_along(annotated)){
  upres[[i]] <- clusterProfiler::simplify(enrichGO(gene = upregulated[[i]], OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH"))
  dores[[i]] <- clusterProfiler::simplify(enrichGO(gene = downregulated[[i]], OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH"))
cat(paste("iter", i),"\r")
  }

for(i in seq_along(annotated)){
  upres[[i]] <- as_tibble(upres[[i]]) %>% mutate(line = trueCoefs[i], regulation = "upregulated") 
  dores[[i]] <- as_tibble(dores[[i]]) %>% mutate(line = trueCoefs[i], regulation = "downregulated")
  cat(paste("iter", i),"\r")
}

res <- rbind(do.call(rbind, upres), do.call(rbind, dores)) %>% mutate(p.sig = -log10(p.adjust))
library(scales)
ggplot(res)+
  geom_point(aes(x = p.sig, size = Count, y = Description, col = regulation))+facet_wrap(~line, nrow = 1)+
  theme_bw()+  scale_y_discrete(labels = label_wrap(40))+xlab("Значимость Р")+
  scale_color_manual("Регуляция", values = setNames(c("firebrick1", "dodgerblue"), c("downregulated", "upregulated")))+
  theme(strip.background.x = element_rect(fill = "white", colour = "black"),
        strip.text.x = element_text(face = "bold")) +ylab("Термины GO Biological Procecces")

#GO Geneset Enrichment Analysis:

#need to think about it more...

fgs <- list()
for(i in seq_along(coefs)){
 temp <-  tibble(LFC = Expressionne[[coefs[i]]], locus = Expressionne$locus, ProbeName = Expressionne$ProbeName) %>% filter(LFC!=0) %>% arrange(-LFC)
  
  fgs[[i]] <- setNames(temp$LFC, temp$locus)
}
names(fgs) <- trueCoefs

#unfortunately not enough RAM for parallel....

# res2 <- list()
# calculate_gsea <- function(i){
#   temp <- gseGO(geneList = fgs[[i]], ont = "BP", OrgDb = org.At.tair.db, keyType = "TAIR", maxGSSize = 10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
#   temp <- as_tibble(temp) %>% mutate(line = trueCoefs[i])
#   return(temp)
#   }
# 
# library(parabar)
# backend <- start_backend(11)
# evaluate(backend, expression = {library(clusterProfiler); library(org.At.tair.db)})
# export(backend, variables = c("trueCoefs", "fgs", "calculate_gsea"))
# gc()
# res2 <- par_lapply(backend, x = seq_along(fgs), fun = calculate_gsea)
# 
# stop_backend(backend)
# rm(backend)



res2 <- list()
for(i in seq_along(fgs)){
  temp <- clusterProfiler::simplify(gseGO(geneList = fgs[[i]], ont = "BP", OrgDb = org.At.tair.db, keyType = "TAIR", maxGSSize = 10000, minGSSize = 1, pvalueCutoff = 0.05, pAdjustMethod = "BH"))
    temp <- as_tibble(temp) %>% mutate(line = trueCoefs[i])
    res2[[i]] <- temp
    }

res3 <- do.call(rbind, res2) %>% mutate(p.sig = -log10(p.adjust))


library(scales)
ggplot(res3)+
  geom_point(aes(x = enrichmentScore, size = rank, y = Description, col = p.sig))+facet_wrap(~line, nrow = 1)+
  scale_color_gradientn("?????????? ?", colours = c("royalblue", "firebrick1"))+
  theme_bw()+  scale_y_discrete(labels = label_wrap(40))+xlab("")+
  theme(strip.background.x = element_rect(fill = "white", colour = "black"),
        strip.text.x = element_text(face = "bold")) +ylab("??????? GO Biological Procecces")+scale_size_continuous("????")

openxlsx::write.xlsx(res3, "full_gsea_go.xlsx")

res4 <- res3 %>%  group_by(line) %>% slice_max(order_by = abs(enrichmentScore), n = 5) %>% ungroup()  
res4$cenr <- NA_integer_
for(i in 1:nrow(res4)){
  res4$cenr[i] <- length(unlist(str_split(res4$core_enrichment[i], pattern = "/")))
}

ggplot(res4)+
  geom_point(aes(x = enrichmentScore, size = cenr, y = Description, col = p.sig))+
  geom_vline(aes(xintercept = 0), col = "firebrick", linetype = "dashed")+
facet_wrap(~line, nrow = 1)+
  scale_color_gradientn("?????????? ?", colours = c("royalblue", "firebrick1"))+
  theme_bw()+  scale_y_discrete(labels = label_wrap(40))+xlab("?????? ??????????")+
  theme(strip.background.x = element_rect(fill = "white", colour = "black"),
        strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +ylab("??????? GO Biological Procecces")+
    scale_size_continuous("???? ??????????")

colnames(Expressionne)[10:17]
trueCoefs
#KEGG enrichment

kres <- list()
kres2 <- list()

for(i in seq_along(trueCoefs)){
temp1 <- enrichKEGG(gene = upregulated[[i]], organism = "ath") %>% as_tibble() %>% mutate(line = trueCoefs[i], regulation = "upregulated")   
temp2 <- enrichKEGG(gene = downregulated[[i]], organism = "ath") %>% as_tibble() %>% mutate(line = trueCoefs[i], regulation = "downregulated")   
kres[[i]] <- temp1
kres2[[i]] <- temp2
cat(paste("iter",i),"\r")}

kres3 <- rbind(do.call(rbind, kres),do.call(rbind, kres2))%>% mutate(p.sig = -log10(p.adjust))


library(scales)
ggplot(kres3)+
  geom_point(aes(x = p.sig, size = Count, y = Description, col = regulation))+facet_wrap(~line, nrow = 1)+
  theme_bw()+  scale_y_discrete(labels = label_wrap(40))+xlab("?????????? P")+
  scale_color_manual("?????????", values = setNames(c("firebrick1", "dodgerblue"), c("downregulated", "upregulated")))+
  theme(strip.background.x = element_rect(fill = "white", colour = "black"),
        strip.text.x = element_text(face = "bold")) +ylab("?????????????? ???? KEGG")

BiocManager::install("pathview")
library(pathview)
data(gse16873.d)
try <-  as.data.frame(fgs[["rpotmp"]]) %>% mutate(locus = names(fgs[["rpotmp"]]))
colnames(try)[1] <- "rpotmp"
try <- try %>% group_by(locus) %>% summarise(rpotmp = mean(rpotmp)) %>% dplyr::arrange(-rpotmp)
try <- as.data.frame(try)
rownames(try) <- try$locus
try2 <- matrix(try[,-1], nrow = nrow(try), dimnames = list(try$locus, "rpotmp"))

pv.out <- pathview(gene.data =try2, pathway.id = "04141",gene.idtype = "TAIR",
                   species = "ath", out.suffix = "rpotmp", low = "firebrick1", mid = "white", high = "dodgerblue")

pv.out <- pathview(gene.data =try2, pathway.id = "00480",gene.idtype = "TAIR",
                   species = "ath", out.suffix = "rpotmp", low = "firebrick1", mid = "white", high = "dodgerblue")

#preparing the heatmap-dataset

  Expressionne$M15[Expressionne$locus %in% setdiff(Expressionne$locus, c(upregulated$OEM15 , downregulated$OEM15))] <- 0
Expressionne$P12[Expressionne$locus %in% setdiff(Expressionne$locus, c(upregulated$OEP12 , downregulated$OEP12))] <- 0
Expressionne$M20[Expressionne$locus %in% setdiff(Expressionne$locus, c(upregulated$OEM20 , downregulated$OEM20))] <- 0
Expressionne$P5[Expressionne$locus %in% setdiff(Expressionne$locus, c(upregulated$OEP5 , downregulated$OEP5))] <- 0
Expressionne$TM3[Expressionne$locus %in% setdiff(Expressionne$locus, c(upregulated$`Tmp-M3` , downregulated$`Tmp-M3`))] <- 0
Expressionne$TM11[Expressionne$locus %in% setdiff(Expressionne$locus, c(upregulated$`Tmp-M11` , downregulated$`Tmp-M11`))] <- 0
Expressionne$Tmp[Expressionne$locus %in% setdiff(Expressionne$locus, c(upregulated$rpotmp , downregulated$rpotmp))] <- 0
Expressionne$TP1[Expressionne$locus %in% setdiff(Expressionne$locus, c(upregulated$`Tmp-P1` , downregulated$`Tmp-P1`))] <- 0
Expressionne <- Expressionne[Expressionne$M15 !=0 |
                               Expressionne$P12 !=0 |
                               Expressionne$M20 !=0 |
                               Expressionne$P5 !=0 |
                               Expressionne$TM3 !=0 |
                               Expressionne$TM11 !=0 |
                               Expressionne$Tmp !=0 |
                               Expressionne$TP1 !=0 ,]

# openxlsx::write.xlsx(Expressionne, "Transcriptome.merged.xlsx")

Expressionne <- openxlsx::read.xlsx("Transcriptome.merged.xlsx")

exp_true <- Expressionne %>% group_by(locus) %>% summarise(OEM15 = mean(M15),
                                                           OEP12 = mean(P12),
                                                           OEM20 = mean(M20),
                                                           OEP5 = mean(P5),
                                                           rpotmp = mean(Tmp),
                                                           `Tmp-M3`= mean(TM3),
                                                           `Tmp-M11`= mean(TM11),
                                                           `Tmp-P1`=mean(TP1)) %>% mutate(symbol = Expressionne$symbol[match(locus, Expressionne$locus)])


                                                           
openxlsx::write.xlsx(exp_true, "exp_true.xlsx")                                                           

#transcription factors

tf <- Expressionne %>% filter(locus %in% TF_list$Gene_ID) %>% mutate(tf_family = TF_list$Family[match(locus, TF_list$Gene_ID)])
colnames(tf)[10:17] <- trueCoefs
range(tf[,c(10:17)])
library(ComplexHeatmap)

tf2 <- tf %>% group_by(locus) %>% summarise(OEM15 = mean(OEM15),
                                      OEP12 = mean(OEP12),
                                      OEM20 = mean(OEM20),
                                      OEP5 = mean(OEP5),
                                      rpotmp = mean(rpotmp),
                                      `Tmp-M3`= mean(`Tmp-M3`),
                                      `Tmp-M11`= mean(`Tmp-M11`),
                                      `Tmp-P1`=mean(`Tmp-P1`)) %>% mutate(symbol = tf$symbol[match(locus, tf$locus)],
                                                                          tf_family = tf$tf_family[match(locus, tf$locus)])

length(unique(tf2$tf_family))
table(tf2$tf_family)

tf2 %>% filter(OEP12 > 0) %>% nrow()
tf2 %>% filter(OEP12 < 0) %>% nrow()


length(intersect(upregulated$OEM15, upregulated$OEP12))
length(intersect(downregulated$OEM15, downregulated$OEP12))

col_fun5 = circlize::colorRamp2(c(4,2,0,-2,-4), c( "gold", "chartreuse", "black", "blue", "firebrick1"))
Heatmap(as.matrix(tf[,c(10:17)]), col = col_fun5, row_split = tf$tf_family,
        row_labels = tf$Description, row_gap = unit(0.5,"mm"),
        row_names_gp = gpar(fontsize = 0), row_title_rot = 0, row_title_gp = gpar(fontsize = 7),  
        right_annotation = HeatmapAnnotation(`Family` = tf$tf_family,
                                             which = "row"
        )
)

openxlsx::write.xlsx(tf, "tf_expression_full.xlsx")
tf <- openxlsx::read.xlsx("tf_expression_full.xlsx")

tf2 <- exp_true %>% filter(locus %in% tf$locus)


sum(tf2$OEM15 > 0)
sum(tf2$OEM15 < 0)

tf3 <- tf2[abs(tf2$OEM15) >1.5 |
            abs(tf2$OEP12) >1.5 |
            abs(tf2$OEM20) >1.5 |
            abs(tf2$OEP5) >1.5 |
            abs(tf2$rpotmp) >1.5 |
            abs(tf2$`Tmp-M3`) >1.5 |
            abs(tf2$`Tmp-M11`) >1.5 |
            abs(tf2$`Tmp-P1`) >1.5 ,]



Heatmap(as.matrix(tf3[,c(2:9)]), col = col_fun5, row_split = tf3$tf_family, name = "LFC",
        row_labels = str_split_fixed(tf3$symbol, ";",2)[,1], row_gap = unit(0.5,"mm"),
        row_names_gp = gpar(fontsize = 5), row_title_rot = 0, row_title_gp = gpar(fontsize = 7),  
        right_annotation = HeatmapAnnotation(`Family` = tf3$tf_family,
                                             which = "row"
        )
)
uptf <- list()
dtf <- list()
for(i in seq_along(trueCoefs)){
temp <- upregulated[[i]]
temp <- temp[temp %in% TF_list$Gene_ID]
uptf[[i]] <- temp

temp <- downregulated[[i]]
temp <- temp[temp %in% TF_list$Gene_ID]
dtf[[i]] <- temp
}
names(uptf) <- trueCoefs
names(dtf) <- trueCoefs
library(venn)
p1 <- venn(uptf[c(1,2,5,7,8)], zcolor = "style", ggplot = T, box = F , ilabels = "counts", ilcs = 1.5, sncs = 1.5)
p2 <- venn(dtf[c(1,2,5,7,8)], zcolor = "style", ggplot = T, box = F , ilabels = "counts", ilcs = 1.5, sncs = 1.5)
ggpubr::ggarrange(p1,p2, ncol = 2, labels = c("UP", "DOWN"))

length(intersect(uptf$OEM15, uptf$OEP12))
length(intersect(dtf$OEM15, dtf$OEP12))

library(igraph)
library(ggraph)
vdf <- exp_true %>% rename("name" = "locus")
tfdb2 <- tfdb[tfdb$V1 %in% exp_true$locus & tfdb$V3 %in% exp_true$locus , ]
tfg <- graph_from_data_frame(tfdb2[,c(1,3)], vertices = exp_true)

m15 <- delete_vertices(tfg, which(abs(V(tfg)$OEM15) <1))
m15 <- delete_vertices(m15, which(degree(m15)<1))
V(m15)$type <- ifelse(V(m15)$name %in% tfdb2$V1, "TF", "Target")
V(m15)$reg <- ifelse(V(m15)$OEM15 > 0, "upregulated", "downregulated")
eddf <- as_data_frame(m15, "edges")

E(m15)$type = ifelse(V(m15)$OEM15[match(eddf$from, V(m15)$name)] * V(m15)$OEM15[match(eddf$to, V(m15)$name)] > 0, "activation", "repression")

# ggraph(m15, layout = "fr")+geom_edge_fan(alpha = 0.3, aes(color = type), edge_width = 0.1)+
#   geom_node_point(aes(col = OEM15, shape = type, size = degree(m15)))+theme_graph()+scale_size_continuous("Degree")+
#   scale_color_gradientn(colors= c("firebrick1", "dodgerblue", "black", "springgreen", "gold"), values = c(0, 0.1, 0.37, 0.5, 1))+
#   scale_edge_color_manual("Regulation", values = setNames(c("firebrick1", "royalblue3"), c("repression", "activation")))

ggraph(m15, layout = "fr")+geom_edge_fan(alpha = 0.3, aes(color = type), edge_width = 0.4)+
  geom_node_point(aes(col = reg, shape = type, size = degree(m15)))+theme_graph()+scale_size_continuous("Degree")+
  scale_color_manual(values= setNames(c("dodgerblue", "firebrick1"), c("upregulated", "downregulated")))+
  scale_edge_color_manual("Regulation", values = setNames(c("firebrick1", "royalblue3"), c("repression", "activation")))

test <- as_data_frame(m15, "vertices") %>% mutate(localization = suba2$location_consensus[match(name, suba2$locus)])

test %>% filter(type == "Target", reg== "upregulated") %>% nrow()
test %>% filter(type == "Target", reg== "downregulated") %>% nrow()

testu <- test %>% filter(type == "Target", reg== "upregulated")
testd <- test %>% filter(type == "Target", reg== "downregulated")

length(grep("mitoc", testu$localization))
length(grep("mitoc", testd$localization))


tfreg <- list()
ver <- as_data_frame(tfg, "vertices")
for(i in seq_along(trueCoefs)){
  
  m15 <- delete_vertices(tfg, which(abs(ver[[trueCoefs[i]]]) <1))
  m15 <- delete_vertices(m15, which(degree(m15)<1))
  V(m15)$type <- ifelse(V(m15)$name %in% tfdb2$V1, "TF", "Target")
  ver2 <- as_data_frame(m15, "vertices")
  V(m15)$reg <- ifelse(ver2[[trueCoefs[i]]] > 0, "upregulated", "downregulated")
  eddf <- as_data_frame(m15, "edges")
  E(m15)$type = ifelse(ver2[[trueCoefs[i]]][match(eddf$from, ver2$name)] * ver2[[trueCoefs[i]]][match(eddf$to, ver2$name)] > 0, "activation", "repression")
  tfreg[[i]] <- m15
}
names(tfreg) <- trueCoefs
# saveRDS(tfreg, "tfreg.rds")

plots <- list()
for(i in seq_along(tfreg)){
  
plots[[i]] <- ggraph(tfreg[[i]], layout = "fr")+geom_edge_fan(alpha = 0.3, aes(color = type), edge_width = 0.1)+
  geom_node_point(aes(col = reg, shape = type), size = 0.5, alpha = 0.8)+theme_graph()+scale_size_continuous("Degree")+
  scale_color_manual(values= setNames(c("chartreuse", "tomato"), c("upregulated", "downregulated")))+
  scale_edge_color_manual("Regulation", values = setNames(c("firebrick1", "royalblue3"), c("repression", "activation")))

}

ggpubr::ggarrange(plotlist = plots, ncol = 4, nrow = 2, common.legend = T, labels = trueCoefs)

which.max(degree(tfreg[[1]]))

sbg <- make_ego_graph(tfreg[[1]], order = 1, nodes = "AT4G38000")[[1]]

table(E(sbg)$type)

as_data_frame(tfreg[[1]], "vertices") %>% filter(name == "AT4G38000")
degree(tfreg[[1]])[which.max(degree(tfreg[[1]]))]

#organellar proteins

mito <- rbind(exp_true[grepl("ATMG", exp_true$locus),], exp_true[grepl("AT2G07", exp_true$locus),])
mito$type <- NA_character_
mito$type[grepl("NAD", mito$symbol)] <- "Complex I"
mito$type[grepl("COB", mito$symbol)] <- "Complex III"
mito$type[grepl("ATP", mito$symbol)] <- "Complex V"
mito$type[grepl("COX", mito$symbol)] <- "Complex IV"
mito$type[grepl("RP", mito$symbol)] <- "Ribosome"
mito$type[is.na(mito$type)==T] <- "Other"


range(mito[,2:9])


tmine <- read.csv("summary.csv", sep = " ")

mito$type[grepl("ORF240A", mito$symbol)] <- "Complex V"
mito$type[grepl("AT2G07711", mito$locus)] <- "Complex I"
mito$symbol[grepl("AT2G07711", mito$locus)] <- "NAD5B"

mito$type[grepl("AT2G07785", mito$locus)] <- "Complex I"
mito$symbol[grepl("AT2G07785", mito$locus)] <- "NAD1_ps"

mito <- mito[-grep("AT2G07785", mito$locus),]

mito$type[grepl("AT2G07786", mito$locus)] <- "Complex I"
mito$symbol[grepl("AT2G07786", mito$locus)] <- "NAD1A"

mito$type[grepl("AT2G07709", mito$locus)] <- "Complex I"
mito$symbol[grepl("AT2G07709", mito$locus)] <- "NAD7"

mito$type[grepl("AT2G07731", mito$locus)] <- "Complex I"
mito$symbol[grepl("AT2G07731", mito$locus)] <- "NAD6"

mito <- mito[-grep("AT2G07731", mito$locus),]

mito$type[grepl("AT2G07733", mito$locus)] <- "Complex I"
mito$symbol[grepl("AT2G07733", mito$locus)] <- "NAD2B"

mito$type[grepl("AT2G07751", mito$locus)] <- "Complex I"
mito$symbol[grepl("AT2G07751", mito$locus)] <- "NAD3"

mito <- mito[-grep("AT2G07751", mito$locus),]

mito$type[grepl("AT2G07733", mito$locus)] <- "Complex I"
mito$symbol[grepl("AT2G07733", mito$locus)] <- "NAD2B"

mito$type[grepl("AT2G07733", mito$locus)] <- "Complex I"
mito$symbol[grepl("AT2G07733", mito$locus)] <- "NAD2B"

mito <- mito[-grep("AT2G07675", mito$locus),] #rpl2 pseudo

mito$type[grepl("AT2G07726", mito$locus)] <- "Ribosome"
mito$symbol[grepl("AT2G07726", mito$locus)] <- "RPS14"

mito$type[grepl("AT2G07725", mito$locus)] <- "Ribosome"
mito$symbol[grepl("AT2G07725", mito$locus)] <- "RPL5P"

mito$type[grepl("AT2G07530", mito$locus)] <- "Other"
mito$symbol[grepl("AT2G07530", mito$locus)] <- "ORF1"

mito <- mito[-grep("AT2G07725", mito$locus),] #rpl5 pseudo
mito <- mito[-grep("AT2G07734", mito$locus),] #rps4 pseudo
mito <- mito[-grep("AT2G07712", mito$locus),] #matr pseudo
mito <- mito[-grep("AT2G07747", mito$locus),] #matr pseudo
mito <- mito[-grep("AT2G07776", mito$locus),] # pseudo

mito <- mito[-grep("AT2G07773", mito$locus),] #orf215A pseudo
mito2 <- mito[-which(is.na(mito$symbol)==T & grepl("AT2G", mito$locus)),]

mito2$type[grepl("ATMG00020", mito2$locus)] <- "rRNA"
mito2$symbol[grepl("ATMG00020", mito2$locus)] <- "RRN26"

mito2$type[grepl("ATMG00516", mito2$locus)] <- "Complex I"
mito2$symbol[grepl("ATMG00516", mito2$locus)] <- "NAD1C"

# ATMG00900 - ccmC
# ATMG00570 - mttB
# ATMG00980 - rpS12

mito2$type[grepl("ATMG00900", mito2$locus)] <- "Other"
mito2$symbol[grepl("ATMG00900", mito2$locus)] <- "ccmC"

mito2$type[grepl("ATMG00570", mito2$locus)] <- "Other"
mito2$symbol[grepl("ATMG00570", mito2$locus)] <- "mttB"

mito2$type[grepl("ATMG00980", mito2$locus)] <- "Ribosome"
mito2$symbol[grepl("ATMG00980", mito2$locus)] <- "rpS12"

openxlsx::write.xlsx(mito2,"mtgenes.xlsx")
mito2 <- openxlsx::read.xlsx("mtgenes.xlsx")

library(ComplexHeatmap)
col_fun5 = circlize::colorRamp2(c(4,1.5,0,-1.5,-4), c( "gold", "chartreuse", "black", "blue", "firebrick1"))
Heatmap(as.matrix(mito2[,c(2:9)]), col = col_fun5, row_split = mito2$type,
        row_labels = str_split_fixed(mito2$symbol, ";",2)[,1], row_gap = unit(0.5,"mm"),
        row_names_gp = gpar(fontsize = 5, fontface = "italic"), row_title_rot = 0, row_title_gp = gpar(fontsize = 7),  
        right_annotation = HeatmapAnnotation(`Type` = mito2$type,
                                             which = "row", col = list(`Type`=setNames(c("grey", "tomato", "brown", "dodgerblue", "darkmagenta", "chartreuse", "palegreen"),
                                                                                       c("Other", "Ribosome", "rRNA", "Complex I", "Complex III", "Complex IV", "Complex V"))
        )
))

suba <- read.csv("suba5/19_8_2022_suba5_summary.csv") %>% select(c(1,3,40)) %>% mutate(locus = str_split_fixed(locus, "\\.",2)[,1]) %>% as_tibble() %>% distinct()

suba2 <- aggregate(suba, location_consensus~locus, FUN = function(x){paste(x, collapse = ";")})
suba2$description <- suba$description[match(suba2$locus, suba$locus)]

msuba <- suba2[grep("mito", suba2$location_consensus),] 
msuba <- msuba[-grep("ATMG", msuba$locus),]
msuba <- msuba[-grep("AT2G07", msuba$locus),]

mne <- exp_true %>% filter(locus %in% msuba$locus) %>% mutate(location = msuba$location_consensus[match(locus, msuba$locus)],
                                                              description =msuba$description[match(locus, msuba$locus)] ,
                                                              sdesc = fdesc$Short_description[match(locus, fdesc$Model_name)],
                                                              family = genefam$Gene_Family[match(locus, genefam$Genomic_Locus_Tag)],
                                                              funk = genefam$Protein_Function[match(locus, genefam$Genomic_Locus_Tag)])

mne$location <- gsub(",", ";",mne$location)
mne$funk[mne$locus %in% c("AT2G20800", "AT4G21490")] <- "Alternative NADH-dehydrogenase B"
mne$funk[which(mne$funk == "NULL")] <- mne$family[which(mne$funk == "NULL")]
table(mne$location)

openxlsx::write.xlsx(mne, "mt-nuclear-encoded.full.xlsx")
mne <- openxlsx::read.xlsx("mt-nuclear-encoded.full.xlsx")

mne %>% filter(OEM15 < 0) %>% nrow()

library(org.At.tair.db)
go_mito <- AnnotationDbi::select(org.At.tair.db, keys = "GO:0005739", keytype = "GO", columns = "TAIR")

exp_true2 <- exp_true[,c(1:3,6,7,9,10)]
exp_true2 <- exp_true2[rowSums(exp_true2[,c(2:6)])>0,]


mne2 <- Expressionne[Expressionne$locus %in% go_mito$TAIR,]
mne2 <- mne2[-grep("ATMG", mne2$locus),]
mne2 <- mne2[-grep("AT2G07", mne2$locus),]
mne2 <- mne2[rowSums(mne2[,c(10:17)])>0,]

mne <- openxlsx::read.xlsx("Supplementary Data S4.xlsx", startRow = 7)
mne2 <- mne[grep("Mito", mne$localization),]
mne2 <- mne2[-grep("ATMG", mne2$locus),]

mne2 <- mne[grep("Chlo", mne$localization),]

mne2$func <- genefam$Gene_Family[match(mne2$locus, genefam$Genomic_Locus_Tag)]
mne2$gene_symbol[is.na(mne2$gene_symbol)==T] <- "-"

library(ComplexHeatmap)
col_fun5 = circlize::colorRamp2(c(4,1.5,0,-1.5,-4), c( "gold", "chartreuse", "black", "dodgerblue", "firebrick1"))

Heatmap(as.matrix(mne2[,c(2:6)]), col = col_fun5, row_split = mne2$func, name = "LFC",
        row_labels = paste(mne2$locus, str_split_fixed(mne2$gene_symbol, ";",2)[,1]), row_gap = unit(0.1,"mm"),
        row_names_gp = gpar(fontsize = 6), row_title_rot = 0, row_title_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 25),
        right_annotation = HeatmapAnnotation(`Localization` = mne2$localization,
                                             # `Family` = mne2$family,
                                              `Family` = mne2$func,
                                             which = "row", col = list(`Localization`= setNames(c ("springgreen", "tomato"), unique(mne2$localization))
                                                                       , `Family` = setNames( viridis::turbo(26), unique(mne2$func)[-4]))
                                                                       )
        )





atsp <- read.delim("E:/PhD Thesis/redo PhD/AtSubP_results.txt", header = F)
tarp <- read.delim("E:/PhD Thesis/redo PhD/TargetP_analysis.txt", header = F)
mitochondrial <- str_split_fixed(c(atsp$V1[atsp$V2 == "Mitochondrion"], tarp$V1[tarp$V2 == "M"]),"\\.",2)[,1] %>% unique()
plastidial <- str_split_fixed(c(atsp$V1[atsp$V2 == "Chloroplast"], tarp$V1[tarp$V2 == "C"]),"\\.",2)[,1] %>% unique()
library(org.At.tair.db)
mitop <- AnnotationDbi::select(org.At.tair.db, keys = "GO:0005739", columns = "TAIR", keytype = "GO")
chlorop <- AnnotationDbi::select(org.At.tair.db, keys = "GO:0009507", columns = "TAIR", keytype = "GO")
mitochondrial <- unique(c(mitochondrial, mitop$TAIR))
plastidial<- unique(c(plastidial, chlorop$TAIR))


mpnc <- exp_true %>% filter(locus %in% mitochondrial) %>% mutate(family = genefam$Gene_Family[match(locus, genefam$Genomic_Locus_Tag)],
                                                                 sdesc = fdesc$Short_description[match(locus, fdesc$Model_name)])
mpnc <- mpnc[-grep("ATMG", mpnc$locus),]
mpnc <- mpnc[-grep("AT2G07", mpnc$locus),]
ppnc <- exp_true %>% filter(locus %in% plastidial)%>% mutate(family = genefam$Gene_Family[match(locus, genefam$Genomic_Locus_Tag)],
                                                             sdesc = fdesc$Short_description[match(locus, fdesc$Model_name)])

library(ComplexHeatmap)
col_fun5 = circlize::colorRamp2(c(6,3,1.5,0,-1.5,-3), c("orange", "gold", "chartreuse", "black", "dodgerblue", "firebrick1"))

Heatmap(as.matrix(mpnc[,c(2:9)]), col = col_fun5, name = "LFC",row_split = mpnc$family,
        row_labels = paste0(mpnc$locus, " \ ",mpnc$symbol), row_gap = unit(0, "mm"),
        row_names_gp = gpar(fontsize = 0),column_names_gp = gpar(fontsize = 15), 
        row_title_rot = 0, row_title_gp = gpar(fontsize = 0),
        right_annotation = HeatmapAnnotation(
        `Family` = mpnc$family,
         which = "row", col = list(`Family`= setNames( viridis::turbo(27), unique(mpnc$family)[-2]))))

Heatmap(as.matrix(ppnc[,c(2:9)]), col = col_fun5, name = "LFC",row_split = ppnc$family,
        row_labels = paste0(ppnc$locus, " \ ",ppnc$symbol), row_gap = unit(0, "mm"),
        row_names_gp = gpar(fontsize = 0),column_names_gp = gpar(fontsize = 15), 
        row_title_rot = 0, row_title_gp = gpar(fontsize = 0),
        right_annotation = HeatmapAnnotation(
          `Family` = ppnc$family,
          which = "row", col = list(`Family`= setNames( viridis::turbo(38), unique(ppnc$family)[-4]))))

         

library(openxlsx)
# write.xlsx(mpnc, "mt-nuclear-encoded.xlsx")
# write.xlsx(ppnc, "pt-nuclear-encoded.xlsx")

library(venn)
uptest <- list(
               OEM15 = upregulated$OEM15[upregulated$OEM15 %in% mpnc$locus],
               OEP12 = upregulated$OEP12[upregulated$OEP12 %in% mpnc$locus],
               `Tmp-P1` = upregulated$`Tmp-P1`[upregulated$`Tmp-P1` %in% mpnc$locus],
               `Tmp-M3` = upregulated$`Tmp-M3`[upregulated$`Tmp-M3` %in% mpnc$locus],
               rpotmp = upregulated$rpotmp[upregulated$rpotmp %in% mpnc$locus])

dotest <- list(
  `OEM15` = downregulated$OEM15[downregulated$OEM15 %in% mpnc$locus],
               `OEP12` = downregulated$OEP12[downregulated$OEP12 %in% mpnc$locus],
  `Tmp-P1` = downregulated$`Tmp-P1`[downregulated$`Tmp-P1` %in% mpnc$locus],
  `Tmp-M3` = downregulated$`Tmp-M3`[downregulated$`Tmp-M3` %in% mpnc$locus],
  rpotmp = downregulated$rpotmp[downregulated$rpotmp %in% mpnc$locus])

uptest <- list(
  OEM15 = upregulated$OEM15[upregulated$OEM15 %in% ppnc$locus],
  OEP12 = upregulated$OEP12[upregulated$OEP12 %in% ppnc$locus],
  `Tmp-P1` = upregulated$`Tmp-P1`[upregulated$`Tmp-P1` %in% ppnc$locus],
  `Tmp-M3` = upregulated$`Tmp-M3`[upregulated$`Tmp-M3` %in% ppnc$locus],
  rpotmp = upregulated$rpotmp[upregulated$rpotmp %in% ppnc$locus])

dotest <- list(
  `OEM15` = downregulated$OEM15[downregulated$OEM15 %in% ppnc$locus],
  `OEP12` = downregulated$OEP12[downregulated$OEP12 %in% ppnc$locus],
  `Tmp-P1` = downregulated$`Tmp-P1`[downregulated$`Tmp-P1` %in% ppnc$locus],
  `Tmp-M3` = downregulated$`Tmp-M3`[downregulated$`Tmp-M3` %in% ppnc$locus],
  rpotmp = downregulated$rpotmp[downregulated$rpotmp %in% ppnc$locus])

p1 <- venn(uptest, zcolor = "style", opacity = 0.15, ggplot  = T, box = F, ilcs=1.2, ilabels = "counts")
p2 <- venn(dotest, zcolor = "style", opacity = 0.15, ggplot  = T, box = F,ilcs=1.2, ilabels = "counts")
ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2)

table(mpnc$family)
unique(mpnc$family)


m15up <- setdiff(uptest$OEM15, c(uptest$OEP12, uptest$`Tmp-P1`, uptest$rpotmp))
m15dn <- setdiff(dotest$OEM15, c(dotest$OEP12, dotest$`Tmp-P1`, dotest$rpotmp))
noup <- setdiff(intersect(uptest$`Tmp-P1`, uptest$rpotmp), c(uptest$OEP12, uptest$OEM15))
nodn <- setdiff(intersect(dotest$`Tmp-P1`, dotest$rpotmp), c(dotest$OEP12, dotest$OEM15))

fg <- tibble(RRG = c(rep("excess mt-RPOTM", 70),rep("No mt-RPOTmp", 25)),
             Regulation = c(rep("enhanced",17), rep("suppressed",53), rep("enhanced",20), rep("suppressed",5)),
             Gene = c(m15up, m15dn, noup, nodn)) %>% mutate(symbol = exp_true$symbol[match(Gene, exp_true$locus)],
                                                            sdesc = fdesc$Short_description[match(Gene, fdesc$Model_name)],
                                                            family = genefam$Gene_Family[match(Gene, genefam$Genomic_Locus_Tag)])


p12up <- setdiff(uptest$OEP12, c(uptest$OEM15, uptest$`Tmp-P1`, uptest$rpotmp))
p12dn <- setdiff(dotest$OEP12, c(dotest$OEM15, dotest$`Tmp-P1`, dotest$rpotmp))
nopup <- setdiff(intersect(uptest$`Tmp-M3`, uptest$rpotmp), c(uptest$OEP12, uptest$OEM15))
nopdn <- setdiff(intersect(dotest$`Tmp-M3`, dotest$rpotmp), c(dotest$OEP12, dotest$OEM15))

fg2 <- tibble(RRG = c(rep("excess pt-RPOTmp", 47),rep("No pt-RPOTmp", 11)),
             Regulation = c(rep("enhanced",23), rep("suppressed",24), rep("enhanced",8), rep("suppressed",3)),
             Gene = c(p12up, p12dn, nopup, nopdn)) %>% mutate(symbol = exp_true$symbol[match(Gene, exp_true$locus)],
                                                            sdesc = fdesc$Short_description[match(Gene, fdesc$Model_name)],
                                                            family = genefam$Gene_Family[match(Gene, genefam$Genomic_Locus_Tag)])
# write.xlsx(fg, "rrg1.xlsx")
# write.xlsx(fg2, "rrg2.xlsx")
fg <- read.xlsx("rrg1.xlsx")
fg2 <- read.xlsx("rrg2.xlsx")

library(igraph)
library(ggraph)
library(scales)
rrg1 <- graph_from_data_frame(fg[,c(1,3,2)])
V(rrg1)$family = fg$family[match(V(rrg1)$name, fg$Gene)]
V(rrg1)$symbol = str_split_fixed(fg$symbol[match(V(rrg1)$name, fg$Gene)], ";",2)[,1]

ggraph(rrg1, layout = "stress")+
  geom_edge_fan(alpha = 0.3, aes(col = Regulation), 
                arrow = arrow(angle = 10, type = "closed", length = unit(5,"mm")), end_cap = circle(3, "mm"))+
  geom_node_point(aes(col = family), size = 3, alpha = 0.6)+
  geom_node_text(aes(label = symbol), size =2, repel = F)+theme_graph()+
  scale_color_manual("?????????", values = setNames(viridis::turbo(12), unique(fg$family)[-1]))+
  scale_edge_color_manual("???????????? ?????", values = setNames(c("green", "firebrick1"), c("enhanced", "suppressed")))

library(clusterProfiler)
enr1 <- enricher(gene = fg$Gene[fg$RRG == "No mt-RPOTmp" ], pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = data.frame(TERM = tfdb$V1, GENE = tfdb$V3), minGSSize = 3, maxGSSize = 10000)
enr1 <- as_tibble(enr1)
t2s <- AnnotationDbi::select(org.At.tair.db, keys = enr1$ID, keytype = "TAIR", columns = "SYMBOL") %>% stats::aggregate(SYMBOL~TAIR, FUN = function(x){paste(x, collapse = ";")}  )
enr1 <- enr1 %>% mutate(symbol = t2s$SYMBOL[match(ID, t2s$TAIR)],
                        p.sig = -log10(p.adjust))
enr1$symbol[is.na(enr1$symbol)==T] <- "NAC-domain"
ggplot(enr1)+
  geom_point(aes(x = p.sig, y = toupper(str_split_fixed(symbol, ";",2)[,1]), size = Count), col = "dodgerblue")+
  scale_size_continuous("????? ?????")+
  ylab("???????????????? ??????")+xlab("?????????? P")+theme_light()+theme(axis.title.x = element_text(size = 15),
                                                                        axis.title.y = element_text(size = 15),
                                                                        axis.text.x = element_text(size = 10),
                                                                        axis.text.y = element_text(size = 15))


rrg2 <- graph_from_data_frame(fg2[,c(1,3,2)])
V(rrg2)$family = fg2$family[match(V(rrg2)$name, fg2$Gene)]
V(rrg2)$symbol = str_split_fixed(fg2$symbol[match(V(rrg2)$name, fg2$Gene)], ";",2)[,1]

ggraph(rrg2, layout = "stress")+
  geom_edge_fan(alpha = 0.3, aes(col = Regulation), 
                arrow = arrow(angle = 10, type = "closed", length = unit(5,"mm")), end_cap = circle(3, "mm"))+
  geom_node_point(aes(col = family), size = 3, alpha = 0.6)+
  geom_node_text(aes(label = symbol), size =2, repel = F)+theme_graph()+
  scale_color_manual("?????????", values = setNames(viridis::turbo(23), unique(fg2$family)[-1]))+
  scale_edge_color_manual("???????????? ?????", values = setNames(c("green", "firebrick1"), c("enhanced", "suppressed")))

library(clusterProfiler)
enr2 <- enricher(gene = fg2$Gene[fg2$RRG == "excess pt-RPOTmp" ], pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = data.frame(TERM = tfdb$V1, GENE = tfdb$V3), minGSSize = 3, maxGSSize = 10000)
enr2 <- as_tibble(enr2)
t2s <- AnnotationDbi::select(org.At.tair.db, keys = enr2$ID, keytype = "TAIR", columns = "SYMBOL") %>% stats::aggregate(SYMBOL~TAIR, FUN = function(x){paste(x, collapse = ";")}  )
enr2 <- enr2 %>% mutate(symbol = t2s$SYMBOL[match(ID, t2s$TAIR)],
                        p.sig = -log10(p.adjust))
# no enrichment.

# van Aken's marker genes
mg <- c("AT2G41730", "AT2G03760", "AT2G04050", "AT2G21640", "AT5G59820", 
        "AT5G09570", "AT1G72910")
mge <- exp_true %>% filter(locus %in% mg)
mge$pertirbations = c("plastid", "mt and pt", "mt and pt", "mt and pt", "mitoch", "mt and pt")
library(ComplexHeatmap)
col_fun5 = circlize::colorRamp2(c(6,3,1.5,0,-1.5,-3), c("orange", "gold", "chartreuse", "black", "dodgerblue", "firebrick1"))

Heatmap(as.matrix(mge[,c(2:9)]), col = col_fun5, name = "LFC",row_split = mge$pertirbations,
        row_labels = paste0(mge$locus, " \ ",str_split_fixed(mge$symbol, ";",2)[,1]), row_gap = unit(0, "mm"),
        row_names_gp = gpar(fontsize =12 ),column_names_gp = gpar(fontsize = 15), 
        row_title_rot = 0, row_title_gp = gpar(fontsize = 0),
        right_annotation = HeatmapAnnotation(
          `?????????` = mge$pertirbations,
          which = "row", col = list(`?????????`= setNames( c("chartreuse", "slateblue", "tomato"), unique(mge$pertirbations)))))





get_interactors <- function(x){
 res <-  unique(c(string[string$TAIR8A %in% x | string$TAIR8B %in% x, c(1,2)][,1],string[string$TAIR8A %in% x | string$TAIR8B %in% x, c(1,2)][,2])) 
return(res)}

#PPR and TPR
PPR <- Expressionne[unique(c(grep("PPR", Expressionne$Description),
         grep("PPR", Expressionne$description),
         grep("PPR", Expressionne$DESCRIPTION),
         grep("PPR", Expressionne$gene_family),
         grep("PPR", Expressionne$gene_subfamily),
         grep("PPR", Expressionne$gene_name2),
         grep("PPR", Expressionne$gene_fam_desc),
         grep("PPR", Expressionne$description2),
         grep("PPR", Expressionne$tair_comp_desc),
         grep("TPR", Expressionne$Description),
         grep("TPR", Expressionne$description),
         grep("TPR", Expressionne$DESCRIPTION),
         grep("TPR", Expressionne$gene_family),
         grep("TPR", Expressionne$gene_subfamily),
         grep("TPR", Expressionne$gene_name2),
         grep("TPR", Expressionne$gene_fam_desc),
         grep("TPR", Expressionne$description2),
         grep("TPR", Expressionne$tair_comp_desc),
         grep("tricopeptid", Expressionne$Description),
         grep("tricopeptid", Expressionne$description),
         grep("tricopeptid", Expressionne$DESCRIPTION),
         grep("tricopeptid", Expressionne$gene_family),
         grep("tricopeptid", Expressionne$gene_subfamily),
         grep("tricopeptid", Expressionne$gene_name2),
         grep("tricopeptid", Expressionne$gene_fam_desc),
         grep("tricopeptid", Expressionne$description2),
         grep("tricopeptid", Expressionne$tair_comp_desc))),]

pppr <- exp_true %>% filter(locus %in% PPR$locus)
t2s <- AnnotationDbi::select(org.At.tair.db, keys = pppr$locus, keytype = "TAIR", columns = "SYMBOL") %>% stats::aggregate(SYMBOL~TAIR, FUN = function(x){paste(x, collapse = ";")}  )
pppr$symbol <- t2s$SYMBOL[match(pppr$locus, t2s$TAIR)]


library(ComplexHeatmap)
pdf <- readxl::read_xls("ppr1.xls") %>% filter(`PPR=AGI`==T)
pdf2 <- readxl::read_xls("ppr2.xls") %>% mutate(TAIR = toupper(str_split_fixed(pdf$`best AGI match`[match(`PPR model`, pdf$`PPR code`)], "\\.",2)[,1]))

liu <- read.xlsx("liuPPR.xlsx") %>% mutate(Gene.ID = str_split_fixed(Gene.ID, "\\.",2)[,1])
pppr3 <- exp_true %>% filter(locus %in% liu$Gene.ID) %>% mutate(subclass = liu$PPR.type[match(locus, liu$Gene.ID)])


pppr <- exp_true[exp_true$locus %in% pdf2$TAIR,] %>% mutate(subclass = pdf2$subclass[match(locus, pdf2$TAIR)])

pppr2 <- exp_true %>% filter(locus %in% PPR$locus) %>% filter(locus %in% setdiff(locus, pppr$locus))
t2s <- AnnotationDbi::select(org.At.tair.db, keys = pppr2$locus, keytype = "TAIR", columns = "SYMBOL") %>% stats::aggregate(SYMBOL~TAIR, FUN = function(x){paste(x, collapse = ";")}  )
pppr2$symbol <- t2s$SYMBOL[match(ppp2r$locus, t2s$TAIR)]
pppr2$subclass <- "Other"
PPR <- rbind(pppr,pppr2)
PPR$locus %in% liu$Gene.ID

PPR$subclass <- liu$PPR.type[match(PPR$locus, liu$Gene.ID)]
PPR$subclass[is.na(PPR$subclass)==T] <- "Other"

suba <- read.csv("suba5/19_8_2022_suba5_summary.csv") %>% select(c(1,3,40)) %>% mutate(locus = str_split_fixed(locus, "\\.",2)[,1]) %>% as_tibble() %>% distinct()

suba2 <- aggregate(suba, location_consensus~locus, FUN = function(x){paste(x, collapse = ";")})
suba2$description <- suba$description[match(suba2$locus, suba$locus)]

PPR$localization <- suba2$location_consensus[match(PPR$locus, suba2$locus)]
PPR$symbol[is.na(PPR$symbol)==T] <- "-"

Heatmap(as.matrix(PPR[,c(2:9)]), col = col_fun5, name = "LFC",row_split = PPR$subclass,
        row_labels = paste(PPR$locus, str_split_fixed(PPR$symbol, ";",2)[,1]), row_gap = unit(0, "mm"),
        row_names_gp = gpar(fontsize =7 ),column_names_gp = gpar(fontsize = 15), 
        row_title_rot = 0, row_title_gp = gpar(fontsize = 0),
        right_annotation = HeatmapAnnotation(
          `????????` = PPR$subclass,
          `???????????` = PPR$localization,
          which = "row", col = list(`????????`= setNames( c("chartreuse", "slateblue", "tomato", "dodgerblue", "gold", "firebrick1"), unique(PPR$subclass)),
          `???????????` = setNames(viridis::turbo(9), unique(PPR$localization)[-9])  )))



#Gene ontology
upGO <- list()
doGO <- list()
pb <- txtProgressBar(max = 8, width = 100, style = 3)
for(i in 1:8){
upGO[[i]] <- enrichGO(upregulated[[i]], OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP")
doGO[[i]] <- enrichGO(downregulated[[i]], OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP")
setTxtProgressBar(pb,i)
}
close(pb)
names(upGO) <- names(upregulated)
names(doGO) <- names(downregulated)

for(i in 1:8){
  upGO[[i]] <- as_tibble(upGO[[i]])
  doGO[[i]] <- as_tibble(doGO[[i]])
}
upGO[[1]]$p.adjust
for(i in 1:8){
  upGO[[i]]$P.significance <- (-1)*log10(upGO[[i]]$p.adjust) 
  doGO[[i]]$P.significance <- (-1)*log10(doGO[[i]]$p.adjust) 
}
for(i in 1:8){
  upGO[[i]]$line <- names(upGO)[i] 
  doGO[[i]]$line <- names(doGO)[i] 
}

updf <- do.call("rbind", upGO)
dodf <- do.call("rbind", doGO)
updf$regulation <- "upregulated"
dodf$regulation <- "downregulated"

rrdf <- rbind(updf,dodf)

rrdf <- pivot_wider(rrdf, names_from = "line",
                    values_from = "P.significance")

updf <- pivot_wider(updf, names_from = "line",
                    values_from = "P.significance")

dodf <- pivot_wider(dodf, names_from = "line",
                    values_from = "P.significance")

updf[is.na(updf)==T] <- 0
dodf[is.na(dodf)==T] <- 0
rrdf[is.na(rrdf)==T] <- 0
summary(updf)
library(ComplexHeatmap)
library(circlize)
col_fun7 = colorRamp2(c(10, 7.5,5,2.5,0), c("red", "orange", "yellow", "green", "black"))

Heatmap(as.matrix(rrdf[,c(11:16)]), col = col_fun7, row_split = rrdf$regulation,
        row_labels = rrdf$Description, 
        row_names_gp = gpar(fontsize = 6), row_title_rot = 0, row_title_gp = gpar(fontsize = 7),  
        right_annotation = HeatmapAnnotation(`Count` = rrdf$Count,
         which = "row"
        )
)

# co-expression

#table

cemdf <- tibble(module = paste0(rep("M",5), c(1:5)),
                Probes = c(2908, 1267, 102, 72, 39 ),
                AGI = c(3149,1312,114,72,39)) %>% pivot_longer(cols = 2:3, names_to = "Parameter", values_to = "Value")
cemdf$Parameter <- factor(cemdf$Parameter, levels = c("Probes", "AGI"))
ggplot(cemdf)+
  geom_col(aes(x = module, y = Value, fill = Parameter), col = "black")+
  geom_text(aes(x = module, y = Value + 500, label = Value),angle = 90, fontface = "bold", size = 7)+facet_wrap(~Parameter)+
  scale_fill_manual("", values = setNames(c("tomato", "dodgerblue"), c("Probes", "AGI")))+theme_light()+ylim(c(0,4000))+
  theme(strip.background = element_rect(fill = "white", colour =  "grey"),
        strip.text = element_text(colour = "black", face = "bold", size = 15),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 15),
        legend.position = "none")+ylab("")+xlab("?????? ??-??????????")


# synthesis
setwd("D:/git/phd-thesis/Dissertation!/")
setwd("E:/ivg-git/phd/phd-thesis/Dissertation!/")
library(igraph)
library(ggraph)
library(openxlsx)
gppi <- readRDS("gppi.rds")
degs <- read.xlsx("Transcriptome.merged.xlsx")
mse <- readRDS("mst_edges.RDS")
msv <- readRDS("mst_ver.RDS")

mst <- graph_from_data_frame(mse, vertices = msv, directed = F)

library(ggraph)
mst
V(mst)$DE <- ifelse(V(mst)$name %in% degs$ProbeName, "DEG","non-DEG")
ggraph(mst, "backbone")+
  geom_edge_link(aes(colour = color), alpha = 0.6, size = 0.7)+geom_node_point(aes(colour = DE),alpha=0.8, size = 0.4)+theme_graph()+
  scale_edge_color_manual("Correlation", values = setNames(c("dodgerblue", "firebrick1"),c("dodgerblue", "firebrick1")))+
  scale_color_manual("DE", values = setNames(c("firebrick1", "darkgrey"),c("DEG", "non-DEG")))

m15st <- subgraph(mst, which(V(mst)$OEM15 !=0)) 
p12st <- subgraph(mst, which(V(mst)$OEP12 !=0))

gppi2 <- igraph::delete_edges(graph = gppi, E(gppi)[which(E(gppi)$lsMR==0)])

gm15 <- purrr::reduce(purrr::map(make_ego_graph(gppi2, order = 1, nodes = which(V(gppi2)$name %in% V(m15st)$TAIR)),
                                 .f = function(x){as_data_frame(x,"edges")}), .f= rbind) %>% graph_from_data_frame(directed = F)
gp12 <- purrr::reduce(purrr::map(make_ego_graph(gppi2, order = 1, nodes = which(V(gppi2)$name %in% V(p12st)$TAIR)),
                                 .f = function(x){as_data_frame(x,"edges")}), .f= rbind)%>% graph_from_data_frame(directed = F)


gm15 <- igraph::simplify(gm15, edge.attr.comb = "first")
E(gm15)$color <- NA_character_
E(gm15)$color[which(E(gm15)$lsMR > 0)] <- "positive"
E(gm15)$color[which(E(gm15)$lsMR < 0)] <- "negative"
E(gm15)$color[which(E(gm15)$lsMR == 0)] <- "n/d"

gp12 <- igraph::simplify(gp12, edge.attr.comb = "first")
E(gp12)$color <- NA_character_
E(gp12)$color[which(E(gp12)$lsMR > 0)] <- "positive"
E(gp12)$color[which(E(gp12)$lsMR < 0)] <- "negative"
E(gp12)$color[which(E(gp12)$lsMR == 0)] <- "n/d"

select_top_n_nodes <- function(g,n){
  deg <-  degree(g)
  bet <- betweenness(g)
  clo <- closeness(g)
  le <- local_efficiency(g)
  hc <- harmonic_centrality(g)
  df <- tibble(deg = deg,
               tair = names(deg))
  db <- tibble(bet = bet,
               tair = names(bet))
  dc <- tibble(clo = clo,
               tair = names(clo))
  dl <- tibble(le = le,
               tair = names(le))
  dh <- tibble(hc = hc,
               tair = names(hc))
  df <- df[order(df$deg, decreasing = T),]
  db <- db[order(db$bet, decreasing = T),]
  dc <- dc[order(dc$clo, decreasing = T),]
  dl <- dl[order(dl$le, decreasing = T),]
  dh <- dh[order(dh$hc, decreasing = T),]
  df <- df[c(1:n),]
  db <- db[c(1:n),]
  dc <- dc[c(1:n),]
  dl <- dl[c(1:n),]
  dh <- dh[c(1:n),]
  df <- pivot_longer(df,cols = 1, names_to = "centrality", values_to = "value")
  db <- pivot_longer(db,cols = 1, names_to = "centrality", values_to = "value")
  dc <- pivot_longer(dc,cols = 1, names_to = "centrality", values_to = "value")
  dl <- pivot_longer(dl,cols = 1, names_to = "centrality", values_to = "value")
  dh <- pivot_longer(dh,cols = 1, names_to = "centrality", values_to = "value")
  result <- rbind(df,db,dc,dl,dh) %>% filter(value > 0) %>% pivot_wider(names_from = "centrality", values_from = "value")
  return(result)
}

gm15 <- delete_edges(gm15, which(E(gm15)$color == "n/d"))
gm15 <- delete_vertices(gm15, which(degree(gm15)==0))
gp12 <- delete_edges(gp12, which(E(gp12)$color == "n/d"))
gp12 <- delete_vertices(gp12, which(degree(gp12)==0))

check <- intersect(V(gm15)$name, V(gp12)$name)
dput(check)
check2 <- setdiff(check, V(mst)$TAIR)
dput(check2)


V(gm15)$goi <- NA_character_
V(gp12)$goi <- NA_character_

V(gm15)$goi[which(V(gm15)$name %in% check2)] <- "PPI"
V(gp12)$goi[which(V(gp12)$name %in% check2)] <- "PPI"

check3 <- intersect(V(m15st)$TAIR, V(p12st)$TAIR)

V(gm15)$goi[which(V(gm15)$name %in% check3)] <- "MST"
V(gp12)$goi[which(V(gp12)$name %in% check3)] <- "MST"


gm15 <- igraph::simplify(gm15, edge.attr.comb = "first")
gp12 <- igraph::simplify(gp12, edge.attr.comb = "first")

library(tidyverse)
mse2 <- mse %>% mutate(from = msv$TAIR[match(from, msv$name)],
                       to = msv$TAIR[match(to, msv$name)], color = "Correlation")
mst2 <- graph_from_data_frame(mse2, directed = F)
fm15 <- subgraph(mst2, V(mst2)[which(V(mst2)$name %in% V(gm15)$name)])
fp12 <- subgraph(mst2, V(mst2)[which(V(mst2)$name %in% V(gp12)$name)])

fp12 <- igraph::simplify(fp12)
fm15 <- igraph::simplify(fm15)


fm15 <- subgraph(fm15,which(degree(fm15)>0))
fp12 <- subgraph(fp12,which(degree(fp12)>0))



gm151 <- add_edges(gm15,as_edgelist(fm15))
E(gm151)$color[is.na(E(gm151)$color)==T] <- "Correlated"

gp121 <- add_edges(gp12,as_edgelist(fp12))
E(gp121)$color[is.na(E(gp121)$color)==T] <- "Correlated"



p1 <- ggraph(gm151, layout = "fr")+geom_edge_fan(aes(color = color), alpha = 0.6, strength = 3)+geom_node_point(aes(colour = goi), alpha = 0.8)+
  scale_edge_color_manual("лоВР", values = setNames(c("gold2", "firebrick1", "dodgerblue"), c("Correlated", "negative", "positive")))+
  scale_color_manual("Общие узлы", values = setNames(c("chartreuse", "deeppink"),c( "PPI", "MST") ))+
  theme_graph()

p2 <- ggraph(gp121, layout = "fr")+geom_edge_fan(aes(color = color), alpha = 0.6, strength = 3)+geom_node_point(aes(colour = goi), alpha = 0.8)+
  scale_edge_color_manual("лоВР", values = setNames(c("gold2", "firebrick1", "dodgerblue"), c("Correlated", "negative", "positive")))+
  scale_color_manual("Общие узлы", values = setNames(c("chartreuse", "deeppink"),c( "PPI", "MST") ))+
  theme_graph()


ggpubr::ggarrange(p1,p2, ncol = 2, labels = "AUTO", font.label = list(size = 25), common.legend = T)

checked <- read.csv("summary3.csv", sep = " ")

#let's add more - the TFtargets!

tft <- readRDS("tft.RDS")
tfdb <- readRDS("tfdb.rds")
tfddb <- readRDS("tfddb.rds")

degs <- read.xlsx("Transcriptome.merged.xlsx")

tfm <- tfddb %>% filter(from %in% V(gm151)$name, to %in% degs$locus[abs(degs$M15)>0]) 
tfp <- tfddb %>% filter(from %in% V(gp121)$name, to %in% degs$locus[abs(degs$P12)>0]) 
tftargets <- intersect(tfm$to, tfp$to)

shared_targets <- exp_true %>% filter(locus %in% tftargets) %>% select(c(1:3, 10))

# 64 genes - targets of TFs found!

shared_targets <- shared_targets %>% mutate(localization = suba2$location_consensus[match(locus, suba2$locus)],
                          desc = suba2$description[match(locus, suba2$locus)],
                          family = genefam$Gene_Family[match(locus, genefam$Genomic_Locus_Tag)],
                          func = genefam$Protein_Function[match(locus, genefam$Genomic_Locus_Tag)])

write.xlsx(shared_targets, "shared_targets.xlsx")

tfm2 <- tfm %>% filter(to %in% tftargets) %>% distinct() 
tfp2 <- tfp %>% filter(to %in% tftargets) %>% distinct()

tfm2 <- graph_from_data_frame(tfm2, directed = F) %>% simplify()
tfp2 <- graph_from_data_frame(tfp2, directed = F) %>% simplify()

gm152 <- igraph::union(gm151, tfm2)
E(gm152)$color[is.na(E(gm152)$color)==T] <- "TF-target"

gp122 <- igraph::union(gp121, tfp2)
E(gp122)$color[is.na(E(gp122)$color)==T] <- "TF-target"

V(gm152)$goi[V(gm152)$name %in% tftargets] <- "TF_target"
V(gp122)$goi[V(gp122)$name %in% tftargets] <- "TF_target"

# saveRDS(gm152, "oem15_combined_network.rds")
# saveRDS(gp122, "oep12_combined_network.rds")

gp122 <- igraph::simplify(gp122, remove.multiple = F, remove.loops = T)

p1 <- ggraph(gm152, layout = "fr")+geom_edge_fan(aes(color = color, linetype = color), alpha = 0.6, strength = 3)+geom_node_point(aes(colour = goi), alpha = 0.8)+
  scale_edge_color_manual("лоВР", values = setNames(c("lightgrey", "gold2", "firebrick1", "dodgerblue"), c("TF-target", "Correlated", "negative", "positive")))+
  scale_edge_linetype_manual("лоВР", values = setNames(c("dashed", "solid", "solid", "solid"), c("TF-target", "Correlated", "negative", "positive")))+
  scale_color_manual("Общие узлы", values = setNames(c("chartreuse", "deeppink", "blue"),c( "PPI", "MST", "TF_target") ))+
  theme_graph()

p2 <- ggraph(gp122, layout = "fr")+geom_edge_fan(aes(color = color, linetype = color), alpha = 0.6, strength = 3)+geom_node_point(aes(colour = goi), alpha = 0.8)+
  scale_edge_color_manual("лоВР", values = setNames(c("lightgrey", "gold2", "firebrick1", "dodgerblue"), c("TF-target", "Correlated", "negative", "positive")))+
  scale_edge_linetype_manual("лоВР", values = setNames(c("dashed", "solid", "solid", "solid"), c("TF-target", "Correlated", "negative", "positive")))+
  scale_color_manual("Общие узлы", values = setNames(c("chartreuse", "deeppink", "blue"),c( "PPI", "MST", "TF_target") ))+
  theme_graph()

ggpubr::ggarrange(p1,p2, ncol = 2, labels = "AUTO", font.label = list(size = 25), common.legend = T)

#summary 3 - shared PPI-nodes,
#shared TF targets:
tftargets <- c("AT3G02100", "AT4G25200", "AT2G11405", "AT5G55020", "AT1G60095", 
  "AT1G61820", "AT2G02490", "AT2G05370", "AT5G02600", "AT5G10200", 
  "AT5G37910", "AT5G48930", "AT1G10070", "AT1G11610", "AT1G19900", 
  "AT1G26773", "AT1G31072", "AT1G35270", "AT1G44085", "AT1G47300", 
  "AT1G48640", "AT1G56300", "AT1G67265", "AT1G76930", "AT1G80390", 
  "AT2G14870", "AT2G20380", "AT2G22750", "AT2G24660", "AT2G37700", 
  "AT2G47550", "AT2G48060", "AT3G16680", "AT3G19140", "AT3G29370", 
  "AT3G50010", "AT3G54802", "AT3G61270", "AT4G01975", "AT4G02830", 
  "AT4G04330", "AT4G08850", "AT4G14465", "AT4G23560", "AT4G31970", 
  "AT4G35280", "AT4G36370", "AT4G36380", "AT4G36580", "AT4G37900", 
  "AT5G07500", "AT5G09795", "AT5G28673", "AT5G35550", "AT5G45090", 
  "AT5G45830", "AT3G06090", "AT5G55430", "AT1G08890", "AT3G02885", 
  "AT3G25490", "AT4G28040", "AT4G33905", "AT5G63085")

degs %>% filter(locus %in% tftargets) %>% View()

library(clusterProfiler)
library(org.At.tair.db)
ego <- enrichGO(fhm$locus, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "ALL", maxGSSize = 10000, minGSSize = 3)
ego

library(tidyverse)
shared <- read.csv("summary3.csv", sep = " ") %>% as_tibble()

degs <- read.xlsx("Transcriptome.merged.xlsx")

fhm <- exp_true %>% filter(locus %in% c(tftargets)) %>% mutate(type = NA_character_) %>% filter(rpotmp ==0, `Tmp-P1` ==0)
fhm$type[fhm$locus %in% shared$input] <- "???"
fhm$type[fhm$locus %in% tftargets] <- "?????? ??"


library(ComplexHeatmap)
col_fun5 = circlize::colorRamp2(c(6,3,1.5,0,-1.5,-3), c("orange", "gold", "chartreuse", "black", "dodgerblue", "firebrick1"))
library(org.At.tair.db)
AnnotationDbi::select(org.At.tair.db, keys = tftargets, keytype = "TAIR", columns = "SYMBOL")

tfta <- read.csv("tftargets.csv", sep = " ")

Heatmap(as.matrix(fhm[,c(2:9)]), col = col_fun5, name = "LFC",row_split = fhm$type,
        row_labels = paste0(fhm$locus, " \ ",str_split_fixed(fhm$symbol, ";",2)[,1]), row_gap = unit(0, "mm"),
        row_names_gp = gpar(fontsize =12 ),column_names_gp = gpar(fontsize = 15), 
        row_title_rot = 0, row_title_gp = gpar(fontsize = 15),
        # right_annotation = HeatmapAnnotation(`?????????` = mge$pertirbations,
          # which = "row", col = list(`?????????`= setNames( c("chartreuse", "slateblue", "tomato"), unique(mge$pertirbations))))
        )

# dark germination

library(tidyverse)
library(openxlsx)
library(rstatix)
setwd("D:/git/phd-thesis/Dissertation!/")
ld <- read.xlsx("lightdark.xlsx") 
ld$Line <- factor(ld$Line, levels = c("Col-0", "rpotmp", "abi4", "OEM1", "OEM15", "OEP12", "OEP15"))
ld$Cond <- factor(ld$Cond, levels = c("Light", "Dark"))

ggplot(ld)+
  geom_point(aes(x = Line, y = Germ))+
  geom_errorbar(aes(x = Line, ymin = Germ-SD, ymax = Germ + SD), width = 0.4)+
  geom_text(aes(x = Line, y = Germ+10, label = SigSelf), col = "dodgerblue", size = 10)+
  geom_text(aes(x = Line, y = Germ+10, label = SigCol), col = "firebrick1", size = 10)+
facet_wrap(~Cond)+theme_bw()+theme(
    axis.text.x = element_text(angle = 90, size = 15, hjust = 1, vjust = 0.5),
    strip.background = element_rect(fill = "white", colour =  "black"),
    strip.text = element_text(face = "bold", colour = "black")
  )
#salt

nacl <- read.xlsx("Salt germination.xlsx") %>% pivot_longer(cols = 2:17, names_to = "Cond", values_to = "Value") 

nacl$Cond <- gsub("8-9", "9", nacl$Cond)
nacl$Cond <- gsub("\\.", "", nacl$Cond)

nacl <- nacl %>% 
  separate_wider_delim(cols = 2, delim = "-", names = c("Cond", "Day"))
nacl$SP <- gsub("Col-M", "OEM", nacl$SP)
nacl$SP <- gsub("Col-P", "OEP", nacl$SP)
nacl$SP <- gsub("TP-1", "Tmp-P1", nacl$SP)
nacl2 <- nacl %>% filter(Cond %in% c("Ctrl", "150mM")) %>% filter(SP %in% c("Сol-0", "OEM1", "OEM14", "OEM15", "OEM20", "OEP8", "OEP12", "OEP15", "OEP18", "rpotmp")) %>% mutate(Value = Value*100 )
nacl2$SP <- factor(nacl2$SP, levels = c("Сol-0", "OEM1", "OEM14", "OEM15", "OEM20", "OEP8", "OEP12", "OEP15", "OEP18", "rpotmp"))

nacl150 <- nacl2 %>% filter(Cond == "150mM", Day %in% c(4,7), is.na(Value)==F) 
naclC <- nacl2 %>% filter(Cond == "Ctrl", Day %in% c(2,3,4,7), is.na(Value)==F) 

nacl150 %>% group_by(Day) %>% wilcox_test(Value~SP, ref.group = "Сol-0", p.adjust.method = "none")
nacl150$Sig <- NA_character_
nacl150$Sig[nacl150$SP == "rpotmp"] <- "*"
sd()
ctrl <- naclC %>% group_by(Day, SP) %>% reframe(Value = mean(Value, na.rm = T))
csd <- naclC %>% group_by(Day, SP) %>% mutate(sd = sd(Value, na.rm = T))
ctrl <- left_join(ctrl, csd[,-4]) %>% distinct()


p1 <- ggplot(ctrl)+
  geom_col(aes(x = SP, y = Value, fill = SP), col = "black", alpha =0.8)+
  geom_errorbar(aes(x = SP, ymin = Value-sd, ymax = Value+sd))+
  scale_fill_manual(values = setNames(c("forestgreen", rep("tomato",4), rep("springgreen",3),"azure2"),as.character(levels(ctrl$SP)[-9])))+
  facet_wrap(~Day, nrow = 1)+theme_bw()+theme(
    axis.text.x = element_text(angle = 90, size = 15, hjust = 1, vjust = 0.5),
    strip.background = element_rect(fill = "white", colour =  "black"),
    strip.text = element_text(face = "bold", colour = "black", size = 15)
  )

p2 <- ggplot(nacl150)+
  geom_boxplot(aes(x = SP, y = Value, fill = SP), alpha = 0.8)+
  geom_text(aes(x = SP, y = 50, label = Sig), size = 12)+
  scale_fill_manual(values = setNames(c("forestgreen", rep("tomato",4), rep("springgreen",4),"azure2"),as.character(levels(ctrl$SP))))+
  facet_wrap(~Day, nrow = 1)+theme_bw()+theme(
    axis.text.x = element_text(angle = 90, size = 15, hjust = 1, vjust = 0.5),
    strip.background = element_rect(fill = "white", colour =  "black"),
    strip.text = element_text(face = "bold", colour = "black", size = 15)
  )
ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2, labels = "AUTO", font.label = list(size = 20))

###conclutions:
suba <- read.csv("suba5/19_8_2022_suba5_summary.csv") %>% select(c(1,3,40)) %>% mutate(locus = str_split_fixed(locus, "\\.",2)[,1]) %>% as_tibble() %>% distinct()

suba2 <- aggregate(suba, location_consensus~locus, FUN = function(x){paste(x, collapse = ";")})
suba2$description <- suba$description[match(suba2$locus, suba$locus)]

Expressionne <- openxlsx::read.xlsx("Transcriptome.merged.xlsx")

exp_true <- Expressionne %>% group_by(locus) %>% summarise(OEM15 = mean(M15),
                                                           OEP12 = mean(P12),
                                                           OEM20 = mean(M20),
                                                           OEP5 = mean(P5),
                                                           rpotmp = mean(Tmp),
                                                           `Tmp-M3`= mean(TM3),
                                                           `Tmp-M11`= mean(TM11),
                                                           `Tmp-P1`=mean(TP1)) %>% mutate(symbol = Expressionne$symbol[match(locus, Expressionne$locus)])



tf <- Expressionne %>% filter(locus %in% TF_list$Gene_ID) %>% mutate(tf_family = TF_list$Family[match(locus, TF_list$Gene_ID)])
colnames(tf)[10:17] <- trueCoefs

tf2 <- tf %>% group_by(locus) %>% summarise(OEM15 = mean(OEM15),
                                            OEP12 = mean(OEP12),
                                            OEM20 = mean(OEM20),
                                            OEP5 = mean(OEP5),
                                            rpotmp = mean(rpotmp),
                                            `Tmp-M3`= mean(`Tmp-M3`),
                                            `Tmp-M11`= mean(`Tmp-M11`),
                                            `Tmp-P1`=mean(`Tmp-P1`)) %>% mutate(symbol = tf$symbol[match(locus, tf$locus)],
                                                                                tf_family = tf$tf_family[match(locus, tf$locus)])

exp_true$localization <- suba2$location_consensus[match(exp_true$locus, suba2$locus)]

plastid <- exp_true[grep("plast", exp_true$localization),]

exp_true %>% filter(OEM15 > 0) %>% nrow()
exp_true %>% filter(OEM15 < 0) %>% nrow()

exp_true %>% filter(OEP12 > 0) %>% nrow()
exp_true %>% filter(OEP12 < 0) %>% nrow()


mne <- read.xlsx("mt-nuclear-encoded.xlsx")
pne <- read.xlsx("pt-nuclear-encoded.xlsx")

pne %>% filter(OEP12 < 0) %>% nrow()

mne %>% filter(OEM15 > 0) %>% nrow()

mse <- readRDS("mst_edges.RDS")
msv <- readRDS("mst_ver.RDS")
library(igraph)
mst <- graph_from_data_frame(mse, vertices =  msv, directed = F)
mst <- igraph::simplify(mst,edge.attr.comb = "first" )
V(mst)$symbol <- str_split_fixed(Expressionne$symbol[match(V(mst)$name, Expressionne$ProbeName)], ";",2)[,1]

library(ggraph)
ggraph(mst, "backbone")+
  geom_edge_link( alpha = 0.8, color = "grey")+ geom_node_point(aes(colour = Tmp), alpha = 0.8, size = 1)+
  theme_graph()+
  scale_color_gradientn(colors = c("firebrick1", "black", "chartreuse"))

sub81548 <- subgraph(mst, V(mst)$bet.com %in% c(8,15,48))
V(sub81548)$symbol[is.na(V(sub81548)$symbol)] <- V(sub81548)$TAIR[is.na(V(sub81548)$symbol)]


ggraph(sub81548, "backbone")+
  geom_edge_link( alpha = 0.8, color = "grey")+ geom_node_point(aes(colour = Tmp), alpha = 0.8, size = 3)+
  geom_node_text(aes(label = symbol), size = 1)+
  theme_graph()+
  scale_color_gradientn(colors = c("firebrick1", "dodgerblue", "black"))

sub4421 <- subgraph(mst, V(mst)$bet.com %in% c(44,21))
V(sub4421)$symbol[is.na(V(sub4421)$symbol)] <- V(sub4421)$TAIR[is.na(V(sub4421)$symbol)]

ggraph(sub4421, "kk")+
  geom_edge_link( alpha = 0.8, color = "grey")+ geom_node_point(aes(colour = Tmp), alpha = 0.8, size = 3)+
  geom_node_text(aes(label = symbol), size = 1)+
  theme_graph()+
  scale_color_gradientn(colors = c("firebrick1", "dodgerblue", "black", "chartreuse"))



sub91 <- subgraph(mst, V(mst)$bet.com %in% c(91))
V(sub91)$symbol[is.na(V(sub91)$symbol)] <- V(sub91)$TAIR[is.na(V(sub91)$symbol)]

ggraph(sub91, "kk")+
  geom_edge_link( alpha = 0.8, color = "grey")+ geom_node_point(aes(colour = Tmp), alpha = 0.8, size = 3)+
  geom_node_text(aes(label = symbol), size = 1)+
  theme_graph()+
  scale_color_gradientn(colors = c("firebrick1", "dodgerblue", "black"))
