setwd("E:/RPOT paper/")
setwd("D:/Postgrad/RPOT paper/")
library(tidyverse)
library(GGally)
library(openxlsx)
library(igraph)
#we can integrate more than two networks! One is expression correlation, 2 - PPI and 3 - TF targets!
tft <- read.delim("Ara_ALL.TFT.txt") %>% as_tibble()
tft <- tft %>% mutate(tf = str_split_fixed(Gene.set.name, "_",2)[,1],
                      from = tf, to = Gene)

# correlation

library(igraph)
library(GGally)
gc()
degs <- read.xlsx("DEGs_separate_with_coex_main.xlsx")
agiarel <- read.delim("agiarel.txt")
logex <- read.xlsx("coex2/log.expr.xlsx")
mods <- read.xlsx("probe2mod2tair.xlsx")

logex <- logex %>% filter(ProbeName %in% mods$genes)
g <- graph_from_adjacency_matrix(
  as.matrix(as.dist(cor(t(logex[,c(1:48)]), method="kendall"))),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)
g <- igraph::simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
V(g)$name <- logex$ProbeName
V(g)$ProbeName <- logex$ProbeName
V(g)$TAIR <- agiarel$locus[match(V(g)$ProbeName, agiarel$array_element_name)]
g <- delete_edges(g, E(g)[which(abs(E(g)$weight)<0.8)])
g <- subgraph(g, which(igraph::degree(g)>0))
E(g)[which(E(g)$weight<0)]$color <- "dodgerblue"
E(g)[which(E(g)$weight>0)]$color <- "firebrick1"
E(g)$weight0 <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
# g <- delete_edges(g, E(g)[which(E(g)$weight<0.8)])
# g <- subgraph(g, which(igraph::degree(g)>0))
range(E(g)$weight)
E(g)$weight <- 1-E(g)$weight
components <- decompose(g)
g <- components[[which.max(sapply(components, vcount))]]
mst <- mst(g, algorithm="prim")
mst.communities <- cluster_edge_betweenness(mst, weights=NULL, directed=FALSE)

#cluster_edge_betweenness() #edges are interpreted as distances, not as connection strengths.
#cluster_fast_greedy()#more weight = stronger connection
# cluster_label_prop() #more weight = stronger connection
# cluster_leiden() # A larger edge weight means a stronger connection for this function.


V(mst)$cluster <- mst.communities$membership + 1

write.table(as_data_frame(mst, "edges"), "kendall_mst_edges.tab")
write.table(as_data_frame(mst, "vertices"), "kendall_mst_vertices.tab")

V(mst)$encoding <- ifelse(grepl("ATMG", V(mst)$TAIR)==T, "Mt", "Nuc")
set.seed(42)
ggnet2(mst,  color = "cluster", edge.color = "color", edge.alpha = 0.6, node.alpha = 0.6,
       size = "degree",shape = "encoding", shape.palette = setNames(c(17,16), c("Mt", "Nuc")),
       palette = setNames(viridis::turbo(length(unique(V(mst)$cluster))), unique(V(mst)$cluster)),
       label = "bet.com", label.size = 1, mode = "kamadakawai",
       edge.size = 1-E(mst)$weight)

 # write.table(as_data_frame(mst, "edges"), "spearman_mst_edges.tab", sep = "\t")
 # write.table(as_data_frame(mst, "vertices"), "spearman_mst_vertices.tab", sep = "\t")

mdf <- tibble(ProbeName = V(mst)$ProbeName,
              TAIR = V(mst)$TAIR,
              Cluster.bet = V(mst)$bet.com)

degs$cluster <- mdf$Cluster.bet[match(degs$locus, mdf$TAIR)]
# degs$cluster.gre <- mdf$cluster.gre[match(degs$locus, mdf$TAIR)]

dme <- degs %>% group_by(locus) %>% summarise(M15 = mean(M15), # this is to omit repeats
                                              P12 = mean(P12),
                                              Tmp = mean(Tmp),
                                              TP1 = mean(TP1),
                                              TM3 = mean(TM3))
degs$symbol2 <- str_split_fixed(degs$symbol,pattern = ";",n = 2)[,1]

dme <- dme %>% mutate(ProbeName = mods$genes[match(dme$locus, mods$TAIR)],
  cluster = mdf$Cluster.bet[match(ProbeName, mdf$ProbeName)],
                      cluster.gre = mdf$cluster.gre[match(ProbeName, mdf$ProbeName)],
               mainhubs = degs$mainhubs[match(dme$locus, degs$locus)],
               main_ppi = degs$main_ppi[match(dme$locus, degs$locus)],
               module = degs$module[match(dme$locus, degs$locus)],
               symbol = degs$symbol2[match(dme$locus, degs$locus)])

dme2 <- dme %>% filter(is.na(cluster)==F)

dme2 <- dme2%>% rename(
  "OEM15" = "M15",
  "OEP12" = "P12",
  "rpotmp" = "Tmp",
  "Tmp-M3" = "TM3",
  "Tmp-P1" = "TP1"
)
col_fun5 = circlize::colorRamp2(c(7,5,2,0,-2,-3, - 7), c("orange", "yellow", "green", "black","blue", "slateblue", "red"))

library(ComplexHeatmap)
Heatmap(as.matrix(dme2[,c(2:6)]), col = col_fun5, name = "LFC", row_split = dme2$cluster,
        na_col = "grey9", row_labels = dme2$symbol,
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 15), 
        row_title_rot = 0, row_title_gp = gpar(fontsize = 7),
        right_annotation = HeatmapAnnotation(#`TF_family` = dme2$tf_family,
                                             `main_hubs` = dme2$mainhubs,
                                             `main_ppi` = dme2$main_ppi,
                                             `module` = dme2$module,
                                              `cluster` = dme2$cluster,
                                             which = "row",
                                             col =  list(#TF_family = structure(viridis::turbo(16),names = unique(dme2$tf_family)[-1]),
                                                         main_hubs = structure(viridis::viridis(5), names = unique(dme2$mainhubs)),
                                                         main_ppi = structure(viridis::viridis(5), names = unique(dme2$main_ppi)),
                                                         module = structure(c("springgreen","gold", "dodgerblue", "tomato", "mediumorchid"), names = unique(dme2$module))
                                                         ,
                                                         cluster = structure(viridis::turbo(87), names = unique(dme2$cluster))
                                             )
        )
)
term2gene <- tibble(TERM = dme2$cluster, GENE = dme2$locus)
genelists <- list()
for(i in c(1:5)){
  genelists[[i]] <- structure(names =dme2$locus, dme2[[colnames(dme2)[i+1]]])
  names(genelists)[i] <- colnames(dme2)[i+1]
  genelists[[i]] <- genelists[[i]][order(genelists[[i]], decreasing = T)]
}
gseres <- list()
for(i in c(1:5)){
  gseres[[i]] <- clusterProfiler::GSEA(genelists[[i]], pAdjustMethod = "none", TERM2GENE = term2gene, pvalueCutoff = 0.2) %>% as_tibble()}

for(i in c(1:5)){
  gseres[[i]]$Line = colnames(dme2)[i+1]
}

gseres <- do.call(rbind,gseres)
write.xlsx(gseres,"coex2/gsea_clusters2.xlsx")

dme3 <- dme2 %>% filter(cluster %in% gseres$ID)
dme3$cluster2 <- dme3$cluster
mst3 <- subgraph(mst, V(mst)$bet.com %in% unique(dme3$cluster))

ggnet2(mst3,  color = "bet.com", edge.color = "color", edge.alpha = 0.8, node.alpha = 0.6,
       size = "degree",
       palette = setNames(viridis::turbo(16), unique(V(mst3)$bet.com)),
       label = "bet.com", label.size = 1, mode= "kamadakawai", max_size = 5,
       edge.size = 1-E(mst3)$weight)

dme3$cluster2[dme3$cluster %in% c(36,43)] <- c("36--43")
dme3$cluster2[dme3$cluster %in% c(8,15,48)] <- c("8--15--48")
dme3$encoding <- ifelse(grepl("ATMG", dme3$locus)==T, "Mt", "Nuc")

Heatmap(as.matrix(dme3[,c(2:6)]), col = col_fun5, name = "LFC", row_split = dme3$cluster2,
        na_col = "grey9", row_labels = paste(dme3$locus, dme3$symbol),
        row_names_gp = gpar(fontsize = 2),
        column_names_gp = gpar(fontsize = 15), 
        row_title_rot = 0, row_title_gp = gpar(fontsize = 7),
        right_annotation = HeatmapAnnotation(#`TF_family` = dme3$tf_family,
                                             `main_hubs` = dme3$mainhubs,
                                             `main_ppi` = dme3$main_ppi,
                                             `module` = dme3$module,
                                             `cluster` = dme3$cluster,
                                             `encoding` = dme3$encoding,
                                             which = "row",
                                             col =  list(#TF_family = structure(viridis::turbo(6),names = unique(dme3$tf_family)[-1]),
                                                         main_hubs = structure(viridis::viridis(4), names = unique(dme3$mainhubs)),
                                                         main_ppi = structure(viridis::viridis(5), names = unique(dme3$main_ppi)),
                                                         module = structure(c("springgreen","gold", "tomato", "mediumorchid", "royalblue"), names = unique(dme3$module))
                                                          ,
                                                         cluster = structure(viridis::turbo(16), names = unique(dme3$cluster)),
                                                         encoding = structure(c("tomato", "royalblue"), names = c("Mt", "Nuc"))
                                             )
        )
)

V(mst3)$OEM15 <- degs$M15[match(V(mst3)$name, degs$ProbeName)]
V(mst3)$OEP12 <- degs$P12[match(V(mst3)$name, degs$ProbeName)]
V(mst3)$rpotmp <- degs$Tmp[match(V(mst3)$name, degs$ProbeName)]
V(mst3)$tmpp1 <- degs$TP1[match(V(mst3)$name, degs$ProbeName)]

V(mst3)$OEM15[is.na(V(mst3)$OEM15)==T] <- 0
V(mst3)$OEP12[is.na(V(mst3)$OEP12)==T] <- 0
V(mst3)$rpotmp[is.na(V(mst3)$rpotmp)==T] <- 0
V(mst3)$tmpp1[is.na(V(mst3)$tmpp1)==T] <- 0
V(mst3)$encoding <- ifelse(grepl("ATMG", V(mst3)$TAIR)==T, "Mt", "Nuc")

# write.table(as_data_frame(mst, "edges"), "mst_edges.tab", sep = "\t")
# write.table(as_data_frame(mst, "vertices"), "mst_vertices.tab", sep = "\t")
# write.table(degs, "degs_w_clusters.tab", sep = "\t")
# write.table(mdf, "mdf.tab", sep = "\t")

set.seed(666)
p1 <- ggnet2(mst3,  color = col_fun5(V(mst3)$OEM15), edge.color = "color", edge.alpha = 0.6, node.alpha = 0.6,
       size = "degree", shape = "encoding", shape.palette = setNames(c(17,16), c("Mt", "Nuc")),
       # palette = setNames(viridis::turbo(16), unique(V(mst3)$bet.com)),
       label = "bet.com", label.size = 1, mode= "kamadakawai", max_size = 5,
       edge.size = 1-E(mst3)$weight)
set.seed(666)
p2 <- ggnet2(mst3,  color = col_fun5(V(mst3)$OEP12), edge.color = "color", edge.alpha = 0.6, node.alpha = 0.6,
             size = "degree", shape = "encoding", shape.palette = setNames(c(17,16), c("Mt", "Nuc")),
             # palette = setNames(viridis::turbo(16), unique(V(mst3)$bet.com)),
             label = "bet.com", label.size = 1, mode= "kamadakawai", max_size = 5,
             edge.size = 1-E(mst3)$weight)
set.seed(666)
p3 <- ggnet2(mst3,  color = col_fun5(V(mst3)$rpotmp), edge.color = "color", edge.alpha = 0.6, node.alpha = 0.6,
             size = "degree", shape = "encoding", shape.palette = setNames(c(17,16), c("Mt", "Nuc")),
             # palette = setNames(viridis::turbo(16), unique(V(mst3)$bet.com)),
             label = "bet.com", label.size = 1, mode= "kamadakawai", max_size = 5,
             edge.size = 1-E(mst3)$weight)
set.seed(666)
p4 <- ggnet2(mst3,  color = col_fun5(V(mst3)$tmpp1), edge.color = "color", edge.alpha = 0.6, node.alpha = 0.6,
             size = "degree", shape = "encoding", shape.palette = setNames(c(17,16), c("Mt", "Nuc")),
             # palette = setNames(viridis::turbo(16), unique(V(mst3)$bet.com)),
             label = "bet.com", label.size = 1, mode= "kamadakawai", max_size = 5,
             edge.size = 1-E(mst3)$weight)

ggpubr::ggarrange(p1,p2,p3,p4, labels = c("OEM15", "OEP12", "rpotmp", "Tmp-P1"), nrow =2 , ncol =2, common.legend = T, legend = "right" )

mst4 <- subgraph(mst3, V(mst3)$bet.com %in% c(36,43))
V(mst4)$TAIR[V(mst4)$TAIR == "no_match"] <- V(mst4)$name[V(mst4)$TAIR == "no_match"]
set.seed(666)
p1 <- ggnet2(mst4,  color = col_fun5(V(mst4)$OEM15), edge.color = "color", edge.alpha = 0.6, node.alpha = 0.6,
             size = "degree", shape = "encoding", shape.palette = setNames(c(17,16), c("Mt", "Nuc")),
             # palette = setNames(viridis::turbo(16), unique(V(mst3)$bet.com)),
             label = "TAIR", label.size = 1, mode= "kamadakawai", max_size = 5,
             edge.size = 1-E(mst4)$weight - 0.5, label.trim = 12)
set.seed(666)
p2 <- ggnet2(mst4,  color = col_fun5(V(mst4)$rpotmp), edge.color = "color", edge.alpha = 0.6, node.alpha = 0.6,
             size = "degree", shape = "encoding", shape.palette = setNames(c(17,16), c("Mt", "Nuc")),
             # palette = setNames(viridis::turbo(16), unique(V(mst3)$bet.com)),
             label = "TAIR", label.size = 1, mode= "kamadakawai", max_size = 5,
             edge.size = 1-E(mst4)$weight - 0.5, label.trim = 12)
ggpubr::ggarrange(p1,p2, labels = c("OEM15", "rpotmp"), nrow =1 , ncol =2, common.legend = T, legend = "right" )
noquote(V(mst4)$name)
writeClipboard(V(mst4)$name)



gseres$psig <- -1*log10(gseres$p.adjust)
ggplot(gseres)+
  geom_point(aes(x = NES, y = ID, col = Line, alpha = psig), size = 4)+facet_wrap(~Line, scales = "free_x")+theme_bw()+xlab("Normalized enrichment score")+ ylab("Network cluster")

mst3 <- subgraph(mst2,V(mst2)$name[which(V(mst2)$color %in% gseres$ID)])
mst3 <- simplify(mst3)
E(mst3)$color <- "grey"

library(org.At.tair.db)
t2s <- AnnotationDbi::select(org.At.tair.db, keys = unique(V(mst3)$TAIR), keytype = "TAIR", columns = "SYMBOL") 
t2s <- t2s %>% distinct()
t2s <- aggregate(SYMBOL~TAIR, t2s, FUN = paste, collapse = ";")
t2s$SYMBOL <- str_split_fixed(t2s$SYMBOL,";",2)[,1]

V(mst3)$symbol <- t2s$SYMBOL[match(V(mst3)$TAIR, t2s$TAIR)]
V(mst3)$symbol[is.na(V(mst3)$symbol)==T] <- V(mst3)$TAIR[is.na(V(mst3)$symbol)==T]
V(mst3)$symbol[V(mst3)$symbol=="no_match"] <- V(mst3)$ProbeName[V(mst3)$symbol=="no_match"]
set.seed(42)
ggnet2(mst3,  node.color = "color", edge.color = "color", edge.alpha = 0.6, node.alpha = 0.6,
       size = "degree",
       palette = setNames(c("gold", "tomato", "springgreen2", "royalblue1", "darkmagenta", "cyan3", "chartreuse"), unique(V(mst3)$color)),
       label = F, label.size = 1, label.color = "label.color",
       edge.size = "weight", mode = "kamadakawai", label.trim = 12, max_size = 5)

mst_edges <-  as_data_frame(mst3)
mst_vertices <- tibble(name = V(mst3)$name,
                       TAIR = V(mst3)$TAIR,
                       ProbeName = V(mst3)$ProbeName,
                       symbol = V(mst3)$symbol,
                       cluster = V(mst3)$color)
# write.xlsx(mst_edges, "mst_edges.xlsx")
# write.xlsx(mst_vertices, "mst_vertices.xlsx")

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

tn <- select_top_n_nodes(mst3, 10)
tn <- tn %>% mutate(cluster = V(mst3)$color[match(tair,V(mst3)$name)],
                    TAIR = V(mst3)$TAIR[match(tair,V(mst3)$name)],
                    symbol = V(mst3)$symbol[match(tair,V(mst3)$name)],
                    ProbeName = V(mst3)$ProbeName[match(tair,V(mst3)$name)],
                    Description = degs$DESCRIPTION[match(TAIR, degs$locus)],
                    in_degs = ifelse(ProbeName %in% degs$ProbeName, "yes", "no"))

mitochondrial <- c(dme3$locus[grep("ATMG", dme3$locus)],dme3$locus[grep("AT2G07", dme3$locus)] )
dme3 <- dme3 %>% mutate(encoding = ifelse(dme3$locus %in% mitochondrial, "Mt", "Nuc"))
Heatmap(as.matrix(dme3[,c(2:6)]), col = col_fun5, name = "LFC", row_split = dme3$cluster2,
        na_col = "grey9", row_labels = paste(dme3$locus, dme3$symbol),
        row_names_gp = gpar(fontsize = 2),
        column_names_gp = gpar(fontsize = 15), 
        row_title_rot = 0, row_title_gp = gpar(fontsize = 7),
        right_annotation = HeatmapAnnotation(`TF_family` = dme3$tf_family,
                                             `main_hubs` = dme3$mainhubs,
                                             `main_ppi` = dme3$main_ppi,
                                             `module` = dme3$module,
                                             `cluster` = dme3$cluster,
                                             `encoding` = dme3$encoding,
                                             which = "row",
                                             col =  list(TF_family = structure(viridis::turbo(6),names = unique(dme3$tf_family)[-1]),
                                                         main_hubs = structure(viridis::viridis(3), names = unique(dme3$mainhubs)),
                                                         main_ppi = structure(viridis::viridis(3), names = unique(dme3$main_ppi)),
                                                         module = structure(c("springgreen","gold", "tomato", "mediumorchid"), names = unique(dme3$module)),
                                                         cluster = structure(viridis::turbo(7), names = unique(dme3$cluster)),
                                                         encoding = structure(c("firebrick1", "royalblue1"), names = c("Mt","Nuc")) ) ))

mitochondrial_mst <- c(V(mst3)$tair[grep("ATMG", V(mst3)$TAIR)],V(mst3)$TAIR[grep("AT2G07", V(mst3)$TAIR)] )
V(mst3)$encoding <- ifelse(V(mst3)$TAIR %in% mitochondrial_mst, "Mt", "Nuc")
V(mst3)$encoding[V(mst3)$TAIR %in% mitochondrial] <- "Mt-DEG"
V(mst3)$encoding[V(mst3)$TAIR %in% dme3$locus[dme3$encoding %in% "Nuc"]] <- "Nuc-DEG"

set.seed(42)
ggnet2(mst3,  node.color = "encoding", edge.color = "color", edge.alpha = 0.6, node.alpha = 0.6,
       size = "degree",
       palette = setNames(c("firebrick4", "firebrick1", "royalblue4", "dodgerblue"), c("Mt","Mt-DEG", "Nuc", "Nuc-DEG")),
       label = "symbol", label.size = 1, label.color = "label.color",
       edge.size = "weight", mode = "kamadakawai", label.trim = 12, max_size = 5)

oems <- dme3 %>% filter(OEM15 !=0 & encoding != "Mt" & cluster2 == "34--19--18--16" )
write.xlsx(oems, "oems.xlsx")
oeps <- dme3 %>% filter(OEP12 !=0 & encoding != "Mt" & cluster2 == "34--19--18--16" )
write.xlsx(oeps, "oeps.xlsx")

degs$Sequence[degs$locus == "AT2G07160"]

V(mst3)$OEM15 <- dme3$OEM15[match(V(mst3)$TAIR, dme3$locus)]
V(mst3)$OEM15[is.na(V(mst3)$OEM15)==T] <-0

V(mst3)$OEP12 <- dme3$OEP12[match(V(mst3)$TAIR, dme3$locus)]
V(mst3)$OEP12[is.na(V(mst3)$OEP12)==T] <-0
set.seed(42)
ggnet2(mst3,  node.color = col_fun5(V(mst3)$OEP12), edge.color = "color", edge.alpha = 0.6, node.alpha = 0.6,
       size = "degree",
       # palette = setNames(c("firebrick4", "firebrick1", "royalblue4", "dodgerblue"), c("Mt","Mt-DEG", "Nuc", "Nuc-DEG")),
       label = "symbol", label.size = 1, label.color = "label.color",
       edge.size = "weight", mode = "kamadakawai", label.trim = 12, max_size = 5)


basel_edges <- read.xlsx("basel interactome/Table S1.xlsx") %>% mutate(tair1 =str_split_fixed(AGI.KEY, " ",3)[,1],
                                                                       tair2 = str_split_fixed(AGI.KEY, " ",3)[,3])
basel_ver <- read.xlsx("basel interactome/Table S2.xlsx")
gc()

basel <- graph_from_data_frame(basel_edges[,c(7,8)], directed = F)
V(basel)$Desc <- basel_ver$Description[match(V(basel)$name, basel_ver$AGI)]
V(basel)$loc <- basel_ver$Localization.consensus[match(V(basel)$name, basel_ver$AGI)]
V(basel)$process <- basel_ver$GO.process[match(V(basel)$name, basel_ver$AGI)]
V(basel)$genename <- basel_ver$Name[match(V(basel)$name, basel_ver$AGI)]

ml <- c("M1","M2", "M3", "M4", "M5")
mlist <- list()
for(i in c(1:5)){
  mlist[[i]] <- mods %>% filter(modules == ml[i]) %>% dplyr::select(TAIR) %>% unique()
}

grl <- list()
for(i in c(1:5)){
  grl[[i]] <- simplify(graph_from_data_frame(basel_edges[basel_edges$tair1 %in% mlist[[i]]$TAIR|
                                                  basel_edges$tair2 %in% mlist[[i]]$TAIR,c(7,8)], directed = F))
  
}

m1 <- select_top_n_nodes(grl[[1]], 10) %>% mutate(module = "M1",
                                                  name = basel_ver$Name[match(tair, basel_ver$AGI)],
                                                  Desc = basel_ver$Description[match(tair, basel_ver$AGI)],
                                                  process = basel_ver$GO.process[match(tair, basel_ver$AGI)],
                                                  localiz = basel_ver$Localization.consensus[match(tair, basel_ver$AGI)])
m2 <- select_top_n_nodes(grl[[2]], 10)%>% mutate(module = "M2",
                                                 name = basel_ver$Name[match(tair, basel_ver$AGI)],
                                                 Desc = basel_ver$Description[match(tair, basel_ver$AGI)],
                                                 process = basel_ver$GO.process[match(tair, basel_ver$AGI)],
                                                 localiz = basel_ver$Localization.consensus[match(tair, basel_ver$AGI)])
m3 <- select_top_n_nodes(grl[[3]], 5)%>% mutate(module = "M3",
                                                name = basel_ver$Name[match(tair, basel_ver$AGI)],
                                                Desc = basel_ver$Description[match(tair, basel_ver$AGI)],
                                                process = basel_ver$GO.process[match(tair, basel_ver$AGI)],
                                                localiz = basel_ver$Localization.consensus[match(tair, basel_ver$AGI)])
m4 <- select_top_n_nodes(grl[[4]], 5)%>% mutate(module = "M4",
                                                name = basel_ver$Name[match(tair, basel_ver$AGI)],
                                                Desc = basel_ver$Description[match(tair, basel_ver$AGI)],
                                                process = basel_ver$GO.process[match(tair, basel_ver$AGI)],
                                                localiz = basel_ver$Localization.consensus[match(tair, basel_ver$AGI)])
m5 <- select_top_n_nodes(grl[[5]], 5)%>% mutate(module = "M5",
                                                name = basel_ver$Name[match(tair, basel_ver$AGI)],
                                                Desc = basel_ver$Description[match(tair, basel_ver$AGI)],
                                                process = basel_ver$GO.process[match(tair, basel_ver$AGI)],
                                                localiz = basel_ver$Localization.consensus[match(tair, basel_ver$AGI)])
unim <- rbind(m1,m2,m3,m4,m5)
# write.xlsx(unim, "top_nodes_basel.xlsx")

tfil <- tft %>% filter((from %in% dme$locus | from %in% dme$symbol) &
                 to %in% dme$locus)
tfg <- graph_from_data_frame(tfil[,c(5,6)])
V(tfg)$OEM15 <- dme$M15[match(V(tfg)$name, dme$locus)]
V(tfg)$OEP12 <- dme$P12[match(V(tfg)$name, dme$locus)]
V(tfg)$rpotmp <- dme$Tmp[match(V(tfg)$name, dme$locus)]
V(tfg)$symbol <- dme$symbol[match(V(tfg)$name, dme$locus)]
tfg <- simplify(tfg)
tfgm15 <- subgraph(tfg, V(tfg)[which(V(tfg)$OEM15 < 0)])
tfgm15 <- subgraph(tfgm15, V(tfgm15)[which(degree(tfgm15) > 0)])
ggnet2(tfgm15,  node.color = col_fun5(V(tfgm15)$OEM15), edge.color = "gold", edge.alpha = 0.6, node.alpha = 0.6,
       size = "degree",
       # palette = setNames(c("firebrick4", "firebrick1", "royalblue4", "dodgerblue"), c("Mt","Mt-DEG", "Nuc", "Nuc-DEG")),
       label = "symbol", label.size = 1, 
        label.trim = 12, max_size = 5)

oemgenes <- dme3 %>% filter(OEM15 !=0) %>% dplyr::select(locus)

oemgenes$locus[grep("AT2G07", oemgenes$locus)]

roots <- oemgenes$locus[grep("AT2G07", oemgenes$locus)]
endees <- setdiff(oemgenes$locus, roots)

test <- all_shortest_paths(basel, from = which(V(basel)$name %in% oemgenes),
                           to = which(V(basel)$name %in% endees),mode = "all"
)
# oh holy motherfucker... at last! it works!
tl <- list()
for(i in seq_along(test$res)){
  tl[[i]] <- make_empty_graph()+ test$res[[i]]$name +igraph::path(test$res[[i]]$name)  
}
ppr <- do.call(igraph::union, tl)
V(ppr)$type <- "path"
V(ppr)$type[which(V(ppr)$name %in% roots)] <- "root"
V(ppr)$type[which(V(ppr)$name %in% endees)] <- "ends"

V(ppr)$OEM15 <- dme$M15[match(V(ppr)$name, dme$locus)]
V(ppr)$OEM15[is.na(V(ppr)$OEM15)==T] <-0 
V(ppr)$symbol <- t2s$SYMBOL[match(V(ppr)$name, t2s$TAIR)]
# V(ppr)$symbol[is.na(V(ppr)$symbol)==T] <- V(ppr)$name[is.na(V(ppr)$symbol)==T]
V(ppr)$symbol[is.na(V(ppr)$symbol)==T] <- basel_ver$Name[match(V(ppr)$name[is.na(V(ppr)$symbol)==T], basel_ver$AGI)]
V(ppr)$symbol[is.na(V(ppr)$symbol)==T] <- V(ppr)$name[is.na(V(ppr)$symbol)==T]
V(ppr)$symbol[c(9,12,26,75)] <- V(ppr)$name[c(9,12,26,75)]
ggnet2(ppr,  node.color = col_fun5(V(ppr)$OEM15), edge.color = "grey", edge.alpha = 0.6, node.alpha = 0.6,
       node.shape = "type", shape.palette = c("root" = 15, "ends" = 17, "path"=16),
       size = "degree",
       # palette = setNames(c("firebrick4", "firebrick1", "royalblue4", "dodgerblue"), c("Mt","Mt-DEG", "Nuc", "Nuc-DEG")),
       label = "symbol", label.size = 1.5,
       label.trim = 20,
       max_size = 8)
m15in <- select_top_n_nodes(ppr,5) %>% mutate(name = basel_ver$Name[match(tair, basel_ver$AGI)],
                                              desc = basel_ver$Description[match(tair, basel_ver$AGI)],
                                              process = basel_ver$GO.process[match(tair, basel_ver$AGI)],
                                              localization = basel_ver$Localization.consensus[match(tair, basel_ver$AGI)])
tofind <- get.edgelist(ppr)%>% as_tibble()
write.xlsx(m15in, "interact_deg_oem15.xlsx")

######################

degs0 <- read.xlsx("DEGS_correct_extended.xlsx")

V(mst)$OEM15 <- degs0$M15[match(V(mst)$name, degs0$ProbeName)]
V(smst)$OEM15 <- degs0$M15[match(V(smst)$name, degs0$ProbeName)]
V(kmst)$OEM15 <- degs0$M15[match(V(kmst)$name, degs0$ProbeName)]

V(mst)$OEM15[is.na(V(mst)$OEM15)==T] <- 0
V(smst)$OEM15[is.na(V(smst)$OEM15)==T] <- 0
V(kmst)$OEM15[is.na(V(kmst)$OEM15)==T] <- 0

V(mst)$OEP12 <- degs0$P12[match(V(mst)$name, degs0$ProbeName)]
V(smst)$OEP12 <- degs0$P12[match(V(smst)$name, degs0$ProbeName)]
V(kmst)$OEP12 <- degs0$P12[match(V(kmst)$name, degs0$ProbeName)]

V(mst)$OEP12[is.na(V(mst)$OEP12)==T] <- 0
V(smst)$OEP12[is.na(V(smst)$OEP12)==T] <- 0
V(kmst)$OEP12[is.na(V(kmst)$OEP12)==T] <- 0

V(mst)$TP1 <- degs0$TP1[match(V(mst)$name, degs0$ProbeName)]
V(smst)$TP1 <- degs0$TP1[match(V(smst)$name, degs0$ProbeName)]
V(kmst)$TP1 <- degs0$TP1[match(V(kmst)$name, degs0$ProbeName)]

V(mst)$TP1[is.na(V(mst)$TP1)==T] <- 0
V(smst)$TP1[is.na(V(smst)$TP1)==T] <- 0
V(kmst)$TP1[is.na(V(kmst)$TP1)==T] <- 0

V(mst)$Tmp <- degs0$Tmp[match(V(mst)$name, degs0$ProbeName)]
V(smst)$Tmp <- degs0$Tmp[match(V(smst)$name, degs0$ProbeName)]
V(kmst)$Tmp <- degs0$Tmp[match(V(kmst)$name, degs0$ProbeName)]

V(mst)$Tmp[is.na(V(mst)$Tmp)==T] <- 0
V(smst)$Tmp[is.na(V(smst)$Tmp)==T] <- 0
V(kmst)$Tmp[is.na(V(kmst)$Tmp)==T] <- 0

V(mst)$TM3 <- degs0$TM3[match(V(mst)$name, degs0$ProbeName)]
V(smst)$TM3 <- degs0$TM3[match(V(smst)$name, degs0$ProbeName)]
V(kmst)$TM3 <- degs0$TM3[match(V(kmst)$name, degs0$ProbeName)]

V(mst)$TM3[is.na(V(mst)$TM3)==T] <- 0
V(smst)$TM3[is.na(V(smst)$TM3)==T] <- 0
V(kmst)$TM3[is.na(V(kmst)$TM3)==T] <- 0

saveRDS(as_data_frame(mst, "edges"), "coex2/Chistovik/mst_edges.RDS")
saveRDS(as_data_frame(smst, "edges"), "coex2/Chistovik/smst_edges.RDS")
saveRDS(as_data_frame(kmst, "edges"), "coex2/Chistovik/kmst_edges.RDS")

saveRDS(as_data_frame(mst, "vertices"), "coex2/Chistovik/mst_ver.RDS")
saveRDS(as_data_frame(smst, "vertices"), "coex2/Chistovik/smst_ver.RDS")
saveRDS(as_data_frame(kmst, "vertices"), "coex2/Chistovik/kmst_vers.RDS")

saveRDS(basel_edges, "coex2/Chistovik/basel_edges.RDS")
saveRDS(basel_ver, "coex2/Chistovik/basel_ver.RDS")

bla <- read.delim("blastres.tab")
bla1 <- read.delim("blast.unknown.probes.tab")
degs0$Sequence[degs0$ProbeName == "A_84_P756485"]
bs <- basel_edges[grep("ATMG", basel_edges$tair1),]
bs2 <- basel_edges[grep("ATMG", basel_edges$tair2),]
bs <- rbind(bs,bs2)
gppi <- graph_from_data_frame(bs[,c(7,8,6)], directed = F)
gppi <- simplify(gppi)
V(gppi)$OEM15 <- degs$M15[match(V(gppi)$name, degs$locus)]
V(gppi)$OEM15[is.na(V(gppi)$OEM15)==T] <- 0
V(gppi)$tmp <- degs$Tmp[match(V(gppi)$name, degs$locus)]
V(gppi)$tmp[is.na(V(gppi)$tmp)==T] <- 0
ggnet2(gppi, node.color = col_fun5(V(gppi)$tmp), edge.color = "grey", mode = "kamadakawai", 
       node.alpha = 0.8, edge.alpha = 0.6, node.size = "degree", label = T, max_size = 3)

BiocManager::install("AnnotationDbi")

saveRDS(tft, "coex2/Chistovik/tft.RDS")
saveRDS(unis, "coex2/Chistovik/unis.RDS")

tnodes <- read.xlsx("top_nodes_basel.xlsx")
degs$mainppi <- "No"
degs$mainppi[degs$locus %in% tnodes$tair[tnodes$module == "M1"]] <- "M1"
degs$mainppi[degs$locus %in% tnodes$tair[tnodes$module == "M2"]] <- "M2"
degs$mainppi[degs$locus %in% tnodes$tair[tnodes$module == "M3"]] <- "M3"
degs$mainppi[degs$locus %in% tnodes$tair[tnodes$module == "M4"]] <- "M4"
degs$mainppi[degs$locus %in% tnodes$tair[tnodes$module == "M5"]] <- "M5"
degs$symbol3
t2s$SYMBOL[match(degs$locus, t2s$TAIR)]
saveRDS(degs, "coex2/Chistovik/degs_separ.RDS")

blan <- read.delim("E:/RPOT paper/coex2/probes.atrtd.tab",header = F)
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
colnames(blan) <- unlist(str_split(readClipboard(), pattern = " "))
blan2 <- blan %>% filter(length == 60, pident == 100)
blan2 <- blan2 %>% mutate(sseqid = str_split_fixed(sseqid, "_",2)[,1],
                 sseqid = str_split_fixed(sseqid, "\\.",2)[,1])
blan3 <- blan2[,c(1,2)] %>% distinct()
blan4 <- aggregate(sseqid~qseqid, blan3, FUN = function(x){paste0(x,collapse = ";")})
blan5 <- data.frame(qseqid = agiarel$array_element_name,sseqid = agiarel$locus ) %>% filter(sseqid != "no_match")

blan4 <- rbind(blan4, blan5)

agidb <- read.csv("coex2/AGIDB2.csv")
agidb$locus <- agiarel$locus[match(agidb$ID, agiarel$array_element_name)]
agidb <- agidb[agidb$locus == "no_match",]
agidb <- agidb[is.na(agidb$SEQUENCE) == F,]
library(msa)
unks <- DNAStringSet(agidb$SEQUENCE)
names(unks) <- agidb$ID
writeXStringSet(unks, "coex2/all_unk.fa")

blan6 <- separate_rows(blan4, sseqid, sep = ";")
library(org.At.tair.db)
t2s <- AnnotationDbi::select(org.At.tair.db, keys = blan6$sseqid, keytype = "TAIR", columns = "SYMBOL") 
t2s <- t2s %>% distinct()
t2s <- aggregate(SYMBOL~TAIR, t2s, FUN = paste, collapse = ";")
# t2s$SYMBOL <- str_split_fixed(t2s$SYMBOL,";",2)[,1]

blan6 <- blan6 %>% mutate(symbols = NA,
                 symbols = t2s$SYMBOL[match(sseqid, t2s$TAIR)],
                 localization = basel_ver$Localization.consensus[match(sseqid, basel_ver$AGI)])
saveRDS(blan4, "coex2/Chistovik/probe2tair.aggr.rds")
topn <- read.xlsx("top_nodes_basel.xlsx")
blan6 <- blan6 %>% mutate(cluster = mdf$Cluster.bet[match(blan6$qseqid, mdf$ProbeName)],
       mainhubs = hubs$module.x[match(blan6$qseqid, hubs$ProbeName)],
       main_ppi = topn$module[match(sseqid, topn$tair)],
       module = mods$module[match(qseqid, mods$genes)]
       )
saveRDS(blan6, "coex2/Chistovik/probe2tair.long.rds")
blan6 <- mods
blan6$cluster1 <- mst_vertices$bet.com[match(blan6$qseqid, mst_vertices$ProbeName)]
blan6$usym <- unis$`Gene.Names.(primary)`[match(blan6$sseqid, unis$TAIR)] 

bla <- read.delim("blastres.tab")
V(mst3)$TAIR[is.na(V(mst3)$symbol)==T]
unks[names(unks) %in% V(mst3)$name[is.na(V(mst3)$symbol)==T]]

mitoz <- read.delim("coex2/probes.ChrM.tab")
colnames(mitoz) <- unlist(str_split(readClipboard(), pattern = " "))
mitoz <- mitoz %>% filter(pident ==100, length == 60)

blan4$mito <- NA
blan4$mito[blan4$qseqid %in% mitoz$qseqid] <- paste0("M:", mitoz$sstart[match(blan4$qseqid[blan4$qseqid %in% mitoz$qseqid], mitoz$qseqid)], "-probe")
agidb[agidb$ID %in% V(mst3)$ProbeName[is.na(V(mst3)$symbol)==T],]

saveRDS(mst3, "coex2/Chistovik/mst36-43.rds")
mst3 <- subgraph(mst2, V(mst2)$bet.com %in% c(36,43))
V(mst3)$symbol[V(mst3)$name == "A_84_P847099"] <- "rrn26"
V(mst3)$encoding[V(mst3)$name == "A_84_P756324"] <- "Nuc"
V(mst3)$symbol[V(mst3)$name == "A_84_P541948"] <- "orf215B"
V(mst3)$symbol[V(mst3)$name == "A_84_P852044"] <- "rrn26"
V(mst3)$symbol[V(mst3)$name == "A_84_P756819"] <- "nad2B"
V(mst3)$symbol[V(mst3)$name == "A_84_P803055"] <- "matR"
V(mst3)$symbol[V(mst3)$name == "A_84_P756806"] <- "orf215A"
V(mst3)$symbol[V(mst3)$name == "A_84_P710040"] <- "nad1C"
V(mst3)$symbol[V(mst3)$name == "A_84_P756465"] <- "matR"
V(mst3)$symbol[V(mst3)$name == "A_84_P789736"] <- "rps4"
V(mst3)$symbol[V(mst3)$name == "A_84_P756797"] <- "AT2G07773"
V(mst3)$encoding[V(mst3)$name == "A_84_P756797"] <- "Nuc"

V(mst3)$symbol[V(mst3)$name == "A_84_P800139"] <- "AT2G07706"
V(mst3)$encoding[V(mst3)$name == "A_84_P800139"] <- "Nuc"

V(mst3)$symbol[V(mst3)$name == "A_84_P756465"] <- "matR"
V(mst3)$symbol[V(mst3)$name == "A_84_P756465"] <- "matR"
V(mst3)$symbol[grepl("M:", V(mst3)$TAIR)] <- V(mst3)$TAIR[grepl("M:", V(mst3)$TAIR)]

V(mst3)$symbol[is.na(V(mst3)$symbol)==T] <- V(mst3)$TAIR[is.na(V(mst3)$symbol)==T]
V(mst3)$symbol <- gsub("AT2G07717", "nad4", V(mst3)$symbol)
V(mst3)$symbol[grepl("ATMG00220", V(mst3)$TAIR)] <- "coB"
V(mst3)$symbol[grepl("ATMG00920", V(mst3)$TAIR)] <- "orf215B"
V(mst3)$symbol[grepl("AT1G64795", V(mst3)$TAIR)] <- "RITA"
V(mst3)$symbol[grepl("A_84_P800386", V(mst3)$ProbeName)] <- "nad4"

V(mst3)$symbol[grepl("AT5G53048", V(mst3)$TAIR)] <- "MhpC"
V(mst3)$symbol[grepl("AT3G58270", V(mst3)$TAIR)] <- "PEARLI 4"
V(mst3)$symbol[grepl("ATMG00516", V(mst3)$TAIR)] <- "nad1C"
V(mst3)$symbol[grepl("ATMG00020", V(mst3)$TAIR)] <- "rrn26"
V(mst3)$symbol[grepl("A_84_P739859", V(mst3)$ProbeName)] <- "matR"

V(mst3)$usym <- mods$usym[match(V(mst3)$name, mods$qseqid)]
V(mst3)$symbol <- str_split_fixed(probe2sym$symbols[match(V(mst3)$name, probe2sym$qseqid)],";",2)[,1]
V(mst3)$mito <- blan4$mito[match(V(mst3)$name, blan4$qseqid)]
V(mst3)$symbol[is.na(V(mst3)$symbol)==T] <- V(mst3)$mito[is.na(V(mst3)$symbol)==T]
V(mst3)$symbol[is.na(V(mst3)$symbol)==T] <- V(mst3)$TAIR[is.na(V(mst3)$symbol)==T]
V(mst3)$symbol[is.na(V(mst3)$symbol)==T] <- V(mst3)$ProbeName[is.na(V(mst3)$symbol)==T]
V(mst3)$encoding[grepl("matR", V(mst3)$symbol)] <- "Mt"
V(mst3)$encoding[grepl("rps", V(mst3)$symbol)] <- "Mt"
V(mst3)$encoding[grepl("M:", V(mst3)$symbol)] <- "Mt"
V(mst3)$encoding[grepl("nad4", V(mst3)$symbol)] <- "Mt"

set.seed(7)

msdf <- igraph::as_data_frame(mst3, "vertices")
write.xlsx(msdf, "mst3df.xlsx")
