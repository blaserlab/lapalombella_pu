res <- read_csv("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/Msgsea/K-Means 10 Features-2.csv")
res
res2 <- res %>%
  dplyr::select(FeatureName, LFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(FeatureName) %>%
  summarize(stat=mean(LFC))
res2



ranks <- deframe(res2)
head(ranks, 20)

library(ggplot2)

#C2
pathways.c2 <- gmtPathways("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/Msgsea/c2.all.v7.5.1.symbols.gmt")

pathways.c2 %>%
  head() %>%
  lapply(head)

fgseaRes <- fgsea(pathways = pathways.c2,
                  stats    = ranks,
                  eps = 0.0,
                  minSize  = 50,
                  maxSize  = 500,
                  nPermSimple = 1000)


head(fgseaRes[order(pval), ])

#GSEAtable
dev.off(); dev.new()
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.c2[topPathways], ranks, fgseaRes,
              gseaParam=0.5)
#waterfall

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj) %>%
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES))  +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="c2 pathways NES from GSEA") +
  theme_minimal()+
  theme(text = element_text(size = 5))

plotEnrichment(pathways.c2[["REACTOME_NEUTROPHIL_DEGRANULATION"]],
               ranks) + labs(title="REACTOME_NEUTROPHIL_DEGRANULATION")

plotEnrichment(pathways.c2[["MOOTHA_MITOCHONDRIA"]],
               ranks) + labs(title="MOOTHA_MITOCHONDRIA")

#Hallmark
pathways.hallmark <- gmtPathways("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/Msgsea/h.all.v7.5.1.symbols.gmt")

pathways.hallmark %>%
  head() %>%
  lapply(head)

fgseaRes <- fgsea(pathways = pathways.hallmark,
                  stats    = ranks,
                  eps = 0.0,
                  minSize  = 10,
                  maxSize  = 500,
                  nPermSimple = 1000)


head(fgseaRes[order(pval), ])

#GSEAtable
dev.off(); dev.new()
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.hallmark[topPathways], ranks, fgseaRes,
              gseaParam=0.5)
#waterfall

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj) %>%
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES))  +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal()+
  theme(text = element_text(size = 5))

plotEnrichment(pathways.hallmark[["HALLMARK_GLYCOLYSIS"]],
               ranks) + labs(title="HALLMARK_GLYCOLYSIS")

plotEnrichment(pathways.hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
               ranks) + labs(title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")

#c5
pathways.c5 <- gmtPathways("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/Msgsea/c5.all.v7.5.1.symbols.gmt")

pathways.c5 %>%
  head() %>%
  lapply(head)

fgseaRes <- fgsea(pathways = pathways.c5,
                  stats    = ranks,
                  eps = 0.0,
                  minSize  = 50,
                  maxSize  = 500,
                  nPermSimple = 1000)


head(fgseaRes[order(pval), ])

#GSEAtable
dev.off(); dev.new()
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.c5[topPathways], ranks, fgseaRes,
              gseaParam=0.5)
#waterfall

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj) %>%
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES))  +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="c5 pathways NES from GSEA") +
  theme_minimal()+
  theme(text = element_text(size = 5))

plotEnrichment(pathways.hallmark[["GOCC_SECRETORY_GRANULE"]],
               ranks) + labs(title="GOCC_SECRETORY_GRANULE")


#c7
pathways.c7 <- gmtPathways("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/Msgsea/c7.immunesigdb.v7.5.1.symbols.gmt")

pathways.c7 %>%
  head() %>%
  lapply(head)

fgseaRes <- fgsea(pathways = pathways.c7,
                  stats    = ranks,
                  eps = 0.0,
                  minSize  = 20,
                  maxSize  = 500,
                  nPermSimple = 1000)


head(fgseaRes[order(pval), ])

#GSEAtable
dev.off(); dev.new()
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.c7[topPathways], ranks, fgseaRes,
              gseaParam=0.5)
#waterfall

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj) %>%
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES))  +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="c5 pathways NES from GSEA") +
  theme_minimal()+
  theme(text = element_text(size = 5))

plotEnrichment(pathways.hallmark[["GOCC_SECRETORY_GRANULE"]],
               ranks) + labs(title="GOCC_SECRETORY_GRANULE")

#K-mean cluster 6

resc6 <- read_csv("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/Msgsea/K-Means 10 Features-cluster6.csv")
resc6
res2c6 <- resc6 %>%
  dplyr::select(FeatureName1, LFC) %>%
  na.omit() %>%
  distinct() %>%
  group_by(FeatureName1) %>%
  summarize(stat=mean(LFC))
res2c6

ranks <- deframe(res2c6)
head(ranks, 20)

#Hallmark
pathways.hallmark <- gmtPathways("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/Msgsea/h.all.v7.5.1.symbols.gmt")

pathways.hallmark %>%
  head() %>%
  lapply(head)

fgseaRes <- fgsea(pathways = pathways.hallmark,
                  stats    = ranks,
                  eps = 0.0,
                  minSize  = 10,
                  maxSize  = 500,
                  nPermSimple = 1000)


head(fgseaRes[order(pval), ])

#GSEAtable
dev.off(); dev.new()
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.hallmark[topPathways], ranks, fgseaRes,
              gseaParam=0.5)
#waterfall

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj) %>%
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES))  +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal()+
  theme(text = element_text(size = 5))

plotEnrichment(pathways.hallmark[["HALLMARK_MYC_TARGETS_V1"]],
               ranks) + labs(title="HALLMARK_MYC_TARGETS_V1")

plotEnrichment(pathways.hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               ranks) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")

#C2
pathways.c2 <- gmtPathways("~/network/X/Labs/Blaser/share/collaborators/lapalombella_pu_network/Msgsea/c2.all.v7.5.1.symbols.gmt")

pathways.c2 %>%
  head() %>%
  lapply(head)

fgseaRes <- fgsea(pathways = pathways.c2,
                  stats    = ranks,
                  eps = 0.0,
                  minSize  = 100,
                  maxSize  = 500,
                  nPermSimple = 1000)


head(fgseaRes[order(pval), ])

#GSEAtable
dev.off(); dev.new()
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.c2[topPathways], ranks, fgseaRes,
              gseaParam=0.5)

#waterfall

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj) %>%
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES))  +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="c2 pathways NES from GSEA") +
  theme_minimal()+
  theme(text = element_text(size = 5))
