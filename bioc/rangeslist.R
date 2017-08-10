library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

# using AnnotationHub:
# library(AnnotationHub)
# ah <- AnnotationHub()
# edb <- query(ah, c("EnsDb", "Hsapiens", "v87"))[[1]]

columns(edb)
txdf <- select(edb,
               keys=keys(edb, "GENEID"),
               columns=c("GENEID","TXID"),
               keytype="GENEID")
head(txdf,20)

table(table(txdf$GENEID))

ebt <- exonsBy(edb, by="tx")

gid <- "ENSG00000196839"
txs <- txdf$TXID[txdf$GENEID == gid]
txs

ebt2 <- ebt[txs]
ebt2

ebt2[1]
ebt2[[1]]
ebt2[[1]][1]
ebt2[[1]][length(ebt2[[1]])]

library(Gviz)

granges2df <- function(x) {
  df <- as(x, "data.frame")
  df <- df[,c("seqnames","start","end","strand","group_name")]
  colnames(df)[1] <- "chromosome"
  colnames(df)[5] <- "transcript"
  df
}

df <- granges2df(ebt2)
df$gene <- gid
grt <- GeneRegionTrack(df)

range(ebt2)
range(unlist(ebt2))

gax <- GenomeAxisTrack()
plotTracks(list(gax,grt), chromosome="20",
           transcriptAnnotation="transcript",
           from=44619522 - 10000, to=44652233 + 10000)

plotTracks(list(gax,grt), chromosome="20",
           transcriptAnnotation="gene",
           from=44619522 - 10000, to=44652233 + 10000)

plotTracks(list(gax,grt), chromosome="20",
           collapseTranscripts=TRUE,
           shape="arrow",
           transcriptAnnotation="gene",
           from=44619522 - 10000, to=44652233 + 10000)

df$feature <- "A"
df$feature[6] <- "B"
grt <- GeneRegionTrack(df)
plotTracks(list(gax,grt), chromosome="20",
           transcriptAnnotation="transcript",
           from=44619522 - 10000, to=44652233 + 10000,
           A="orange", B="dodgerblue")
