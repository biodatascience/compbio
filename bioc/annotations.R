library("AnnotationHub")
ah <- AnnotationHub()
display(ah)
# search for species "Homo sapiens" and description "ChIP-seq"
ah["AH28856"]
peaks <- ah[["AH28856"]]
peaks
seqinfo(peaks)

# let's look up human-specific information
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
orgdb

columns(orgdb)
# one of my example genes: SFTPB
# https://gtexportal.org/home/gene/SFTPB
select(orgdb, "SFTPB", "GENENAME", "SYMBOL")
go <- select(orgdb, "SFTPB", "GO", "SYMBOL")
go <- go[go$ONTOLOGY == "BP",] # biological processes

library(GO.db)
columns(GO.db)
go2 <- select(GO.db, go$GO, c("TERM","DEFINITION"), "GOID")
split(go2, seq_len(nrow(go2)))

select(orgdb, "NR3C1", "GENENAME", "SYMBOL")
go <- select(orgdb, "NR3C1", "GO", "SYMBOL")
go <- go[go$ONTOLOGY == "BP",] # biological processes
go2 <- select(GO.db, go$GO, c("TERM","DEFINITION"), "GOID")
split(go2, seq_len(nrow(go2)))

