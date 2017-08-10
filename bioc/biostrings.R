library(Biostrings)
dna <- "AAACGCG"

dna <- DNAString("AAACGCG")
reverseComplement(dna)

dna <- DNAStringSet(c("AAA","CGC","TCT"))
dna

letterFrequency(dna, "C")
letterFrequency(dna, "C", as.prob=TRUE)
letterFrequency(dna, "CG", as.prob=TRUE)

dna <- DNAStringSet(c("AACTCAA","CTCTAAA","AAAGAG"))
vmatchPattern("CTC", dna)

dna <- DNAString("AAACTCAAAGAGAAATTTAAA")
pd <- PDict(c("CTC","GAG","TTT","AAA"))
matchPDict(pd, dna)

library(AnnotationHub)
ah <- AnnotationHub()
display(ah)

res <- query(ah, c("ensembl","GRCh38","dna.primary_assembly"))
res$sourceurl
genome <- query(ah, c("ensembl","GRCh38","dna.primary","release-87"))[[1]]

# Note! this downloads ~1 Gb to ~/.AnnotationHub

seqinfo(genome)
gr <- GRanges("1", IRanges(1e6 + 1, 1e6 + 100))
dna <- getSeq(genome, gr)
dna

txome <- query(ah, c("ensembl","GRCh38","cdna.all","release-87"))[[1]]
seqinfo(txome)
txs <- seqnames(seqinfo(txome))

edb <- query(ah, c("EnsDb", "Hsapiens", "v87"))[[1]]
ebg <- exonsBy(edb, by="tx")

ebg[1]
len <- sum(width(ebg[[1]]))

tx <- grep(names(ebg[1]), txs, value=TRUE)
tx.seq <- getSeq(txome, GRanges(tx, IRanges(1, len)))

exon.seq <- getSeq(genome, ebg[[1]])
unlist(exon.seq)

tx.seq == unlist(exon.seq)
