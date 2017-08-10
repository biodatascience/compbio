library(GenomicRanges)

ir <- IRanges(start=100 + 1:10, width=10)
ir

ir <- IRanges(start=100 + 1:10, end=100 + 10:19)
ir

st <- rep(c("+","-"),5)
st

gr <- GRanges(seqnames="chr1", ranges=ir, strand=st)
gr

gr + 10

shift(gr, 10)

resize(gr, width=1)

flank(gr, width=2, both=TRUE)

seqnames(gr)
seqinfo(gr)

library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
g <- genes(edb)
length(g)

seqnames(g)
g <- keepStandardChromosomes(g, pruning.mode="coarse")
seqnames(g)
seqinfo(g)
genome(g)

table(seqnames(g))

