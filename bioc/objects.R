library(Biobase)
phenoData <- data.frame(sample=factor(1:6),
                        condition=factor(c("A","A","B","B","C","C")),
                        treated=factor(rep(0:1,3)))
phenoData

featureData <- data.frame(geneID=1:10, geneSymbol=letters[1:10])
featureData

exprs <- matrix(rnorm(6 * 10), ncol=6, nrow=10)

eset <- ExpressionSet(exprs,
                      AnnotatedDataFrame(phenoData),
                      AnnotatedDataFrame(featureData))
eset

pData(eset)
fData(eset)

idx <- eset$condition %in% c("A","B")
eset.sub <- eset[,idx]
pData(eset.sub)

levels(eset.sub$condition)

library(magrittr)
eset.sub$condition %<>% droplevels
levels(eset.sub$condition)

exprs(eset.sub)
exprs(eset)[,idx]

library(GEOquery)
epi <- getGEO("GSE4302")[[1]]

pData(epi)[1,]
epi$condition <- epi$characteristics_ch1
levels(epi$condition)
levels(epi$condition) %<>% (function(x) sub("sample type: ","",x))
levels(epi$condition)

exprs(epi)[1:5,1:5]
boxplot(exprs(epi),range=0)

library(matrixStats)
rv <- rowVars(exprs(epi))

library(pheatmap)
pheatmap(exprs(epi)[head(order(rv, decreasing=TRUE),100),],
         annotation_col=pData(epi)["condition"])

epi.sub <- epi[,epi$condition %in% c("Healthy control","Smoker")]
rv <- rowVars(exprs(epi.sub))
pheatmap(exprs(epi.sub)[head(order(rv, decreasing=TRUE),100),],
         annotation_col=pData(epi.sub)["condition"])


