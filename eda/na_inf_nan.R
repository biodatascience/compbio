x <- c(1,2,3,NA)
mean(x)
sd(x)
mean(x, na.rm=TRUE)
sd(x, na.rm=TRUE)
is.na(x)
mean(x[!is.na(x)])

m <- matrix(c(0,0,0,1,1,1,1:9),ncol=3,byrow=TRUE)
rowMeans(m)
apply(m, 1, sd)
z <- apply(m, 1, function(x) mean(x)/sd(x))
z

z[!is.nan(z)]
z[!is.nan(z) & is.finite(z)]

apply(m, 1, function(x) if (sd(x) == 0) 0 else mean(x)/sd(x))
