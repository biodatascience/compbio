library(readr)
orgs <- read_tsv("organisms_genes.tsv")
head(orgs)
table(orgs$type)

library(magrittr)
orgs$type %<>% factor

library(rafalib)
bigpar()
with(orgs, stripplot(bp ~ type))
with(orgs, stripplot(log10(bp) ~ type))
orgs <- orgs[orgs$name != "Wheat",]
with(orgs, plot(bp, coding, col=type, pch=20, cex=2))
legend("topright", levels(orgs$type), col=1:nlevels(orgs$type), pch=20)
# with(orgs, identify(bp, coding, labels=orgs$name)) # right click to esc

orgs[order(orgs$bp, decreasing=TRUE)[1:5],c("name","bp")]
colnames(orgs)

library(dplyr)
orgs %>% mutate(genes = coding + noncoding) -> orgs
orgs$genes

orgs %>% group_by(type) %>% summarize(n=n(),
                                      cr.avg=mean(coding/genes),
                                      cr.sd=sd(coding/genes))

orgs %>% group_by(type) %>% summarize(n=n(),
                                      cr.avg=mean(coding/genes),
                                      cr.sd=sd(coding/genes)) -> tab

library(xtable)
xtable(tab)

library(knitr)
kable(tab)

library(ggplot2) # see: http://ggplot2.tidyverse.org/reference/
ggplot(orgs, aes(bp, coding)) + geom_point()
ggplot(orgs, aes(bp, coding, col=type)) + geom_point()

# here we have to name the median something other than `bp`
orgs %>%
  group_by(type) %>%
  summarize(basepairs=median(bp),
            min=min(bp),
            max=max(bp)) %>%
  ggplot(aes(type, basepairs, ymin=min, ymax=max)) +
  geom_pointrange()

ggplot(orgs, aes(bp, coding, col=type)) +
  geom_point() + geom_smooth(se=FALSE, method="lm")

ggplot(orgs, aes(bp, coding)) +
  geom_point() + facet_wrap(~ type, scales="free")

# overlapping densities and histograms
# a really useful plot, here with simulated data:
ns <- c(1000,500,200)
dat <- data.frame(type=rep(letters[1:3], ns),
                  z=rnorm(1700, mean=rep(c(0,2,4),ns), sd=rep(c(1,1,1),ns)))
ggplot(dat, aes(z, col=type, fill=type)) + geom_density(alpha=0.1)

ggplot(dat, aes(z, col=type, fill=type)) +
  geom_histogram(alpha=0.1, position="identity")

ggplot(dat, aes(z, col=type, fill=type)) +
  geom_histogram(position="dodge")
