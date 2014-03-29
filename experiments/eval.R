library(ggplot2)

d <- read.csv(file="trec-2005-and-profile-IDX_SAWIT.csv",sep=";")
d <- cbind(d,qry_and="RANKED-AND")
f <- read.csv(file="trec-2005-profile-IDX_SAWIT.csv",sep=";")
f <- cbind(f,qry_and="RANKED-OR")

i <- read.csv(file="trec-2006-and-profile-IDX_SAWIT.csv",sep=";")
i <- cbind(i,qry_and="RANKED-AND")
j <- read.csv(file="trec-2006-profile-IDX_SAWIT.csv",sep=";")
j <- cbind(j,qry_and="RANKED-OR")

d2 <- read.csv(file="trec-2005-and-profile-IDX_SAWIT2.csv",sep=";")
d2 <- cbind(d2,qry_and="RANKED-AND")
f2 <- read.csv(file="trec-2005-profile-IDX_SAWIT2.csv",sep=";")
f2 <- cbind(f2,qry_and="RANKED-OR")

i2 <- read.csv(file="trec-2006-and-profile-IDX_SAWIT2.csv",sep=";")
i2 <- cbind(i2,qry_and="RANKED-AND")
j2 <- read.csv(file="trec-2006-profile-IDX_SAWIT2.csv",sep=";")
j2 <- cbind(j2,qry_and="RANKED-OR")

g <- rbind(d,f)
g <- cbind(g,qryfile="trec2005")
h <- rbind(i,j)
h <- cbind(h,qryfile="trec2006")
l <- rbind(g,h)

g2 <- rbind(d2,f2)
g2 <- cbind(g2,qryfile="trec2005")
h2 <- rbind(i2,j2)
h2 <- cbind(h2,qryfile="trec2006")
l2 <- rbind(g2,h2)

q2 <- rbind(l,l2)


p <- ggplot(q2,aes(factor(k),qry_time/1000,fill=index)) 
p <- p + geom_boxplot()
p <- p + facet_grid(qryfile ~ qry_and)
p <- p + scale_y_log10(limits=c(0.1, 10000),breaks=c(1,10,100,1000,10000))
p <- p + annotation_logticks(sides = "lr") 
print(p)