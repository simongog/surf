data <- read.csv2("trec8_wtdup_stat.txt",sep=",",header=F)

pdf("trec8_wtdup_stat.pdf")

plot(data$V1,cumsum(data$V2*data$V1)/crossprod(data$V2,data$V1),ylim=c(0,1),xlab="Node size",ylab="Ratio of covered elements in DUP")

#dev.copy2pdf(file="trec8_wtdup_stat.pdf")

dev.off()
