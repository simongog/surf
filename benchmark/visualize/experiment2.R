require(tikzDevice) 
source("../../external/sdsl-lite/benchmark/basic_functions.R")


mypaste <- function(...){
    paste(...,sep="")
}

my_format2 <- function(...){
   format(..., nsmall=2,digits=2, big.mark=",")
}

tikz("experiment2.tex", width = 6.3, height = 3, standAlone = F)

par(las=1)

raw <- data_frame_from_key_value_pairs("../results/experiment2.txt")
raw2 <- raw[,c("collection_file","index_name","time_per_query")]
raw3 <- raw2[order(raw2$collection_file),]
m <- matrix(raw3[["time_per_query"]], nrow=2)
#            nrow=length(unique(raw2[["collection_file"]])))
comp_bv_col <- "gray80"
bv_col <- "gray20"

perm <- c(3,1,4,2,5)
m <- m[,perm]
names <- unique(gsub("../collections/","",raw3[["collection_file"]]))

barplot(m,
        beside=T,
        names.arg=names[perm],
        col=c(comp_bv_col,bv_col),
        ylab="Avg time per query [$\\mu$secs]",cex.names=0.8)

legend("top", legend=c("rrr\\_vector<63>","bit\\_vector"),
        title="$K^2$-treap bitvector", 
        fill=c(comp_bv_col, bv_col),
        bty="n")

dev.off()


