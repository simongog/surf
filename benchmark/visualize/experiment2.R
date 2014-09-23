library(RColorBrewer)
if ( !exists( "tikzDeviceLoaded" ) ){  
    require(tikzDevice) 
    options(tikzLatexPackages = 
              c(getOption('tikzLatexPackages'), 
                "\\input{/Users/sgog/Downloads/topk/macros.tex}"
               ) 
           )
	tikzDeviceLoaded = T
}
source("../../external/sdsl-lite/benchmark/basic_functions.R")

mypaste <- function(...){
    paste(...,sep="")
}

my_format2 <- function(...){
   format(..., nsmall=2,digits=2, big.mark=",")
}

tikz("experiment2.tex", width = 5.9, height = 1.6, standAlone = F)

par(las=1)
par(oma=c(0,0,0,0))
par(mar=c(2,4,1,0) + 0.1)
par(mgp=c(2,0.5,0))

raw <- data_frame_from_key_value_pairs("../results/experiment2.txt")
raw2 <- raw[,c("collection_file","index_name","time_per_query")]
raw3 <- raw2[order(raw2$collection_file),]
m <- matrix(raw3[["time_per_query"]], nrow=2)
#            nrow=length(unique(raw2[["collection_file"]])))
mycol <- brewer.pal(5, "PuBu")
comp_bv_col <- mycol[2] #"gray80"
bv_col <- mycol[4] # "gray20"

perm <- c(3,1,4,2,5)
m <- m[,perm]
names <- unique(gsub("../collections/","",raw3[["collection_file"]]))
names <- gsub("ENWIKI(.*)","\\\\ENWIKI\\1",names)
names <- gsub("gov2","\\\\GOVII",names)

barplot(m, ylim=c(0,140),
        beside=T,
        names.arg=names[perm],
        col=c(comp_bv_col,bv_col),
        ylab="Avg. time per query [$\\mu s$]",
        cex.names=1)
abline(h=seq(0,140,10), col="gray")
barplot(m,ylim=c(0,140),
        beside=T,
        names.arg=names[perm],
        col=c(comp_bv_col,bv_col),
        ylab="Avg. time per query [$\\mu s$]",
        cex.names=1,add=T)

legend("top", legend=c("rrr\\_vector<63>","bit\\_vector"),
#        title="$K^2$-treap bitvector", 
        fill=c(comp_bv_col, bv_col),
        bty="n")

#box("outer", col="blue") 
#box("figure", col="green")  
#box("plot", col="red")

dev.off()
