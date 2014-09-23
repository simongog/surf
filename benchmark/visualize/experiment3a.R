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

tikz("experiment3a.tex", width = 5.9, height = 2.5, standAlone = F)


raw <- data_frame_from_key_value_pairs("../results/experiment3.txt")
raw2 <- data_frame_from_key_value_pairs("../../external/sdsl-lite/benchmark/document_retrieval/results/all.txt")
raw4 <- data_frame_from_key_value_pairs("../results/experiment4.txt")
raw5 <- data_frame_from_key_value_pairs("../results/experiment5.txt")
raw4 <- cbind(raw4, TLE=rep(0, nrow(raw4)))
raw4 <- cbind(raw4, collection_file=rep("ENWIKISML", nrow(raw4)))
colnames(raw4) <- gsub("idx","index_name",colnames(raw4))

dcc13data <- raw4[c("collection_file","index_name","query_len","time_per_query","TLE")]
fastdata <- raw5[c("collection_file","index_name","query_len","time_per_query","TLE")]
fastdata <- subset(fastdata,fastdata$index_name=="IDX_NNXU")

names(raw2) <- gsub("TC_ID","collection_file",names(raw2))
names(raw2) <- gsub("IDX_ID","index_name",names(raw2))
raw3 <- data_frame_from_key_value_pairs("../../external/sdsl-lite/benchmark/document_retrieval/results/all_int.txt")
names(raw3) <- gsub("TC_ID","collection_file",names(raw3))
names(raw3) <- gsub("IDX_ID","index_name",names(raw3))

data <- rbind( 
    raw[c("collection_file","index_name","query_len","time_per_query","TLE")],
    raw2[c("collection_file","index_name","query_len","time_per_query","TLE")],
    raw3[c("collection_file","index_name","query_len","time_per_query","TLE")])
data <- subset(data,data[["TLE"]]==0)

data[["collection_file"]] <- gsub("../collections/","",data[["collection_file"]])
data[["index_name"]] <- as.character( data[["index_name"]] )

d <- split(data, data$collection_file)
collections <- names(d)

idx2pch <- list()
idx2pch[["IDX_NNX"]] <- 1
idx2pch[["IDX_NNX_INT"]] <- 1
idx2pch[["GREEDY"]] <- 2 
idx2pch[["GREEDYINT"]] <- 2 
idx2pch[["SORT"]] <- 3
idx2pch[["SORTINT"]] <- 3
idx2pch[["DCC13KN"]] <- 4 
idx2pch[["IDX_NNXU"]] <- 1

idx2name <- list()
idx2name[["IDX_NNX"]] <- "\\IDXNNX"
idx2name[["IDX_NNXU"]] <- "\\IDXNNXU"
idx2name[["IDX_NNX_INT"]] <- "\\IDXNNX"
idx2name[["GREEDY"]] <- "\\GREEDY"
idx2name[["GREEDYINT"]] <- "\\GREEDY"
idx2name[["SORT"]] <- "\\SORT"
idx2name[["SORTINT"]] <- "\\SORT"
idx2name[["DCC13KN"]] <- "\\IDXKN" 

idx2col <- list()
idx2col[["IDX_NNX"]] <- "red"
idx2col[["IDX_NNXU"]] <- "blue"
idx2col[["IDX_NNX_INT"]] <- "red"
idx2col[["GREEDY"]] <- "gray30"
idx2col[["GREEDYINT"]] <- "gray30"
idx2col[["SORT"]] <- "gray30"
idx2col[["SORTINT"]] <- "gray30"
idx2col[["DCC13KN"]] <- "blue" 


#par(mfrow=c(2,2))
par(las=1)
  # length of tick mark as a fraction of the height of a line of text, default=-0.5
par(tcl=-0.2) 
par(oma=c(3,4,0,0.2)) # outer margin (bottom,left,top,right)
par(mar=c(1,2,1.5,0.5)) # inner margin (bottom,left,top,right)



collection="gov2"
    cat(collection)
    cat("\n")
    d2 <- d[[collection]]
    if ( collection == "ENWIKISML" ) {
        d2 <- rbind(d2, dcc13data)
        d2 <- rbind(d2, fastdata)
    }
    d3 <- split(d2,d2[["index_name"]])
    if ( length(d3) >= 1 ) {

        colname <- gsub("ENWIKI(.*)","\\\\ENWIKI\\1",collection)
        colname <- gsub("gov2","\\\\GOVII",colname)

        xxx <- "s"
        yyy <- "s"

        plot(NA, xlim=c(1,20),ylim=c(10,1000), log="y",
             main=colname,
             xlab="Pattern length",
             ylab="",
             xaxt=xxx, yaxt="n"  )
        if ( xxx == "n" ){
            axis(1, labels=F)
        } else {
            mtext(1, text = "Pattern length $m$", line = 2.5, las=0)
        }
        if ( yyy == "n" ){
            axis(2, labels=F)
        } else {
            mtext(2, text = "Avg. time per query [$\\mu s$]", line = 4.5, las=0)
            axis(2, at=c(10,100,1000,10000,100000), labels=c("10","100","1,000","10,000","100,000"),las=1)
        }
        ff <- seq(1,9)
        abline(h=c(ff*10,ff*100,ff*1000,ff*10000,100000),col="gray")
        
        for ( idx in names(d3) ) {
            lines(d3[[idx]][["query_len"]], d3[[idx]][["time_per_query"]], 
                  pch=idx2pch[[idx]], type="b",col=idx2col[[idx]])
        }
        x <- c("IDX_NNX")

        legend("topright",
                legend=as.character(unlist(idx2name[x])),
                pch=as.vector(unlist(idx2pch[x])),
                col=as.vector(unlist(idx2col[x])),
                bty="n", cex=1
              )

    }
#box("outer", col="blue") 
#box("figure", col="green")  
#box("plot", col="red")

dev.off()
