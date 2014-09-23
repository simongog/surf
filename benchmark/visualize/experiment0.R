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

tikz("experiment0.tex", width = 5.9, height = 1.7, standAlone = F)

par(las=1)
par(oma=c(0,0,0,0))
par(mar=c(2,4,0.3,0) + 0.1)

raw <- read.csv2("../results/experiment0.txt", header = F)
raw <- raw[,c(3,5,6,8,9,10)]
colnames(raw) <- c("comp","comp_size","idx_size","col","col_size","idx")
raw[["frac_size"]] <- raw[["comp_size"]]/raw[["col_size"]]

data <- subset(raw, raw$idx == "IDX_NN")
tc <- c("ENWIKISML","ENWIKIBIG")

create_matrix <- function(data ,tc){
    m <- c()
    for ( c in  tc ) {
        d <- subset(data, data$col==c)[,c("frac_size","comp")]
        x <- c( subset(d,d$comp=="CSA")[[1]], 
                subset(d,d$comp=="H")[[1]]+subset(d,d$comp=="H_SELECT")[[1]],
                subset(d,d$comp=="W_AND_P")[[1]],
                subset(d,d$comp=="DOC")[[1]],
                subset(d, d$comp=="RMQ_C")[[1]],
                subset(d,d$comp=="BORDER")[[1]]+subset(d,d$comp=="BORDER_RANK")[[1]]
             )
        m <- cbind(m,x)
    }
    m
}

m <- create_matrix(data, tc)
tc_int <- c("ENWIKISMLINT","ENWIKIBIGINT","gov2")
m2 <- create_matrix(subset(raw, raw$idx == "IDX_NN_INT"), tc_int)
tcm <- c(tc,tc_int)
m <- cbind(m,m2)

comp_col <- brewer.pal(5, "PuBu")
barplot(m, ylim=c(0,3.0),names.arg=c("\\ENWIKISML","\\ENWIKIBIG","\\ENWIKISMLINT","\\ENWIKIBIGINT","\\GOVII"), cex.names=1,col=rev(comp_col),
        ylab="Space [fraction of input]", space=c(rep(0.2,length(tc)),0.5,rep(0.2,length(tc_int)-1))
        )
abline(h=seq(0.0, 3, 0.1), col="gray")

barplot(m, ylim=c(0,3.0), names.arg=c("\\ENWIKISML","\\ENWIKIBIG","\\ENWIKISMLINT","\\ENWIKIBIGINT","\\GOVII"), cex.names=1,col=rev(comp_col),
        ylab="Space [fraction of input]", space=c(rep(0.2,length(tc)),0.5,rep(0.2,length(tc_int)-1)), add=T
        )

legend("topright", legend=rev(c("CSA","H","$K^2$-treap","DOC","RMQC")),
        fill=comp_col,
        bty="n")
#box("outer", col="blue") 
#box("figure", col="green")  
#box("plot", col="red")

dev.off()
