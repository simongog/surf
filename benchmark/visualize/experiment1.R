require(tikzDevice) 
library(RColorBrewer)
source("../../external/sdsl-lite/benchmark/basic_functions.R")

raw <- data_frame_from_key_value_pairs("../results/experiment1.txt")

mypaste <- function(...){
    paste(...,sep="")
}

my_format2 <- function(...){
   format(..., nsmall=2,digits=2, big.mark=",")
}

tikz("experiment1.tex", width = 5.9, height = 2.0, standAlone = F)

par(mfrow=c(1,2))
par(las=1)
par(yaxs="i")
par(oma=c(0,0,0,0))
par(mar=c(3,4,1,0.5) + 0.1)
par(mgp=c(2,0.5,0))


raw[["total_time"]] <- raw[["time_per_query"]]
raw2 <- split(raw, raw$multi_occ)

time_all <- raw2[["0"]][["total_time"]] 
time_ex_singltons <- raw2[["1"]][["total_time"]] 
match_time <- rep(raw2[["0"]][["total_time"]][[1]], length(time_all))
multi_time <- time_ex_singltons - match_time
singlton_time <- time_all - multi_time - match_time
ks <- raw2[["0"]][["k"]][-1]

data <- rbind(match_time, multi_time, singlton_time)[,-1]

mycol <- brewer.pal(6, "PuBu")
csa_match_col <-  mycol[6]#"gray10"
k2_treap_col <- mycol[4] #"gray30"
rmq_csa_col <- mycol[2] #"gray80"

partcol <- c(csa_match_col, k2_treap_col, rmq_csa_col)

barplot(data,  ylab="Avg. time per query [$\\mu s$]",
        xlab="k", names.arg= ks,  col=partcol, cex.axis=0.8, cex.names=0.8)
abline(h=seq(0,800,100), col="gray")
barplot(data, ylab="Avg. time per query [$\\mu s$]",
        xlab="k", names.arg= ks,  col=partcol, cex.axis=0.8, cex.names=0.8, add=T)



par(xpd=T)
legend(list(x = 0,y = 1300), # "topleft",
       legend=rev(c("CSA matching", "$K^2$-treap search","RMQC + CSA access")), fill=rev(partcol),
       bty="n", cex=0.9)
par(xpd=F)

all_occ    <-  (raw2[["0"]][["check_sum"]])[-1]
multi_occ  <-  (raw2[["1"]][["check_sum"]])[-1]
single_occ <-  all_occ - multi_occ 

queries1 <- raw2[["1"]][["queries"]][-1]
csa1_time <- rep(raw2[["1"]][["time_per_query"]][[1]],nrow(raw2[["1"]])-1)*queries1
multi_occ_time <- (raw2[["1"]][["time_per_query"]])[-1]*queries1
multi_occ_ex_csa_time <- multi_occ_time-csa1_time
avg_multi_occ_time <- multi_occ_time/multi_occ
avg_multi_occ_ex_csa_time <- multi_occ_ex_csa_time/multi_occ
avg_multi_occ_csa_time <- avg_multi_occ_time - avg_multi_occ_ex_csa_time


queries0 <- raw2[["0"]][["queries"]][-1]
csa0_time <- rep(raw2[["0"]][["time_per_query"]][[1]],nrow(raw2[["0"]])-1)*queries0
all_ex_csa_time <- (raw2[["0"]][["time_per_query"]])[-1]-csa0_time
all_time <- (raw2[["0"]][["time_per_query"]])[-1]*queries0

single_occ_time <- all_time - multi_occ_time
avg_single_occ_time <- single_occ_time/single_occ
single_occ_ex_csa_time <- all_ex_csa_time - multi_occ_ex_csa_time 
avg_single_occ_ex_csa_time <- single_occ_ex_csa_time/single_occ
avg_single_occ_csa_time <- avg_single_occ_time-avg_single_occ_ex_csa_time

avg_single_occ_ex_csa_time <- avg_single_occ_time - (csa0_time-csa1_time)/single_occ

time_per_doc <- raw2[["0"]][["time_per_query"]][-1]*queries0/all_occ

plot(NA, xlab="k", ylab="Avg. time per doc. [$\\mu s$]",
     ylim=c(1,1000), xlim=c(0.5, max(length(ks))), xaxt='n', log="y", bty="n",
     cex.axis=0.8)
abline(h=c(1,2,3,4,5,8,7,8,9,10,20,30,40,50,60,70,80,90,100), col="gray90", lwd=0.8)

axis(1, at=seq(1,length(ks)), labels=as.character(ks), tick=F, cex.axis=0.8)

mixed_col = "black"

par(xpd=T)
legend(list(x = 2,y = 3000),#  "topright",
        legend=c("$K^2$-treap retrieved","RMQC + CSA retrieved","weighted average"),
        fill=c(k2_treap_col, rmq_csa_col, mixed_col),
        bty="n", cex=0.9
        )

par(xpd=F)
for ( i in 1:nrow(raw2[["0"]][-1,]) ){

    lines( rep(i,2)-0.15,
           c(1, avg_multi_occ_time[i]),
           lwd=5, col=k2_treap_col, lend=1 ) 
    lines( rep(i,2)+0.15,
           c(1, time_per_doc[i]),
           lwd=5, col=mixed_col, lend=1 ) 
    lines( rep(i,2),
           c(1, avg_single_occ_time[i]),
           lwd=5, col=rmq_csa_col, lend=1 ) 
}

#box("outer", col="blue") 
#box("figure", col="green")  
#box("plot", col="red")

dev.off()
