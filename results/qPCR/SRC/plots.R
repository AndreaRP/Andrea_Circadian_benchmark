setwd("/data3/arubio/projects/Andrea_Circadian_benchmark/RESULTS/qPCR/")
library("dplyr")
library("ggplot2")
library("cowplot")
library("ArpyLib")
library("reshape2")
results <- read.csv("./qPCR_results_clean.csv", stringsAsFactors = F)
results$Sample.Name[!grepl("^SAMPLE1", results$Sample.Name)] <- gsub("SAMPLE", "SAMPLE0", results$Sample.Name[!grepl("^SAMPLE1", results$Sample.Name)])
results$Sample.Name[grep("^SAMPLE1$", results$Sample.Name)] <- "SAMPLE01"
# Assign ZT
results$ZT <- plyr::revalue(results$Sample.Name, c(SAMPLE01 = 10, SAMPLE02 = 11, SAMPLE03 = 13, SAMPLE04 = 15
                                                 , SAMPLE05 = 16, SAMPLE06 = 18, SAMPLE07 = 19, SAMPLE08 = 21
                                                 , SAMPLE09 = 22, SAMPLE10 = 24, SAMPLE11 = 1, SAMPLE12 = 3
                                                 , SAMPLE13 = 4, SAMPLE14 = 6, SAMPLE15 = 7, SAMPLE16 = 9)) 
results <- results[order(results$Gene, results$Sample.Name),]


# Avg
r <- as.data.frame(results %>% group_by(Gene, Sample.Name) %>% summarise(average = mean(Ct))) #(dplyr)
r$ZT <- plyr::revalue(r$Sample.Name, c(SAMPLE01 = 10, SAMPLE02 = 11, SAMPLE03 = 13, SAMPLE04 = 15
                                     , SAMPLE05 = 16, SAMPLE06 = 18, SAMPLE07 = 19, SAMPLE08 = 21
                                     , SAMPLE09 = 22, SAMPLE10 = 24, SAMPLE11 = 1, SAMPLE12 = 3
                                     , SAMPLE13 = 4, SAMPLE14 = 6, SAMPLE15 = 7, SAMPLE16 = 9)) 
r <- r[order(r$Gene, r$Sample.Name),]
ctrl <- r[which(r$Gene=="Hmbs"),]
ctrl$ind <- 1
norm.vector <- numeric()
for (gene in unique(results$Gene)){
  r1 <- r[which(r$Gene==gene),]
  norm.vector <- c(norm.vector, ctrl$average-r1$average)
}
# Norm
r$norm <- norm.vector
# Log
r$log <- 2^r$norm

# CircaN
s2c <- data.frame(sample=ctrl$Sample.Name
                  , time=as.numeric(ctrl$ZT)
                  , ind=ctrl$ind
                  , stringsAsFactors = F
)

circan.data <- dcast(data = r, formula = Gene~Sample.Name, fun.aggregate = sum, value.var = "norm")
circan.results <- circan(data=circan.data, s2c=s2c[,c("sample", "ind", "time")], mode="port")
rownames(circan.results) <- circan.results$feature

# Plot
candidates <- unique(results$Gene)
for (g in candidates){
  # Get the gene data
  gene.data <- r[which(r$Gene == g),]
  # Get the fit data
  g1 <- plyr::round_any(as.numeric(circan.results[g, c("estimate.amp", "estimate.per", "estimate.phase", "r")])
                      , accuracy = .01, f = floor)
  
  data.long <- data.frame(data = as.numeric(gene.data$norm)
                        , time = gene.data$ZT
                        , ind = 1
                        , condition = "WT"
                        , stringsAsFactors = F
                          )
  data.long <- data.long[order(as.numeric(data.long$time)),]
  p <- data.long
  p$time <- as.numeric(p$time)
  
  temp <- ggplot(p, aes(x = time, y = data, group = condition)) +
    stat_function(fun = function(t) g1[1]*sin(2*pi*t/g1[2] + g1[3]) + mean(p$data), size=0.2, colour="red") + # regression function
    geom_line(data=p, aes(group = interaction(condition, ind))
              , size=0.4, color="#1e8449") + # mean line
    geom_point(data=p, aes(group = interaction(condition, ind), color=condition)
               , size=0.7, colour="#2ecc71") +
    scale_x_continuous(name="Time", breaks=as.numeric(unique(results$ZT))) +
    ylab("Expression") +
    ggtitle(g) +
    theme_minimal()
  
  assign(paste("gene", g, sep="_"), plot_to_gtable(temp))
}

plot_list <- list()
for(i in paste("gene", candidates, sep="_")){
  plot_list[[i]] <- get(i)
}

plot_grid(plotlist=plot_list , ncol = 3)
