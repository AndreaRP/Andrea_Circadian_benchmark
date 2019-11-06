'%!in%' <- function(x,y)!('%in%'(x,y))

nls_test <- function(dataset, s2c){
  dataset$X <- rownames(dataset)
  datafr <- dataset %>% dplyr::select(ncol(dataset), everything())
  results <- nls_regression(data=datafr, s2c=s2c, mode="port")
  return(results)
}

jtk_test <- function(dataset, s2c_s){
  # dataset <- subsampled_exprs
  data <- dataset
  # s2c_s <- subsampled_s2c
  annot <- data.frame(Probeset=rownames(data))
  timepoints <- length(unique(s2c_s$time))
  nrep <- length(unique(s2c_s$ind))
  lag <- unique(diff(as.numeric(unique(s2c_s$time))))

  jtkdist(timepoints, nrep) # total time points, # replicates per time point
  
  periods <- ceiling(8/lag):round(30/lag) # number of time points per cycle. (10/6=; 20/6)
  # cat(timepoints, nrep, lag, periods)
  jtk.init(periods, lag)  # 4 is the number of hours between time points
  
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results.jtk <- cbind(annot,res,data)
  results.jtk <- results.jtk[order(res$ADJ.P,-res$AMP),]
  return(results.jtk)
}

bc_test <- function(dataset, s2c){
  data <- dataset
  timepoints <- length(unique(s2c$time))
  lag <- unique(diff(as.numeric(unique(s2c$time))))
  len <- timepoints*lag
  write.table(data
            , file=paste("./RAW/Hughes/temp_",len,"x", lag,".tsv",sep="")
            , quote=FALSE, sep='\t', row.names=T, col.names=NA)
  command <- paste("bash SW/BioCycle_run.sh -i '",paste("./RAW/Hughes/temp_",len,"x", lag,".tsv", sep="")
                   ,"' -o './RESULTS/Hughes/files/' -s 20 -e 30", sep="")
  system(command)
  file.remove(paste("./RAW/Hughes/temp_",len,"x", lag,".tsv",sep=""))
  file.rename(paste("./RESULTS/Hughes/files/BC_temp_", len, "x", lag, "_results.csv", sep="")
            , paste("./RESULTS/Hughes/files/GSE11923_BC_", len, "_", lag, ".csv", sep=""))
  bc_results <- read.csv(paste("./RESULTS/Hughes/files/GSE11923_BC_", len, "_", lag, ".csv", sep=""), sep="\t")
  return(bc_results)
}

ls_test <- function(dataset,s2c){
  write.csv(dataset, "./dataset_temp.csv")
  # write.csv(s2c, "./s2c_temp.csv")
  meta2d(infile="./dataset_temp.csv"
         , outdir = "./"
         , filestyle="csv"
         , timepoints = as.numeric(s2c$time)
         , minper = 8
         , maxper = 30
         , cycMethod = "LS"
         , analysisStrategy = "auto"
         , outputFile = TRUE
         , outIntegration = "both"
         , adjustPhase = "predictedPer"
         , combinePvalue = "fisher"
         , weightedPerPha = FALSE
         , ARSmle = "auto"
         , ARSdefaultPer = 24
         , outRawData = FALSE
         , releaseNote = TRUE
         , outSymbol = "")
  # remove temps
  file.remove("./dataset_temp.csv")
  # file.remove("./s2c_temp.csv")
  # Read results
  ls_result <- read.csv("./meta2d_dataset_temp.csv")
  file.remove("./LSresult_dataset_temp.csv")
  file.remove("./meta2d_dataset_temp.csv")
  return(ls_result)
}

ars_test <- function(dataset, s2c){
  write.csv(dataset, "./dataset_temp.csv")
  write.csv(s2c, "./s2c_temp.csv")
  meta3d(datafile="./dataset_temp.csv"
         , designfile="./s2c_temp.csv"
         , outdir = "./"
         , filestyle="csv"
         , design_libColm = 2 # sample ID
         , design_subjectColm = 4 # Subject ID
         , minper = 8, maxper = 30, cycMethodOne = "ARS"
         , timeUnit = "hour"
         , design_hrColm = 3 # order index of time point
         # , design_groupColm = NULL
         , design_libIDrename = NULL, adjustPhase = "predictedPer",
         combinePvalue = "fisher", weightedMethod = TRUE,
         outIntegration = "both", ARSmle = "auto", ARSdefaultPer = 24,
         dayZeroBased = FALSE, outSymbol = "")
  # remove temps
  file.remove("./dataset_temp.csv")
  file.remove("./s2c_temp.csv")
  # Read results
  ars_result <- read.csv("meta3dGroupID_AllSubjects_s2c_temp.csv")
  file.remove("./meta3dGroupID_AllSubjects_s2c_temp.csv")
  file.remove(list.files("./", pattern ="meta3dSubjectID_", full.names = TRUE))
  return(ars_result)
}

gen_samples <- function(total_genes, duration, interval, per_circ_genes, n_replicates, sdev=1, rseed=50){
  # total_genes=18800
  # duration=48
  # interval=1
  # per_circ_genes=0.2
  # n_replicates=3
  # sdev=1
  # rseed=50
  
  # Array to store the data
  temp <- data.frame(replicates = rep(c(1:n_replicates), duration/interval), timepoints = rep(seq(1,duration,by=interval), each=n_replicates))
  colnames_array <- apply(temp[,c("timepoints", "replicates")], 1, paste, collapse = "_rep")
  colnames_array <- paste("S", colnames_array, sep="")
  in_silico_data <- setNames(data.frame(matrix(ncol = length(colnames_array)+3, nrow = 0))
                             , c(colnames_array,"amp","per","phase"))
  
  # Non-circadian genes
  # set parameters settings
  n_points <- duration/interval*n_replicates
  fparams <- data.frame(m1 = n_points, m2 = 1, shape2 = 4, lb = 4, ub = 14, pde = 0.02, sym = 0.5)
  dparams <- data.frame(lambda1 = 0.13, lambda2 = 2, muminde = 1, sdde = 0.5)
  not_circ_n <- as.integer(total_genes*(1-per_circ_genes))
  # generate synthetic data without using real microarray data as seed
  mydata <- madsim(mdata=NULL, n=not_circ_n, ratio=0, fparams, dparams, sdev, rseed)
  data <- mydata$xdata[,1:n_points]
  # data <- paste("S",rep(seq(1,duration,by=interval),3), rep(c(1:3),each=duration/interval), sep="_")
  params <- data.frame("amp" = rep(NA,not_circ_n), "per" = rep(NA,not_circ_n), "phase" = rep(NA,not_circ_n))
  data <- cbind(data,params)
  # Assign data
  in_silico_data[1:not_circ_n,] <- data[1:not_circ_n,]
  
  # Circa genes
  # Get range of values
  min_val <- min(as.matrix(data[,1:(ncol(data)-3)]))
  max_val <- max(as.matrix(data[,1:(ncol(data)-3)]))
  mean_val <- mean(as.matrix(data[,1:(ncol(data)-3)]))
  # Create random values
  amps <- runif(total_genes*per_circ_genes, min=0.5, max=4)
  pers <- runif(total_genes*per_circ_genes, min=8, max=30)
  phases <- runif(total_genes*per_circ_genes, min=-6, max=6)
  # Define function
  circ = function(t,amp,per,phase){
    y = amp*sin(2*pi*t/per + phase)
    return(y)
  }
  
  circ_n <- as.integer(total_genes*per_circ_genes)
  for (r in 1:circ_n){
    # Seed (True) value
    res=as.numeric(unlist(lapply(list(seq(1,duration,by=interval)),circ,amps[r],pers[r],phases[r])))
    # Add a random value around the mean of the array data generated
    r_add <- rnorm(1, mean=mean_val, sd=sdev)
    res_1 <- res+r_add
    # Generate replicates with deviation from the true value
    for (z in 1:n_replicates){assign(paste("rep",z, sep=""), numeric())}
    # Generate the replicates around the seed
    for (i in 1:length(res)){
      for (rep in 1:n_replicates){
        assign(paste("rep",rep, sep=""), c(get(paste("rep",rep, sep="")),res_1[i]+rnorm(1,mean=res_1[i],sd=sdev)))
        # assign(paste("rep",rep, sep=""), c(get(paste("rep",rep, sep="")),res[i] + rnorm(1,mean=0,sd=sdev)))
      }
    }
    all_reps=numeric()
    for(z in 1:n_replicates){
      all_reps <- c(all_reps, get(paste("rep",z, sep="")))
    }
    long <- data.frame(data=all_reps #c(rep3,rep1,rep2)
                       ,ind=rep(c(1:n_replicates),each=duration/interval)
                       ,time=rep(seq(1,duration,by=interval),n_replicates))
    # unlist(unstack(long, data~ind+time))
    in_silico_data[nrow(in_silico_data)+1,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),amps[r],pers[r],phases[r])
    
  }
  
  # Shuffle rows
  in_silico_data <- in_silico_data[sample(nrow(in_silico_data)),]
  rownames(in_silico_data) <- c(1:nrow(in_silico_data))
  return(in_silico_data)
}

statistics_calc <- function(detected_genes, real_genes, real.n){
  # Check how many genes overlap between predicted and real
  real.p <- real_genes
  detected <- detected_genes
  
  tp <-  Reduce(intersect, list(real.p,detected))
  fp <- detected[detected %!in% real.p]
  fn <- real.p[real.p %!in% detected]
  tn <- real.n[real.n %!in% detected]
  
  df <- data.frame(rbind(tp=length(tp), fp=length(fp),fn=length(fn), tn=length(tn)))
  return(df)
}



