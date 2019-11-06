library("nlme")
library("madsim")
library("EnvStats")
source("/data3/arubio/src/My_library.R")
setClass("harmony_data", slots=list(G1="data.frame", G2="data.frame", Comparison="data.frame", run.info="list"))


gen_samples <- function(total_genes, duration, interval, per_circ_genes, n_replicates, sdev=1, rseed=50){
  # total_genes=15000
  # duration=24
  # interval=4
  # per_circ_genes=0.2
  # n_replicates=3
  # sdev=1.5
  # rseed=50

  # total_genes = 15000
  # duration = 24
  # interval = 4
  # per_circ_genes = 0.2
  # n_replicates = 4
  # sdev = 1.5
  # rseed = 10

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
  amps <- c(runif(total_genes/2*per_circ_genes, min=-1.5, max=-0.5)
            , runif(total_genes/2*per_circ_genes, min=0.5, max=1.5))
  pers <- runif(total_genes*per_circ_genes, min=10, max=30)
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

gen_samples_bimodal <- function(total_genes, duration, interval, per_circ_genes, n_replicates, sdev=1){
  # We generate randomly ~14k amplitudes between -1,5 to -0.5 and 0.5 to 1.5 and n periods from 10 to 30
  # total_genes=15000
  # duration=24
  # interval=4
  # per_circ_genes=20
  # n_replicates=3
  # sdev=1.5
  
  amps <- c(runif(total_genes/2*per_circ_genes, min=-1.5, max=-0.5)
            , runif(total_genes/2*per_circ_genes, min=0.5, max=1.5))
  pers <- runif(total_genes*per_circ_genes, min=10, max=30)
  phases <- runif(total_genes*per_circ_genes, min=-6, max=6)
  
  circ = function(t,amp,per,phase){
    y = amp*sin(2*pi*t/per + phase)
    return(y)
  }
  
  # Array to store the data
  temp <- data.frame(replicates = rep(c(1:3), duration/interval), timepoints = rep(seq(1,duration,by=interval), each=3))
  colnames_array <- apply(temp[,c("timepoints", "replicates")], 1, paste, collapse = "_rep")
  colnames_array <- paste("S", colnames_array, sep="")
  in_silico_data <- setNames(data.frame(matrix(ncol = length(colnames_array)+3, nrow = 0))
                             , c(colnames_array,"amp","per","phase"))
  
  # Circa genes
  for (r in 1:as.integer(total_genes*per_circ_genes)){
    # Seed (True) value
    res=as.numeric(unlist(lapply(list(seq(1,duration,by=interval)),circ,amps[r],pers[r],phases[r])))
    # Generate replicates with deviation from the true value
    rep1=rep2=rep3=numeric()
    # Generate the replicates around the seed
    for (i in 1:length(res)){
      rep1=c(rep1,res[i]+rnormMix(1, mean1 = res[i], sd1 = sdev, mean2 = res[1]+6, sd2 = sdev*2, p.mix = 0.35))
      rep2=c(rep2,res[i]+rnormMix(1, mean1 = res[i], sd1 = sdev, mean2 = res[1]+6, sd2 = sdev*2, p.mix = 0.35))
      rep3=c(rep3,res[i]+rnormMix(1, mean1 = res[i], sd1 = sdev, mean2 = res[1]+6, sd2 = sdev*2, p.mix = 0.35))
    }
    long <- data.frame(data=c(rep3,rep1,rep2)
                       ,ind=rep(c(1:3),each=duration/interval)
                       ,time=rep(seq(1,duration,by=interval),3))
    # unlist(unstack(long, data~ind+time))
    in_silico_data[r,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),amps[r],pers[r],phases[r])
    
  }
  
  # Non-circa genes
  for (r in 1:as.integer(total_genes*(1-per_circ_genes))){
    # Seed (True); after standarization all the means are centered ~0 and all the sd are = 1
    res=rep(0, times=duration/interval)
    # Generate replicates with deviation from the true value
    rep1=rep2=rep3=numeric()
    # Generate the replicates around the seed
    for (i in 1:length(res)){
      rep1=c(rep1,res[i]+rnormMix(1, mean1 = res[i], sd1 = sdev, mean2 = res[1]+6, sd2 = sdev*2, p.mix = 0.35))
      rep2=c(rep2,res[i]+rnormMix(1, mean1 = res[i], sd1 = sdev, mean2 = res[1]+6, sd2 = sdev*2, p.mix = 0.35))
      rep3=c(rep3,res[i]+rnormMix(1, mean1 = res[i], sd1 = sdev, mean2 = res[1]+6, sd2 = sdev*2, p.mix = 0.35))
    }
    long <- data.frame(data=c(rep3,rep1,rep2)
                       ,ind=rep(c(1:3),each=duration/interval)
                       ,time=rep(seq(1,duration,by=interval),3))
    in_silico_data[nrow(in_silico_data)+1,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),NA,NA,NA)
  }
  
  # Shuffle rows
  in_silico_data <- in_silico_data[sample(nrow(in_silico_data)),]
  rownames(in_silico_data) <- c(1:nrow(in_silico_data))
  return(in_silico_data)
}

nls_regression <- function(data, s2c, shiny=FALSE, mode="default", init_value=24, max_per=Inf, min_per=-Inf){
  # Format s2c so nothing does weird things
  # s2c <- subsampled_s2c
  # data <- datafr
  # mode="port"
  
  s2c$time <- as.numeric(as.character(s2c$time))
  s2c$ind <- as.integer(as.character(s2c$ind))
  s2c$sample <- as.character(s2c$sample)
  # metadata must have sample, ind and time columns. Data is a data frame or matrix
  # with the ids in the first column
  t <- unique(as.numeric(as.character(s2c$time)))
  data0 <- data
  total_genes <- nrow(data)
  rownames(data) <- as.character(data[,1])
  data <- data[,-1]
  data <- as.matrix(data)
  results.cols <- c("feature"
                    ,"estimate.amp","std.error.amp","statistic.amp","p.value.amp"
                    ,"estimate.phase","std.error.phase","statistic.phase","p.value.phase"
                    ,"estimate.per" ,"std.error.per","statistic.per", "p.value.per"
                    , "r")
  results <- setNames(data.frame(matrix(ncol = 14, nrow = 0)), results.cols)
  for(gene in 1:nrow(data)){
    # Standarize data
    data.st <- (data[gene,]-mean(data[gene,]))/sqrt(var(data[gene,]))
    # Group data by timepoint
    df0 <- merge(data.st, s2c, by.x = 0, by.y = "sample")
    df0 <- dplyr::rename(df0, data = x)
    df <- as.data.frame(df0[,c("data", "ind", "time")])
    df$data <- as.numeric(df$data)
    df$time <- as.numeric(df$time)
    gd <- groupedData(data~time|ind,data=df)
    # plot(gd)
    result = tryCatch({
      # Fit
      nls.model=nls(  data~amp*sin(2*pi*time/per+phase)
                      , start=list(amp=1, phase=0, per=init_value)
                      , algorithm = mode
                      , lower = list(per=min_per)
                      , upper = list(per=max_per)
                      , data = gd
                      , control = list(maxiter = 200
                                       # , warnOnly = TRUE
                      )
      )
      # statistic values
      # stats <- tidy(nls.model)
      stats <- as.data.frame(tidy(nls.model))
      rownames(stats) <- stats[,1]
      stats <- stats[,-1]
      vec <- as.numeric(c(t(stats)))
      names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
      # R (goodness of fit)
      r <- cor(gd$data, predict(nls.model))
    }, error = function(e) {
      # print(e)
      result2 = tryCatch({
        # Try with fixed period
        nls.model=nls(  data~amp*sin(2*pi*time/init_value+phase)
                        , start=list(amp=1,phase=0)
                        , algorithm = mode
                        , data=gd
                        , control = list(maxiter = 200
                                         # , warnOnly = TRUE
                        )
        )
        # cat(stats)
        # statistic values
        # stats <- tidy(nls.model)
        stats <- as.data.frame(tidy(nls.model))
        rownames(stats) <- stats[,1]
        stats <- stats[,-1]
        vec <- as.numeric(c(t(stats)))
        names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
        # As we haven't estimated the period we add the data to vec
        temp <- c(24, "NA", "NA", "NA")
        names(temp) <- c("estimate.period", "std.error.period"
                         , "statistic.period", "p.value.period")
        vec <<- c(vec, temp)
        # R (goodness of fit)
        r <<- cor(gd$data, predict(nls.model))
      }, error = function(e) {
        # Can't be adjusted to a circadian curve
        vec <<- rep("NA", times = ncol(results)-1)
        r <<- "NA"
      })
    }, finally = {
      res <- rbind(c(feature=rownames(data)[gene], vec, r = r))
      results[gene,] <- res
    })
    if(shiny) incProgress(amount=1/total_genes)
  }
  # Calc adjusted p.val for period
  results$BH.q.value.per <- p.adjust(as.numeric(results$p.value.per), method="BH")
  results$BH.q.value.amp <- p.adjust(as.numeric(results$p.value.amp), method="BH")
  results$BH.q.value.phase <- p.adjust(as.numeric(results$p.value.phase), method="BH")
  return(results)
}

nls_regression_centr <- function(data, s2c, shiny=FALSE, mode="default", init_value=24, max_per=Inf, min_per=-Inf){
  # Format s2c so nothing does weird things
  s2c$time <- as.numeric(as.character(s2c$time))
  s2c$ind <- as.integer(as.character(s2c$ind))
  s2c$sample <- as.character(s2c$sample)
  # metadata must have sample, ind and time columns. Data is a data frame or matrix
  # with the ids in the first column
  t <- unique(as.numeric(as.character(s2c$time)))
  data0 <- data
  total_genes <- nrow(data)
  rownames(data) <- as.character(data[,1])
  data <- data[,-1]
  data <- as.matrix(data)
  results.cols <- c("feature"
                    ,"estimate.amp","std.error.amp","statistic.amp","p.value.amp"
                    ,"estimate.phase","std.error.phase","statistic.phase","p.value.phase"
                    ,"estimate.per" ,"std.error.per","statistic.per", "p.value.per"
                    , "r")
  results <- setNames(data.frame(matrix(ncol = 14, nrow = 0))
                      , results.cols)
  for(gene in 1:nrow(data)){
    # Standarize data
    data.st <- data[gene,]-mean(data[gene,])
    # Group data by timepoint
    df0 <- merge(data.st, s2c, by.x = 0, by.y = "sample")
    df0 <- dplyr::rename(df0, data = x)
    df <- as.data.frame(df0[,c("data", "ind", "time")])
    df$data <- as.numeric(df$data)
    df$time <- as.numeric(df$time)
    gd <- groupedData(data~time|ind,data=df)
    # plot(gd)
    result = tryCatch({
      # Fit
      nls.model=nls(data~amp*sin(2*pi*time/per+phase)
                  , start=list(amp=1, phase=0, per=init_value)
                  , algorithm = mode
                  , lower = list(per=min_per)
                  , upper = list(per=max_per)
                  , data = gd
                  , control = list(maxiter = 200
               # , warnOnly = TRUE
                  )
      )
      # statistic values
      # stats <- tidy(nls.model)
      stats <- as.data.frame(tidy(nls.model))
      rownames(stats) <- stats[,1]
      stats <- stats[,-1]
      vec <- as.numeric(c(t(stats)))
      names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
      # R (goodness of fit)
      r <- cor(gd$data, predict(nls.model))
    }, error = function(e) {
      # print(e)
      result2 = tryCatch({
        # Try with fixed period
        nls.model=nls(  data~amp*sin(2*pi*time/init_value+phase)
                        , start=list(amp=1,phase=0)
                        , algorithm = mode
                        , data=gd
                        , control = list(maxiter = 200
                                         # , warnOnly = TRUE
                        )
        )
        # cat(stats)
        # statistic values
        # stats <- tidy(nls.model)
        stats <- as.data.frame(tidy(nls.model))
        rownames(stats) <- stats[,1]
        stats <- stats[,-1]
        vec <- as.numeric(c(t(stats)))
        names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
        # As we haven't estimated the period we add the data to vec
        temp <- c(24, "NA", "NA", "NA")
        names(temp) <- c("estimate.period", "std.error.period"
                         , "statistic.period", "p.value.period")
        vec <<- c(vec, temp)
        # R (goodness of fit)
        r <<- cor(gd$data, predict(nls.model))
      }, error = function(e) {
        # Can't be adjusted to a circadian curve
        vec <<- rep("NA", times = ncol(results)-1)
        r <<- "NA"
      })
    }, finally = {
      res <- rbind(c(feature=rownames(data)[gene], vec, r = r))
      results[gene,] <- res
    })
    if(shiny) incProgress(amount=1/total_genes)
  }
  # Calc adjusted p.val for period
  results$BH.q.value.per <- p.adjust(as.numeric(results$p.value.per), method="BH")
  results$BH.q.value.amp <- p.adjust(as.numeric(results$p.value.amp), method="BH")
  results$BH.q.value.phase <- p.adjust(as.numeric(results$p.value.phase), method="BH")
  
  # if(shiny) incProgress(total_genes/total_genes)
  return(results)
}

gen_samples_paired <- function(total_genes=15000, duration=24, interval=4, per_circ_genes=0.2
                               , n_replicates=3, sdev=1, per_gained_genes=0.05, per_lost_genes=0.2){
  # total_genes = Total amount of genes to simulate
  # duration = Amount of hours of the experiment
  # interval = Interval between timepoints
  # per_circ_genes = % of TOTAL genes that will be rhythmic in G1
  # n_replicates = Number of replicates per timepoint
  # sdev = Standard deviation of the replicates
  # per_gained_genes = % of TOTAL genes that will gain rhythmicity in G2
  # per_lost_genes = % of RHYTHMIC genes that will lose rhythmicity in G2
  
  # We generate randomly ~14k amplitudes between -1,5 to -0.5 and 0.5 to 1.5 and n periods from 10 to 30
  amps <- c(runif(total_genes/2*per_circ_genes, min=-1.5, max=-0.5)
            , runif(total_genes/2*per_circ_genes, min=0.5, max=1.5))
  pers <- runif(total_genes*per_circ_genes, min=10, max=30)
  phases <- runif(total_genes*per_circ_genes, min=-6, max=6)
  
  circ = function(t,amp,per,phase){
    y = amp*sin(2*pi*t/per + phase)
    return(y)
  }
  
  # Array to store the data
  temp <- data.frame(replicates = rep(c(1:3), duration/interval), timepoints = rep(seq(1,duration,by=interval), each=3))
  colnames_array <- apply(temp[,c("timepoints", "replicates")], 1, paste, collapse = "_rep")
  colnames_array <- paste("S", colnames_array, "_G1", sep="")
  in_silico_data <- setNames(data.frame(matrix(ncol = length(colnames_array)+3, nrow = 0))
                             , c(colnames_array,"amp1","per1","phase1"))
  
  # Circa genes
  for (r in 1:as.integer(total_genes*per_circ_genes)){
    # Seed (True) value
    res=as.numeric(unlist(lapply(list(seq(1,duration,by=interval)),circ,amps[r],pers[r],phases[r])))
    # Generate replicates with deviation from the true value
    rep1=rep2=rep3=numeric()
    # after standarization all the means are centered ~0 and all the sd are = 1
    # Generate the replicates around the seed
    for (i in 1:length(res)){
      rep1=c(rep1,res[i]+rnorm(1,mean=res[i],sd=sdev))
      rep2=c(rep2,res[i]+rnorm(1,mean=res[i],sd=sdev))
      rep3=c(rep3,res[i]+rnorm(1,mean=res[i],sd=sdev))
    }
    long <- data.frame(data=c(rep3,rep1,rep2)
                       ,ind=rep(c(1:3),each=duration/interval)
                       ,time=rep(seq(1,duration,by=interval),3))
    # unlist(unstack(long, data~ind+time))
    in_silico_data[r,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),amps[r],pers[r],phases[r])
    
  }
  # Randomly select an 80% of the original genes that will share rhythimicity. 
  # The remaining 20% will lose rhythmicity
  shared <- sample(rownames(in_silico_data), size=total_genes*per_circ_genes*(1-per_lost_genes), replace=FALSE)
  lose <- rownames(in_silico_data)[-which(rownames(in_silico_data) %!in% shared)]
  # Non-circa genes
  for (r in 1:as.integer(total_genes*(1-per_circ_genes))){
    # Seed (True); after standarization all the means are centered ~0 and all the sd are = 1
    res=rep(0, times=duration/interval)
    # Generate replicates with deviation from the true value
    rep1=rep2=rep3=numeric()
    # Generate the replicates around the seed
    for (i in 1:length(res)){
      rep1=c(rep1,res[i]+rnorm(1,mean=res[i],sd=sdev))
      rep2=c(rep2,res[i]+rnorm(1,mean=res[i],sd=sdev))
      rep3=c(rep3,res[i]+rnorm(1,mean=res[i],sd=sdev))
    }
    long <- data.frame(data=c(rep3,rep1,rep2)
                       ,ind=rep(c(1:3),each=duration/interval)
                       ,time=rep(seq(1,duration,by=interval),3))
    in_silico_data[nrow(in_silico_data)+1,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),NA,NA,NA)
  }
  # Form the non rhythmic genes, a 5% will gain oscillation
  gain <- sample(rownames(in_silico_data[which(is.na(in_silico_data$amp1)),])
                 , size=total_genes*(1-per_circ_genes)*per_gained_genes, replace=FALSE)
    
  
  # We modify the original amps, pers and phases of the shared genes
  amp_changes <- sample(c(1, 2, 3), size = length(shared), replace = TRUE, prob = c(0.6, 0.3, 0.1))
  # amps2 <- amps*amp_changes
  per_changes <- sample(c(0, 2, 5, 10), size = length(shared), replace = TRUE, prob = c(0.6, 0.2, 0.1, 0.1))
  # pers2 <- pers+per_changes
  phase_changes <- sample(c(0, 2, 5, 10), size = length(shared), replace = TRUE, prob = c(0.6, 0.2, 0.1, 0.1))
  # phases2 <- phases+phase_changes
  # generate randomly the new gene stats
  new_amps <- c(runif(length(gain)/2, min=-1.5, max=-0.5)
            , runif(length(gain)/2, min=0.5, max=1.5))
  new_pers <- runif(length(gain), min=10, max=30)
  new_phases <- runif(length(gain), min=-6, max=6)
  
  # Create second group
  # Array to store the data
  temp <- data.frame(replicates = rep(c(1:3), duration/interval), timepoints = rep(seq(1,duration,by=interval), each=3))
  colnames_array <- apply(temp[,c("timepoints", "replicates")], 1, paste, collapse = "_rep")
  colnames_array <- paste("S", colnames_array, "_G2", sep="")
  in_silico_data2 <- setNames(data.frame(matrix(ncol = length(colnames_array)+3, nrow = 0))
                             , c(colnames_array,"amp2","per2","phase2"))
  
  shared_counter <- 1
  gain_counter <- 1
  for(i in 1:total_genes){
    if(i %in% shared){
      amp2 <- amps[as.numeric(i)]*amp_changes[shared_counter]
      per2 <- pers[as.numeric(i)]+per_changes[shared_counter]
      phase2 <- phases[as.numeric(i)]+phase_changes[shared_counter]
      shared_counter <- shared_counter+1
      res=as.numeric(unlist(lapply(list(seq(1,duration,by=interval)),circ,amp2,per2,phase2)))
      # Generate replicates with deviation from the true value
      rep1=rep2=rep3=numeric()
      # after standarization all the means are centered ~0 and all the sd are = 1
      # Generate the replicates around the seed
      for (j in 1:length(res)){
        rep1=c(rep1,res[j]+rnorm(1,mean=res[j],sd=sdev))
        rep2=c(rep2,res[j]+rnorm(1,mean=res[j],sd=sdev))
        rep3=c(rep3,res[j]+rnorm(1,mean=res[j],sd=sdev))
      }
      long <- data.frame(data=c(rep3,rep1,rep2)
                         ,ind=rep(c(1:3),each=duration/interval)
                         ,time=rep(seq(1,duration,by=interval),3))
      # unlist(unstack(long, data~ind+time))
      in_silico_data2[i,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),amp2,per2,phase2)
      # print("shared: ",i)
    }else{
      if(i %in% lose){
        # Generate non-oscilating gene
        # Seed (True); after standarization all the means are centered ~0 and all the sd are = sdev
        res=rep(0, times=duration/interval)
        # Generate replicates with deviation from the true value
        rep1=rep2=rep3=numeric()
        # Generate the replicates around the seed
        for (j in 1:length(res)){
          rep1=c(rep1,res[j]+rnorm(1,mean=res[j],sd=sdev))
          rep2=c(rep2,res[j]+rnorm(1,mean=res[j],sd=sdev))
          rep3=c(rep3,res[j]+rnorm(1,mean=res[j],sd=sdev))
        }
        long <- data.frame(data=c(rep3,rep1,rep2)
                           ,ind=rep(c(1:3),each=duration/interval)
                           ,time=rep(seq(1,duration,by=interval),3))
        in_silico_data2[i,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),NA,NA,NA)
      }else{
        if(i %in% gain){
          res=as.numeric(unlist(lapply(list(seq(1,duration,by=interval)),circ,new_amps[gain_counter]
                                       ,new_pers[gain_counter],new_phases[gain_counter])))
          
          gain_counter <- gain_counter+1
          # Generate replicates with deviation from the true value
          rep1=rep2=rep3=numeric()
          # after standarization all the means are centered ~0 and all the sd are = 1
          # Generate the replicates around the seed
          for (j in 1:length(res)){
            rep1=c(rep1,res[j]+rnorm(1,mean=res[j],sd=sdev))
            rep2=c(rep2,res[j]+rnorm(1,mean=res[j],sd=sdev))
            rep3=c(rep3,res[j]+rnorm(1,mean=res[j],sd=sdev))
          }
          long <- data.frame(data=c(rep3,rep1,rep2)
                             ,ind=rep(c(1:3),each=duration/interval)
                             ,time=rep(seq(1,duration,by=interval),3))
          # unlist(unstack(long, data~ind+time))
          in_silico_data2[i,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),new_amps[gain_counter]
                                   ,new_pers[gain_counter],new_phases[gain_counter])
        }else{
          # Just a non-rhythmic gene
          res=rep(0, times=duration/interval)
          # Generate replicates with deviation from the true value
          rep1=rep2=rep3=numeric()
          # Generate the replicates around the seed
          for (j in 1:length(res)){
            rep1=c(rep1,res[j]+rnorm(1,mean=res[j],sd=sdev))
            rep2=c(rep2,res[j]+rnorm(1,mean=res[j],sd=sdev))
            rep3=c(rep3,res[j]+rnorm(1,mean=res[j],sd=sdev))
          }
          long <- data.frame(data=c(rep3,rep1,rep2)
                             ,ind=rep(c(1:3),each=duration/interval)
                             ,time=rep(seq(1,duration,by=interval),3))
          in_silico_data2[i,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),NA,NA,NA)
        }
      }
    }
  }
   
  in_silico_data_whole <- cbind(in_silico_data, in_silico_data2)
  # Shuffle rows
  in_silico_data_whole <- in_silico_data_whole[sample(nrow(in_silico_data_whole)),]
  # Rename rows
  rownames(in_silico_data_whole) <- c(1:nrow(in_silico_data_whole))
  return(in_silico_data_whole)
  }
  
gen_samples_paired_simple <- function(total_genes=15000, duration=24, interval=4, per_circ_genes=0.2
                               , n_replicates=3, sdev=1, per_gained_genes=0.05, per_lost_genes=0.2){
  # total_genes = Total amount of genes to simulate
  # duration = Amount of hours of the experiment
  # interval = Interval between timepoints
  # per_circ_genes = % of TOTAL genes that will be rhythmic in G1
  # n_replicates = Number of replicates per timepoint
  # sdev = Standard deviation of the replicates
  # per_gained_genes = % of TOTAL genes that will gain rhythmicity in G2
  # per_lost_genes = % of RHYTHMIC genes that will lose rhythmicity in G2
  
  # total_genes = 15000
  # duration = 24
  # interval = 4
  # per_circ_genes = 0.2
  # n_replicates = 3
  # sdev = 1
  # per_gained_genes = 0.05
  # per_lost_genes = 0.2
  
  # We generate randomly ~14k amplitudes between -1,5 to -0.5 and 0.5 to 1.5 and n periods from 10 to 30
  amps <- c(runif(total_genes/2*per_circ_genes, min=-1.5, max=-0.5)
            , runif(total_genes/2*per_circ_genes, min=0.5, max=1.5))
  pers <- runif(total_genes*per_circ_genes, min=10, max=30)
  phases <- runif(total_genes*per_circ_genes, min=-6, max=6)
  
  circ = function(t,amp,per,phase){
    y = amp*sin(2*pi*t/per + phase)
    return(y)
  }
  
  # Array to store the data
  temp <- data.frame(replicates = rep(c(1:3), duration/interval), timepoints = rep(seq(1,duration,by=interval), each=3))
  colnames_array <- apply(temp[,c("timepoints", "replicates")], 1, paste, collapse = "_rep")
  colnames_array <- paste("S", colnames_array, "_G1", sep="")
  in_silico_data <- setNames(data.frame(matrix(ncol = length(colnames_array)+3, nrow = 0))
                             , c(colnames_array,"amp1","per1","phase1"))
  
  # Circa genes
  for (r in 1:as.integer(total_genes*per_circ_genes)){
    # Seed (True) value
    res=as.numeric(unlist(lapply(list(seq(1,duration,by=interval)),circ,amps[r],pers[r],phases[r])))
    # Generate replicates with deviation from the true value
    rep1=rep2=rep3=numeric()
    # after standarization all the means are centered ~0 and all the sd are = 1
    # Generate the replicates around the seed
    for (i in 1:length(res)){
      rep1=c(rep1,res[i]+rnorm(1,mean=res[i],sd=sdev))
      rep2=c(rep2,res[i]+rnorm(1,mean=res[i],sd=sdev))
      rep3=c(rep3,res[i]+rnorm(1,mean=res[i],sd=sdev))
    }
    long <- data.frame(data=c(rep3,rep1,rep2)
                       ,ind=rep(c(1:3),each=duration/interval)
                       ,time=rep(seq(1,duration,by=interval),3))
    # unlist(unstack(long, data~ind+time))
    in_silico_data[r,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),amps[r],pers[r],phases[r])
    
  }
  # Randomly select an 80% of the original genes that will share rhythimicity. 
  # The remaining 20% will lose rhythmicity
  shared <- sample(rownames(in_silico_data), size=total_genes*per_circ_genes*(1-per_lost_genes), replace=FALSE)
  lose <- rownames(in_silico_data)[-which(rownames(in_silico_data) %!in% shared)]
  # Non-circa genes
  for (r in 1:as.integer(total_genes*(1-per_circ_genes))){
    # Seed (True); after standarization all the means are centered ~0 and all the sd are = 1
    res=rep(0, times=duration/interval)
    # Generate replicates with deviation from the true value
    rep1=rep2=rep3=numeric()
    # Generate the replicates around the seed
    for (i in 1:length(res)){
      rep1=c(rep1,res[i]+rnorm(1,mean=res[i],sd=sdev))
      rep2=c(rep2,res[i]+rnorm(1,mean=res[i],sd=sdev))
      rep3=c(rep3,res[i]+rnorm(1,mean=res[i],sd=sdev))
    }
    long <- data.frame(data=c(rep3,rep1,rep2)
                       ,ind=rep(c(1:3),each=duration/interval)
                       ,time=rep(seq(1,duration,by=interval),3))
    in_silico_data[nrow(in_silico_data)+1,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),NA,NA,NA)
  }
  # Form the non rhythmic genes, a 5% will gain oscillation
  gain <- sample(rownames(in_silico_data[which(is.na(in_silico_data$amp1)),])
                 , size=total_genes*(1-per_circ_genes)*per_gained_genes, replace=FALSE)
  
  
  # We modify the original amps, pers and phases of the shared genes (but now only one thing will be changed)
  amp_changes <- sample(c(2, 3, 4), size = length(shared)/3, replace = TRUE, prob = c(0.6, 0.3, 0.1))
  amp_changes <- c(amp_changes, rep(1, (length(shared)*2/3)))

  per_changes <- sample(c(2, 5, 10), size = length(shared)/3, replace = TRUE, prob = c(0.6, 0.3, 0.1))
  per_changes <- c(rep(0, (length(shared)/3)), per_changes, rep(0, (length(shared)/3)))

  phase_changes <- sample(c(2, 5, 10), size = length(shared)/3, replace = TRUE, prob = c(0.6, 0.5, 0.1))
  phase_changes <- c( rep(0, (length(shared)*2/3)), phase_changes)

  # generate randomly the new gene stats
  new_amps <- c(runif(length(gain)/2, min=-1.5, max=-0.5)
                , runif(length(gain)/2, min=0.5, max=1.5))
  new_pers <- runif(length(gain), min=10, max=30)
  new_phases <- runif(length(gain), min=-6, max=6)
  
  # Create second group
  # Array to store the data
  temp <- data.frame(replicates = rep(c(1:3), duration/interval), timepoints = rep(seq(1,duration,by=interval), each=3))
  colnames_array <- apply(temp[,c("timepoints", "replicates")], 1, paste, collapse = "_rep")
  colnames_array <- paste("S", colnames_array, "_G2", sep="")
  in_silico_data2 <- setNames(data.frame(matrix(ncol = length(colnames_array)+3, nrow = 0))
                              , c(colnames_array,"amp2","per2","phase2"))
  
  shared_counter <- 1
  gain_counter <- 1
  for(i in 1:total_genes){
    if(i %in% shared){
      amp2 <- amps[as.numeric(i)]*amp_changes[shared_counter]
      per2 <- pers[as.numeric(i)]+per_changes[shared_counter]
      phase2 <- phases[as.numeric(i)]+phase_changes[shared_counter]
      shared_counter <- shared_counter+1
      res=as.numeric(unlist(lapply(list(seq(1,duration,by=interval)),circ,amp2,per2,phase2)))
      # Generate replicates with deviation from the true value
      rep1=rep2=rep3=numeric()
      # after standarization all the means are centered ~0 and all the sd are = 1
      # Generate the replicates around the seed
      for (j in 1:length(res)){
        rep1=c(rep1,res[j]+rnorm(1,mean=res[j],sd=sdev))
        rep2=c(rep2,res[j]+rnorm(1,mean=res[j],sd=sdev))
        rep3=c(rep3,res[j]+rnorm(1,mean=res[j],sd=sdev))
      }
      long <- data.frame(data=c(rep3,rep1,rep2)
                         ,ind=rep(c(1:3),each=duration/interval)
                         ,time=rep(seq(1,duration,by=interval),3))
      # unlist(unstack(long, data~ind+time))
      in_silico_data2[i,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),amp2,per2,phase2)
      # print("shared: ",i)
    }else{
      if(i %in% lose){
        # Generate non-oscilating gene
        # Seed (True); after standarization all the means are centered ~0 and all the sd are = sdev
        res=rep(0, times=duration/interval)
        # Generate replicates with deviation from the true value
        rep1=rep2=rep3=numeric()
        # Generate the replicates around the seed
        for (j in 1:length(res)){
          rep1=c(rep1,res[j]+rnorm(1,mean=res[j],sd=sdev))
          rep2=c(rep2,res[j]+rnorm(1,mean=res[j],sd=sdev))
          rep3=c(rep3,res[j]+rnorm(1,mean=res[j],sd=sdev))
        }
        long <- data.frame(data=c(rep3,rep1,rep2)
                           ,ind=rep(c(1:3),each=duration/interval)
                           ,time=rep(seq(1,duration,by=interval),3))
        in_silico_data2[i,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),NA,NA,NA)
      }else{
        if(i %in% gain){
          res=as.numeric(unlist(lapply(list(seq(1,duration,by=interval)),circ,new_amps[gain_counter]
                                       ,new_pers[gain_counter],new_phases[gain_counter])))
          
          gain_counter <- gain_counter+1
          # Generate replicates with deviation from the true value
          rep1=rep2=rep3=numeric()
          # after standarization all the means are centered ~0 and all the sd are = 1
          # Generate the replicates around the seed
          for (j in 1:length(res)){
            rep1=c(rep1,res[j]+rnorm(1,mean=res[j],sd=sdev))
            rep2=c(rep2,res[j]+rnorm(1,mean=res[j],sd=sdev))
            rep3=c(rep3,res[j]+rnorm(1,mean=res[j],sd=sdev))
          }
          long <- data.frame(data=c(rep3,rep1,rep2)
                             ,ind=rep(c(1:3),each=duration/interval)
                             ,time=rep(seq(1,duration,by=interval),3))
          # unlist(unstack(long, data~ind+time))
          in_silico_data2[i,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),new_amps[gain_counter]
                                   ,new_pers[gain_counter],new_phases[gain_counter])
        }else{
          # Just a non-rhythmic gene
          res=rep(0, times=duration/interval)
          # Generate replicates with deviation from the true value
          rep1=rep2=rep3=numeric()
          # Generate the replicates around the seed
          for (j in 1:length(res)){
            rep1=c(rep1,res[j]+rnorm(1,mean=res[j],sd=sdev))
            rep2=c(rep2,res[j]+rnorm(1,mean=res[j],sd=sdev))
            rep3=c(rep3,res[j]+rnorm(1,mean=res[j],sd=sdev))
          }
          long <- data.frame(data=c(rep3,rep1,rep2)
                             ,ind=rep(c(1:3),each=duration/interval)
                             ,time=rep(seq(1,duration,by=interval),3))
          in_silico_data2[i,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),NA,NA,NA)
        }
      }
    }
  }
  
  in_silico_data_whole <- cbind(in_silico_data, in_silico_data2)
  # Shuffle rows
  in_silico_data_whole <- in_silico_data_whole[sample(nrow(in_silico_data_whole)),]
  # Rename rows
  rownames(in_silico_data_whole) <- c(1:nrow(in_silico_data_whole))
  return(in_silico_data_whole)
}


nls_regression_groups <- function(data, s2c, id_col, sample_col, time_col
                                , group_col, ind_col, shiny=FALSE, mode="default"
                                , init_value=24, max_per=Inf, min_per=-Inf){
  
  # data=datafr
  # s2c=s2c.is
  # id_col=1
  # sample_col=1
  # time_col=3
  # group_col=4
  # ind_col=2
  # shiny=FALSE
  # mode="port"
  # init_value=24
  # max_per=Inf
  # min_per=-Inf
  
  # Format s2c so nothing does weird things
  s2c[,time_col] <- as.numeric(as.character(s2c[,time_col]))
  s2c[,ind_col] <- as.integer(as.character(s2c[,ind_col]))
  s2c[,sample_col]<- as.character(s2c[,sample_col])
  s2c[,group_col] <- as.character(s2c[,group_col])
  # Formatting data
  rown <- data[,id_col]
  data <- data[,-id_col]
  data[,1:ncol(data)] <- sapply(data[,1:ncol(data)], as.numeric)
  rownames(data) <- rown
  
  G1_name <- unique(s2c[,group_col])[1]
  G2_name <- unique(s2c[,group_col])[2]
  data1 <- data[,s2c[which(s2c[,group_col] == G1_name),sample_col]]
  data2 <- data[,s2c[which(s2c[,group_col] == G2_name),sample_col]]
  
  # data1 <- data[,s2c[which(s2c[,group_col] == unique(s2c[,group_col])[1]),sample_col]]
  # data2 <- data[,s2c[which(s2c[,group_col] == unique(s2c[,group_col])[2]),sample_col]]
  
  results.cols <- c("estimate.amp","std.error.amp","statistic.amp","p.value.amp"
                    ,"estimate.phase","std.error.phase","statistic.phase","p.value.phase"
                    ,"estimate.per" ,"std.error.per","statistic.per", "p.value.per"
                    , "r")
  group.cols <- c("feature", "estimate.Amp2_Amp1", "std.error.Amp2_Amp1", "statistic.Amp2_Amp1", "p.value.Amp2_Amp1"
                  , "estimate.Ph2_Ph1", "std.error.Ph2_Ph1", "statistic.Ph2_Ph1", "p.value.Ph2_Ph1"
                  , "estimate.Per2_Per1", "std.error.Per2_Per1", "statistic.Per2_Per1", "p.value.Per2_Per1")
  
  G1 <- setNames(data.frame(matrix(ncol = length(results.cols)+1, nrow = 0))
               , c("feature", paste(results.cols, "1", sep="")))
  G2 <- setNames(data.frame(matrix(ncol = length(results.cols)+1, nrow = 0))
               , c("feature", paste(results.cols, "2", sep="")))
  comp <- setNames(data.frame(matrix(ncol = length(group.cols), nrow = 0))
                   , group.cols)
  # results <- new("harmony_data", G1=G1, G2=G2, Comparison=comp)
  results <- new("harmony_data", G1=G1, G2=G2, Comparison=comp
                , run.info=list(G1=G1_name, G2=G2_name, id_col=id_col, sample_col=sample_col
                , time_col=time_col, group_col=group_col
                , ind_col=ind_col, mode=mode, init_value=init_value
                , max_per=max_per, min_per=min_per))
  for(gene in 1:nrow(data)){
    # Center data
    data.cent.1 <- (data1[gene,]-mean(as.numeric(data1[gene,])))
    # names(data.cent.1) <- names(data1[gene,])
    data.cent.2 <- (data2[gene,]-mean(as.numeric(data2[gene,])))
    # names(data.cent.2) <- names(data2[gene,])
    # Merged data
    merged_data <- cbind(data.cent.1, data.cent.2)
    rownames(merged_data) <- "data"
    # Construct df0 (ind must be different between groups)
    df0 <- merge(t(merged_data), s2c, by.x = 0, by.y = sample_col)
    # Change names to match formula
    df0 <- dplyr::rename(df0, group = colnames(s2c)[group_col]
                         , time = colnames(s2c)[time_col]
                         , ind = colnames(s2c)[ind_col])
    df0 <- df0[order(df0$group, df0$time),]
    df0$data <- as.numeric(df0$data)
    df0$time <- as.numeric(df0$time)
    # Group
    gd <- nlme::groupedData(data~time|ind, data=df0)
    
    
    # Possible try
    result = tryCatch({
      nls.model=stats::nls(data~
                          as.numeric(group==unique(s2c$group)[1])*amp1*sin(2*pi*time/per1+phase1) + 
                          as.numeric(group==unique(s2c$group)[2])*amp2*sin(2*pi*time/per2+phase2)
                        , algorithm = mode
                        , lower = list(per=min_per)
                        , upper = list(per=max_per)
                        , data = gd
                        , start = c(amp1=1, phase1=0, per1=init_value, amp2=1, phase2=0, per2=init_value)
      )
      
      # Each group' stats separately
      # summary(m.str)
      # stats <- broom::tidy(nls.model)
      stats <- as.data.frame(broom::tidy(nls.model))
      rownames(stats) <- stats[,1]
      coeff <- stats$term
      stats <- stats[,-1]
      # Populate G1 and G2
      vec <- as.numeric(c(t(stats)))
      names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
      vec1 <- vec[grep("1", names(vec))]
      vec2 <- vec[grep("2", names(vec))]
      # R (goodness of fit)
      r1 <- cor(gd[which(gd$group==unique(s2c$group)[1]),"data"]
                , predict(nls.model)[1:length(which(gd$group==unique(s2c$group)[1]))])
      r2 <- cor(gd[which(gd$group==unique(s2c$group)[2]),"data"]
                , predict(nls.model)[(length(which(gd$group==unique(s2c$group)[1]))+1):length(gd$group)])
      
      # Groups comparison
      K <- rbind(c(-1, 0, 0, 1 , 0, 0)
               , c( 0, -1, 0, 0, 1, 0)
               , c( 0, 0, -1, 0, 0, 1)
      )
      rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1", "Per2_Per1")
      colnames(K) <- names(coef(nls.model))
      
      group.stats <- as.data.frame(broom::tidy(summary(multcomp::glht(nls.model, linfct = K))))
      rownames(group.stats) <- group.stats[,1]
      coeff.group <- group.stats$lhs
      group.stats <- group.stats[,-c(1:2)]
      gs.vec <- as.numeric(c(t(group.stats)))
      names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
    }, error = function(err1) {
      convergence_error1 <- grepl("Convergence failure", as.character(err1))
      if (convergence_error1){
        result2 = tryCatch({ ## Per 1 fixed
          # Try with fixed period1
          nls.model=stats::nls(data~
                          as.numeric(group==unique(s2c$group)[1])*amp1*sin(2*pi*time/init_value+phase1) +
                          as.numeric(group==unique(s2c$group)[2])*amp2*sin(2*pi*time/per2+phase2)
                        , algorithm = mode
                        , lower = list(per=min_per)
                        , upper = list(per=max_per)
                        , data = gd
                        , start = c(amp1=1, phase1=0, amp2=1, phase2=0, per2=init_value)
          )
          # cat(stats)
          # statistic values
          stats <- broom::tidy(nls.model)
          rownames(stats) <- stats[,1]
          stats <- stats[,-1]
          vec <- as.numeric(c(t(stats)))
          names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
          # As we haven't estimated period 1 we add the data to vec1
          vec1 <- vec[grep("1", names(vec))]
          temp <- c(init_value, "NA", "NA", "NA")
          names(temp) <- c("estimate.per1", "std.error.per1"
                           , "statistic.per1", "p.value.per1")
          vec1 <- c(vec1, temp)
          mode(vec1) <- "numeric"
          vec2 <- vec[grep("2", names(vec))]
          # R (goodness of fit)
          r1 <- cor(gd[which(gd$group==unique(s2c$group)[1]),"data"]
                    , predict(nls.model)[1:length(which(gd$group==unique(s2c$group)[1]))])
          r2 <- cor(gd[which(gd$group==unique(s2c$group)[2]),"data"]
                    , predict(nls.model)[(length(which(gd$group==unique(s2c$group)[1]))+1):length(gd$group)])
          
          # Groups comparison
          K <- rbind(c(-1, 0, 1, 0, 0)
                     , c(0, -1, 0, 1, 0)
          )
          rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1")
          colnames(K) <- names(coef(nls.model))
          
          # group.stats <- tidy(summary(glht(nls.model, linfct = K)))
          group.stats <- as.data.frame(broom::tidy(summary(multcomp::glht(nls.model, linfct = K))))
          rownames(group.stats) <- group.stats[,1]
          coeff.group <- group.stats$lhs
          group.stats <- group.stats[,-c(1:2)]
          gs.vec <- as.numeric(c(t(group.stats)))
          names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
          gs.vec <- c(gs.vec, setNames(as.numeric(c(vec2["estimate.per2"]-vec1["estimate.per1"], NA, NA, NA))
                                       , c("estimate.Per2_Per1", "std.error.Per2_Per1"
                                           , "statistic.Per2_Per1", "p.value.Per2_Per1")))
          
        }, error = function(err2) {
          convergence_error2 <- grepl("Convergence failure", as.character(err2))
          if(convergence_error2){
            result = tryCatch({
              # Fixed per2
              nls.model=stats::nls(data~
                              as.numeric(group==unique(s2c$group)[1])*amp1*sin(2*pi*time/per1+phase1) +
                              as.numeric(group==unique(s2c$group)[2])*amp2*sin(2*pi*time/init_value+phase2)
                            , algorithm = mode
                            , lower = list(per=min_per)
                            , upper = list(per=max_per)
                            , data = gd
                            , start = c(amp1=1, phase1=0, amp2=1, phase2=0, per1=init_value)
              )
              # cat(stats)
              # statistic values
              stats <- tidy(nls.model)
              rownames(stats) <- stats[,1]
              stats <- stats[,-1]
              vec <- as.numeric(c(t(stats)))
              names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
              
              vec1 <- vec[grep("1", names(vec))]
              # As we haven't estimated period 2 we add the data to vec2
              vec2 <- vec[grep("2", names(vec))]
              temp <- c(init_value, "NA", "NA", "NA")
              names(temp) <- c("estimate.per2", "std.error.per2"
                               , "statistic.per2", "p.value.per2")
              vec2 <- c(vec2, temp)
              mode(vec2) <- "numeric"
              
              # R (goodness of fit)
              r1 <- cor(gd[which(gd$group==unique(s2c$group)[1]),"data"]
                        , predict(nls.model)[1:length(which(gd$group==unique(s2c$group)[1]))])
              r2 <- cor(gd[which(gd$group==unique(s2c$group)[2]),"data"]
                        , predict(nls.model)[(length(which(gd$group==unique(s2c$group)[1]))+1):length(gd$group)])
              
              # Groups comparison
              K <- rbind(c(-1, 0, 1, 0, 0)
                         , c(0, -1, 0, 1, 0)
              )
              rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1")
              colnames(K) <- names(coef(nls.model))
              
              # group.stats <- tidy(summary(glht(nls.model, linfct = K)))
              group.stats <- as.data.frame(broom::tidy(summary(multcomp::glht(nls.model, linfct = K))))
              rownames(group.stats) <- group.stats[,1]
              coeff.group <- group.stats$lhs
              group.stats <- group.stats[,-c(1:2)]
              gs.vec <- as.numeric(c(t(group.stats)))
              names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
              gs.vec <- c(gs.vec, setNames(as.numeric(c(vec2["estimate.per2"]-vec1["estimate.per1"], NA, NA, NA))
                                           , c("estimate.Per2_Per1", "std.error.Per2_Per1"
                                               , "statistic.Per2_Per1", "p.value.Per2_Per1")))
            }, error = function(err3) {
              convergence_error3 <- grepl("Convergence failure", as.character(err3))
              if(convergence_error3){
                # Fix per 1 and 2
                result = tryCatch({
                  # Fixed per2
                  nls.model=stats::nls(data~
                                  as.numeric(group==unique(s2c$group)[1])*amp1*sin(2*pi*time/init_value+phase1) +
                                  as.numeric(group==unique(s2c$group)[2])*amp2*sin(2*pi*time/init_value+phase2)
                                , algorithm = mode
                                , lower = list(per=min_per)
                                , upper = list(per=max_per)
                                , data = gd
                                , start = c(amp1=1, phase1=0, amp2=1, phase2=0)
                  )
                  # cat(stats)
                  # statistic values
                  stats <- tidy(nls.model)
                  rownames(stats) <- stats[,1]
                  stats <- stats[,-1]
                  vec <- as.numeric(c(t(stats)))
                  names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
                  # As we haven't estimated period 1 we add the data to vec2
                  vec1 <- vec[grep("1", names(vec))]
                  temp <- c(init_value, "NA", "NA", "NA")
                  names(temp) <- c("estimate.per1", "std.error.per1"
                                   , "statistic.per1", "p.value.per1")
                  vec1 <- c(vec1, temp)
                  mode(vec1) <- "numeric"
                  # As we haven't estimated period 2 we add the data to vec2
                  vec2 <- vec[grep("2", names(vec))]
                  temp <- c(init_value, "NA", "NA", "NA")
                  names(temp) <- c("estimate.per2", "std.error.per2"
                                   , "statistic.per2", "p.value.per2")
                  vec2 <- c(vec2, temp)
                  mode(vec2) <- "numeric"
                  
                  # R (goodness of fit)
                  r1 <- cor(gd[which(gd$group==unique(s2c$group)[1]),"data"]
                            , predict(nls.model)[1:length(which(gd$group==unique(s2c$group)[1]))])
                  r2 <- cor(gd[which(gd$group==unique(s2c$group)[2]),"data"]
                            , predict(nls.model)[(length(which(gd$group==unique(s2c$group)[1]))+1):length(gd$group)])
                  
                  # Groups comparison
                  K <- rbind(c(-1, 0, 1, 0)
                             , c(0, -1, 0, 1)
                  )
                  rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1")
                  colnames(K) <- names(coef(nls.model))
                  
                  # group.stats <- tidy(summary(glht(nls.model, linfct = K)))
                  group.stats <- as.data.frame(broom::tidy(summary(multcomp::glht(nls.model, linfct = K))))
                  rownames(group.stats) <- group.stats[,1]
                  coeff.group <- group.stats$lhs
                  group.stats <- group.stats[,-c(1:2)]
                  gs.vec <- as.numeric(c(t(group.stats)))
                  names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
                  gs.vec <- c(gs.vec, setNames(as.numeric(c(vec2["estimate.per2"]-vec1["estimate.per1"], NA, NA, NA))
                                               , c("estimate.Per2_Per1", "std.error.Per2_Per1"
                                                   , "statistic.Per2_Per1", "p.value.Per2_Per1")))
                }, error = function(err4) {
                  convergence_error4 <- grepl("Convergence failure", as.character(err4))
                  if(convergence_error4){
                    # Can't be adjusted to a circadian curve
                    vec1 = vec2 <<- rep("NA", times = length(results.cols))
                    r1 = r2 <<- "NA"
                    gs.vec <<- rep("NA", times = length(group.cols)-1)
                  }# end error 4
                })
              } # end error 3
            })
        }# end error 2
      })
    } # end error 1
    }, finally = {
        res1 <- rbind(c(feature=rownames(data)[gene], as.numeric(vec1), r1 = as.numeric(r1)))
        results@G1[gene,] <- res1
        res2 <- as.numeric(rbind(c(feature=rownames(data)[gene], as.numeric(vec2), r2 = as.numeric(r2))))
        results@G2[gene,] <- res2
        results@Comparison[gene,] <- as.matrix(rbind(c(feature=rownames(data)[gene], as.numeric(gs.vec))))
    })
    if(shiny) incProgress(amount=1/total_genes)
  }
  
  # G1
  # Calc adjusted p.val for period
  results@G1$BH.q.value.per1 <- p.adjust(as.numeric(results@G1$p.value.per1), method="BH")
  results@G1$BH.q.value.amp1 <- p.adjust(as.numeric(results@G1$p.value.amp1), method="BH")
  results@G1$BH.q.value.phase1 <- p.adjust(as.numeric(results@G1$p.value.phase1), method="BH")
  # G2
  results@G2$BH.q.value.per2 <- p.adjust(as.numeric(results@G2$p.value.per2), method="BH")
  results@G2$BH.q.value.amp2 <- p.adjust(as.numeric(results@G2$p.value.amp2), method="BH")
  results@G2$BH.q.value.phase2 <- p.adjust(as.numeric(results@G2$p.value.phase2), method="BH")
  # Comparison
  results@Comparison$BH.q.value.Per2_Per1 <- p.adjust(as.numeric(results@Comparison$p.value.Per2_Per1), method="BH")
  results@Comparison$BH.q.value.Amp2_Amp1 <- p.adjust(as.numeric(results@Comparison$p.value.Amp2_Amp1), method="BH")
  results@Comparison$BH.q.value.Ph2_Ph1 <- p.adjust(as.numeric(results@Comparison$p.value.Ph2_Ph1), method="BH")
  
  
  # if(shiny) incProgress(total_genes/total_genes)
  return(results)
}


# gen_samples <- function(total_genes, duration, interval, per_circ_genes, n_replicates, sdev=1){
#   # We generate randomly ~14k amplitudes between -1,5 to -0.5 and 0.5 to 1.5 and n periods from 10 to 30
#   # total_genes=15000
#   # duration=24
#   # interval=4
#   # per_circ_genes=20
#   # n_replicates=3
#   # sdev=1.5
# 
#   amps <- c(runif(total_genes/2*per_circ_genes, min=-1.5, max=-0.5)
#             , runif(total_genes/2*per_circ_genes, min=0.5, max=1.5))
#   pers <- runif(total_genes*per_circ_genes, min=10, max=30)
#   phases <- runif(total_genes*per_circ_genes, min=-6, max=6)
# 
#   circ = function(t,amp,per,phase){
#     y = amp*sin(2*pi*t/per + phase)
#     return(y)
#   }
# 
#   # Array to store the data
#   temp <- data.frame(replicates = rep(c(1:3), duration/interval), timepoints = rep(seq(1,duration,by=interval), each=3))
#   colnames_array <- apply(temp[,c("timepoints", "replicates")], 1, paste, collapse = "_rep")
#   colnames_array <- paste("S", colnames_array, sep="")
#   in_silico_data <- setNames(data.frame(matrix(ncol = length(colnames_array)+3, nrow = 0))
#                              , c(colnames_array,"amp","per","phase"))
# 
#   # Circa genes
#   for (r in 1:as.integer(total_genes*per_circ_genes)){
#     # Seed (True) value
#     res=as.numeric(unlist(lapply(list(seq(1,duration,by=interval)),circ,amps[r],pers[r],phases[r])))
#     # Generate replicates with deviation from the true value
#     rep1=rep2=rep3=numeric()
#     # Generate the replicates around the seed
#     for (i in 1:length(res)){
#       rep1=c(rep1,res[i]+rnorm(1,mean=res[i],sd=sdev))
#       rep2=c(rep2,res[i]+rnorm(1,mean=res[i],sd=sdev))
#       rep3=c(rep3,res[i]+rnorm(1,mean=res[i],sd=sdev))
#     }
#     long <- data.frame(data=c(rep3,rep1,rep2)
#                        ,ind=rep(c(1:3),each=duration/interval)
#                        ,time=rep(seq(1,duration,by=interval),3))
#     # unlist(unstack(long, data~ind+time))
#     in_silico_data[r,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),amps[r],pers[r],phases[r])
# 
#   }
# 
#   # Non-circa genes
#   for (r in 1:as.integer(total_genes*(1-per_circ_genes))){
#     # Seed (True); after standarization all the means are centered ~0 and all the sd are = 1
#     res=rep(0, times=duration/interval)
#     # Generate replicates with deviation from the true value
#     rep1=rep2=rep3=numeric()
#     # Generate the replicates around the seed
#     for (i in 1:length(res)){
#       rep1=c(rep1,res[i]+rnorm(1,mean=res[i],sd=sdev))
#       rep2=c(rep2,res[i]+rnorm(1,mean=res[i],sd=sdev))
#       rep3=c(rep3,res[i]+rnorm(1,mean=res[i],sd=sdev))
#     }
#     long <- data.frame(data=c(rep3,rep1,rep2)
#                        ,ind=rep(c(1:3),each=duration/interval)
#                        ,time=rep(seq(1,duration,by=interval),3))
#     in_silico_data[nrow(in_silico_data)+1,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),NA,NA,NA)
#   }
# 
#   # Shuffle rows
#   in_silico_data <- in_silico_data[sample(nrow(in_silico_data)),]
#   rownames(in_silico_data) <- c(1:nrow(in_silico_data))
#   return(in_silico_data)
# }

# gen_samples <- function(total_genes, duration, interval, per_circ_genes, n_replicates, sdev=1){
#   # We generate randomly ~14k amplitudes between -1,5 to -0.5 and 0.5 to 1.5 and n periods from 10 to 30
#   # total_genes=15000
#   # duration=24
#   # interval=4
#   # per_circ_genes=20
#   # n_replicates=3
#   # sdev=1.5
# 
#   amps <- c(runif(total_genes/2*per_circ_genes, min=-1.5, max=-0.5)
#             , runif(total_genes/2*per_circ_genes, min=0.5, max=1.5))
#   pers <- runif(total_genes*per_circ_genes, min=10, max=30)
#   phases <- runif(total_genes*per_circ_genes, min=-6, max=6)
# 
#   circ = function(t,amp,per,phase){
#     y = amp*sin(2*pi*t/per + phase)
#     return(y)
#   }
# 
#   # Array to store the data
#   temp <- data.frame(replicates = rep(c(1:3), duration/interval), timepoints = rep(seq(1,duration,by=interval), each=3))
#   colnames_array <- apply(temp[,c("timepoints", "replicates")], 1, paste, collapse = "_rep")
#   colnames_array <- paste("S", colnames_array, sep="")
#   in_silico_data <- setNames(data.frame(matrix(ncol = length(colnames_array)+3, nrow = 0))
#                              , c(colnames_array,"amp","per","phase"))
# 
#   # Circa genes
#   for (r in 1:as.integer(total_genes*per_circ_genes)){
#     # Seed (True) value
#     res=as.numeric(unlist(lapply(list(seq(1,duration,by=interval)),circ,amps[r],pers[r],phases[r])))
#     # Generate replicates with deviation from the true value
#     for (z in 1:n_replicates){assign(paste("rep",z, sep=""), numeric())}
#     # rep1=rep2=rep3=numeric()
#     # Generate the replicates around the seed
#     for (i in 1:length(res)){
#       for (rep in 1:n_replicates){
#         assign(paste("rep",rep, sep=""), c(get(paste("rep",rep, sep="")),res[i]+rnorm(1,mean=res[i],sd=sdev)))
#       }
#     }
#     long <- data.frame(data=c(rep3,rep1,rep2)
#                        ,ind=rep(c(1:3),each=duration/interval)
#                        ,time=rep(seq(1,duration,by=interval),3))
#     # unlist(unstack(long, data~ind+time))
#     in_silico_data[r,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),amps[r],pers[r],phases[r])
# 
#   }
# 
#   # Non-circa genes
#   for (r in 1:as.integer(total_genes*(1-per_circ_genes))){
#     # Seed (True); after standarization all the means are centered ~0 and all the sd are = 1
#     res=rep(0, times=duration/interval)
#     # Generate replicates with deviation from the true value
#     # rep1=rep2=rep3=numeric()
#     for (z in 1:n_replicates){assign(paste("rep",z, sep=""), numeric())}
#     # Generate the replicates around the seed
#     for (i in 1:length(res)){
#       for (rep in 1:n_replicates){
#         assign(paste("rep",rep, sep=""), c(get(paste("rep",rep, sep="")),res[i]+rnorm(1,mean=res[i],sd=sdev)))
#       }
#       # rep1=c(rep1,res[i]+rnorm(1,mean=res[i],sd=sdev))
#       # rep2=c(rep2,res[i]+rnorm(1,mean=res[i],sd=sdev))
#       # rep3=c(rep3,res[i]+rnorm(1,mean=res[i],sd=sdev))
#     }
#     long <- data.frame(data=c(rep3,rep1,rep2)
#                        ,ind=rep(c(1:3),each=duration/interval)
#                        ,time=rep(seq(1,duration,by=interval),3))
#     in_silico_data[nrow(in_silico_data)+1,] <- c(as.numeric(unlist(unstack(long, data~ind+time))),NA,NA,NA)
#   }
# 
#   # Shuffle rows
#   in_silico_data <- in_silico_data[sample(nrow(in_silico_data)),]
#   rownames(in_silico_data) <- c(1:nrow(in_silico_data))
#   return(in_silico_data)
# }


# nls_regression_groups <- function(data, s2c, id_col, sample_col, time_col, group_col, ind_col
#                            , shiny=FALSE, mode="default", init_value=24
#                            , max_per=Inf, min_per=-Inf){
#   
#   # data=as.data.frame(datafr[1:100,])
#   # s2c=s2c
#   # id_col=ncol(data)
#   # sample_col=1
#   # time_col=3
#   # group_col=4
#   # ind_col=2
#   # shiny=FALSE
#   # mode="port"
#   # init_value=24
#   # max_per=Inf
#   # min_per=-Inf
#   
#   
#   # Format s2c so nothing does weird things
#   s2c[,time_col] <- as.numeric(as.character(s2c[,time_col]))
#   s2c[,ind_col] <- as.integer(as.character(s2c[,ind_col]))
#   s2c[,sample_col]<- as.character(s2c[,sample_col])
#   s2c[,group_col] <- as.character(s2c[,group_col])
#   
#   rown <- data[,id_col]
#   data <- data[,-id_col]
#   data[,1:ncol(data)] <- sapply(data[,1:ncol(data)], as.numeric)
#   rownames(data) <- rown
#   
#   G1_name <- unique(s2c[,group_col])[1]
#   G2_name <- unique(s2c[,group_col])[2]
#   data1 <- data[,s2c[which(s2c[,group_col] == G1_name),sample_col]]
#   data2 <- data[,s2c[which(s2c[,group_col] == G2_name),sample_col]]
#   
#   results.cols <- c("estimate.amp","std.error.amp","statistic.amp","p.value.amp"
#                     ,"estimate.phase","std.error.phase","statistic.phase","p.value.phase"
#                     ,"estimate.per" ,"std.error.per","statistic.per", "p.value.per"
#                     , "r")
#   group.cols <- c("feature", "estimate.Amp2_Amp1", "std.error.Amp2_Amp1", "statistic.Amp2_Amp1", "p.value.Amp2_Amp1"
#                   , "estimate.Ph2_Ph1", "std.error.Ph2_Ph1", "statistic.Ph2_Ph1", "p.value.Ph2_Ph1"
#                   , "estimate.Per2_Per1", "std.error.Per2_Per1", "statistic.Per2_Per1", "p.value.Per2_Per1")
#   G1 <- setNames(data.frame(matrix(ncol = length(results.cols)+1, nrow = 0))
#                       , c("feature", paste(results.cols, "1", sep="")))
#   G2 <- setNames(data.frame(matrix(ncol = length(results.cols)+1, nrow = 0))
#                       , c("feature", paste(results.cols, "2", sep="")))
#   comp <- setNames(data.frame(matrix(ncol = length(group.cols), nrow = 0))
#                        , group.cols)
#   results <- new("harmony_data", G1=G1, G2=G2, Comparison=comp
#                  , run.info=list(G1=G1_name, G2=G2_name, id_col=id_col, sample_col=sample_col
#                                  , time_col=time_col, group_col=group_col
#                                  , ind_col=ind_col, mode=mode, init_value=init_value
#                                  , max_per=max_per, min_per=min_per))
#   
#   for(gene in 1:nrow(data)){
#     # Center data
#     data.cent.1 <- (data1[gene,]-mean(as.numeric(data1[gene,])))
#     # names(data.cent.1) <- names(data1[gene,])
#     data.cent.2 <- (data2[gene,]-mean(as.numeric(data2[gene,])))
#     # names(data.cent.2) <- names(data2[gene,])
#     # Merged data
#     merged_data <- cbind(data.cent.1, data.cent.2)
#     rownames(merged_data) <- "data"
#     # Construct df0 (ind must be different between groups)
#     df0 <- merge(t(merged_data), s2c, by.x = 0, by.y = sample_col)
#     # Change names to match formula
#     df0 <- dplyr::rename(df0, group = colnames(s2c)[group_col]
#                          , time = colnames(s2c)[time_col]
#                          , ind = colnames(s2c)[ind_col])
#     
#     df0 <- df0[order(df0$group, df0$time),]
#     df0$data <- as.numeric(df0$data)
#     df0$time <- as.numeric(df0$time)
#     # Group
#     gd <- groupedData(data~time|ind, data=df0)
#     
#     
#     # Possible try
#     result = tryCatch({
#         nls.model=nls(data~
#                         as.numeric(group==unique(s2c$group)[1])*amp1*sin(2*pi*time/per1+phase1) + 
#                         as.numeric(group==unique(s2c$group)[2])*amp2*sin(2*pi*time/per2+phase2)
#                       , algorithm = mode
#                       , lower = list(per=min_per)
#                       , upper = list(per=max_per)
#                       , data = gd
#                       , start = c(amp1=1, phase1=0, per1=init_value, amp2=1, phase2=0, per2=init_value)
#         )
# 
#         # Each group' stats separately
#         # summary(m.str)
#         stats <- tidy(nls.model)
#         rownames(stats) <- stats[,1]
#         coeff <- stats$term
#         stats <- stats[,-1]
#         # Populate G1 and G2
#         vec <- as.numeric(c(t(stats)))
#         names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
#         vec1 <- vec[grep("1", names(vec))]
#         vec2 <- vec[grep("2", names(vec))]
#         # R (goodness of fit)
#         r1 <- cor(gd[which(gd$group==unique(s2c$group)[1]),"data"]
#                   , predict(nls.model)[1:length(which(gd$group==unique(s2c$group)[1]))])
#         r2 <- cor(gd[which(gd$group==unique(s2c$group)[2]),"data"]
#                   , predict(nls.model)[(length(which(gd$group==unique(s2c$group)[1]))+1):length(gd$group)])
#         
#         # Groups comparison
#         K <- rbind(c(-1, 0, 0, 1 , 0, 0)
#                    , c(0, -1, 0, 0, 1, 0)
#                    , c(0, 0, -1, 0, 0, 1)
#         )
#         rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1", "Per2_Per1")
#         colnames(K) <- names(coef(nls.model))
#         
#         # summary(glht(nls.model, linfct = K), test=adjusted("none"))
#         # standard error is just the standard deviation divided by the square root of the sample size.
#         
#         
#         
#         mine <- glht(nls.model, linfct = K)
#         group.stats <- tidy(summary(mine))
#         rownames(group.stats) <- group.stats[,1]
#         coeff.group <- group.stats$lhs
#         group.stats <- group.stats[,-c(1:2)]
#         gs.vec <<- as.numeric(c(t(group.stats)))
#         names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
#     }, error = function(err1) {
#           result2 = tryCatch({ ## Per 1 fixed
#                 # Try with fixed period1
#                 nls.model=nls(data~
#                                 as.numeric(group==unique(s2c$group)[1])*amp1*sin(2*pi*time/init_value+phase1) +
#                                 as.numeric(group==unique(s2c$group)[2])*amp2*sin(2*pi*time/per2+phase2)
#                               , algorithm = mode
#                               , lower = list(per=min_per)
#                               , upper = list(per=max_per)
#                               , data = gd
#                               , start = c(amp1=1, phase1=0, amp2=1, phase2=0, per2=init_value)
#                 )
#                 # cat(stats)
#                 # statistic values
#                 stats <- tidy(nls.model)
#                 rownames(stats) <- stats[,1]
#                 stats <- stats[,-1]
#                 vec <- as.numeric(c(t(stats)))
#                 names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
#                 # As we haven't estimated period 1 we add the data to vec1
#                 vec1 <- vec[grep("1", names(vec))]
#                 temp <- c(init_value, "NA", "NA", "NA")
#                 names(temp) <- c("estimate.per1", "std.error.per1"
#                                  , "statistic.per1", "p.value.per1")
#                 vec1 <<- c(vec1, temp)
#                 mode(vec1) <- "numeric"
#                 vec2 <<- vec[grep("2", names(vec))]
#                 # R (goodness of fit)
#                 r1 <- cor(gd[which(gd$group==unique(s2c$group)[1]),"data"]
#                           , predict(nls.model)[1:length(which(gd$group==unique(s2c$group)[1]))])
#                 r2 <- cor(gd[which(gd$group==unique(s2c$group)[2]),"data"]
#                           , predict(nls.model)[(length(which(gd$group==unique(s2c$group)[1]))+1):length(gd$group)])
#                 
#                 # Groups comparison
#                 K <- rbind(c(-1, 0, 1, 0, 0)
#                            , c(0, -1, 0, 1, 0)
#                 )
#                 rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1")
#                 colnames(K) <- names(coef(nls.model))
#                 
#                 group.stats <- tidy(summary(glht(nls.model, linfct = K)))
#                 rownames(group.stats) <- group.stats[,1]
#                 coeff.group <- group.stats$lhs
#                 group.stats <- group.stats[,-c(1:2)]
#                 gs.vec <- as.numeric(c(t(group.stats)))
#                 names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
#                 gs.vec <- c(gs.vec, setNames(as.numeric(c(vec2["estimate.per2"]-vec1["estimate.per1"], NA, NA, NA))
#                                    , c("estimate.Per2_Per1", "std.error.Per2_Per1"
#                                      , "statistic.Per2_Per1", "p.value.Per2_Per1")))
#                 
#               }, error = function(e) {
#                 result = tryCatch({
#                   # Fixed per2
#                   nls.model=nls(data~
#                                   as.numeric(group==unique(s2c$group)[1])*amp1*sin(2*pi*time/per1+phase1) +
#                                   as.numeric(group==unique(s2c$group)[2])*amp2*sin(2*pi*time/init_value+phase2)
#                                 , algorithm = mode
#                                 , lower = list(per=min_per)
#                                 , upper = list(per=max_per)
#                                 , data = gd
#                                 , start = c(amp1=1, phase1=0, amp2=1, phase2=0, per1=init_value)
#                   )
#                   # cat(stats)
#                   # statistic values
#                   stats <- tidy(nls.model)
#                   rownames(stats) <- stats[,1]
#                   stats <- stats[,-1]
#                   vec <- as.numeric(c(t(stats)))
#                   names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
#                   
#                   vec1 <<- vec[grep("1", names(vec))]
#                   # As we haven't estimated period 2 we add the data to vec2
#                   vec2 <- vec[grep("2", names(vec))]
#                   temp <- c(init_value, "NA", "NA", "NA")
#                   names(temp) <- c("estimate.per2", "std.error.per2"
#                                    , "statistic.per2", "p.value.per2")
#                   vec2 <- c(vec2, temp)
#                   mode(vec2) <- "numeric"
#                  
#                   # R (goodness of fit)
#                   r1 <- cor(gd[which(gd$group==unique(s2c$group)[1]),"data"]
#                             , predict(nls.model)[1:length(which(gd$group==unique(s2c$group)[1]))])
#                   r2 <- cor(gd[which(gd$group==unique(s2c$group)[2]),"data"]
#                             , predict(nls.model)[(length(which(gd$group==unique(s2c$group)[1]))+1):length(gd$group)])
#                   
#                   # Groups comparison
#                   K <- rbind(c(-1, 0, 1, 0, 0)
#                              , c(0, -1, 0, 1, 0)
#                   )
#                   rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1")
#                   colnames(K) <- names(coef(nls.model))
#                   
#                   group.stats <- tidy(summary(glht(nls.model, linfct = K)))
#                   rownames(group.stats) <- group.stats[,1]
#                   coeff.group <- group.stats$lhs
#                   group.stats <- group.stats[,-c(1:2)]
#                   gs.vec <- as.numeric(c(t(group.stats)))
#                   names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
#                   gs.vec <<- c(gs.vec, setNames(as.numeric(c(vec2["estimate.per2"]-vec1["estimate.per1"], NA, NA, NA))
#                                                , c("estimate.Per2_Per1", "std.error.Per2_Per1"
#                                                    , "statistic.Per2_Per1", "p.value.Per2_Per1")))
#                 }, error = function(e) {
#                   # Fix per 1 and 2
#                   result = tryCatch({
#                     # Fixed per2
#                     nls.model=nls(data~
#                                     as.numeric(group==unique(s2c$group)[1])*amp1*sin(2*pi*time/init_value+phase1) +
#                                     as.numeric(group==unique(s2c$group)[2])*amp2*sin(2*pi*time/init_value+phase2)
#                                   , algorithm = mode
#                                   , lower = list(per=min_per)
#                                   , upper = list(per=max_per)
#                                   , data = gd
#                                   , start = c(amp1=1, phase1=0, amp2=1, phase2=0)
#                     )
#                     # cat(stats)
#                     # statistic values
#                     stats <- tidy(nls.model)
#                     rownames(stats) <- stats[,1]
#                     stats <- stats[,-1]
#                     vec <- as.numeric(c(t(stats)))
#                     names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
#                     # As we haven't estimated period 1 we add the data to vec2
#                     vec1 <- vec[grep("1", names(vec))]
#                     temp <- c(init_value, "NA", "NA", "NA")
#                     names(temp) <- c("estimate.per1", "std.error.per1"
#                                      , "statistic.per1", "p.value.per1")
#                     vec1 <- c(vec1, temp)
#                     mode(vec1) <- "numeric"
#                     # As we haven't estimated period 2 we add the data to vec2
#                     vec2 <- vec[grep("2", names(vec))]
#                     temp <- c(init_value, "NA", "NA", "NA")
#                     names(temp) <- c("estimate.per2", "std.error.per2"
#                                      , "statistic.per2", "p.value.per2")
#                     vec2 <- c(vec2, temp)
#                     mode(vec2) <- "numeric"
#                     
#                     # R (goodness of fit)
#                     r1 <- cor(gd[which(gd$group==unique(s2c$group)[1]),"data"]
#                               , predict(nls.model)[1:length(which(gd$group==unique(s2c$group)[1]))])
#                     r2 <- cor(gd[which(gd$group==unique(s2c$group)[2]),"data"]
#                               , predict(nls.model)[(length(which(gd$group==unique(s2c$group)[1]))+1):length(gd$group)])
#                     
#                     # Groups comparison
#                     K <- rbind(c(-1, 0, 1, 0)
#                                , c(0, -1, 0, 1)
#                     )
#                     rownames(K) <- c("Amp2_Amp1", "Ph2_Ph1")
#                     colnames(K) <- names(coef(nls.model))
#                     
#                     group.stats <- tidy(summary(glht(nls.model, linfct = K)))
#                     rownames(group.stats) <- group.stats[,1]
#                     coeff.group <- group.stats$lhs
#                     group.stats <- group.stats[,-c(1:2)]
#                     gs.vec <- as.numeric(c(t(group.stats)))
#                     names(gs.vec) <- c(outer(colnames(group.stats), rownames(group.stats), paste, sep = "."))
#                     gs.vec <- c(gs.vec, setNames(as.numeric(c(vec2["estimate.per2"]-vec1["estimate.per1"], NA, NA, NA))
#                                                  , c("estimate.Per2_Per1", "std.error.Per2_Per1"
#                                                      , "statistic.Per2_Per1", "p.value.Per2_Per1")))
#                   }, error = function(e) {
#                     # Can't be adjusted to a circadian curve
#                     vec1 = vec2 <<- rep("NA", times = length(results.cols))
#                     r1 = r2 <<- "NA"
#                     gs.vec <<- rep("NA", times = length(group.cols)-1)
#                   })
#                 })
#               })
#     }, finally = {
#       res1 <- rbind(c(feature=rownames(data)[gene], vec1, r1 = r1))
#       results@G1[gene,] <- res1
#       res2 <- rbind(c(feature=rownames(data)[gene], vec2, r2 = r2))
#       results@G2[gene,] <- res2
#       results@Comparison[gene,] <- as.matrix(rbind(c(feature=rownames(data)[gene], gs.vec)))
#     })
#     if(shiny) incProgress(amount=1/total_genes)
#   }
#   
#   # G1
#   # Calc adjusted p.val for period
#   results@G1$BH.q.value.per1 <- p.adjust(as.numeric(results@G1$p.value.per1), method="BH")
#   results@G1$BH.q.value.amp1 <- p.adjust(as.numeric(results@G1$p.value.amp1), method="BH")
#   results@G1$BH.q.value.phase1 <- p.adjust(as.numeric(results@G1$p.value.phase1), method="BH")
#   # G2
#   results@G2$BH.q.value.per2 <- p.adjust(as.numeric(results@G2$p.value.per2), method="BH")
#   results@G2$BH.q.value.amp2 <- p.adjust(as.numeric(results@G2$p.value.amp2), method="BH")
#   results@G2$BH.q.value.phase2 <- p.adjust(as.numeric(results@G2$p.value.phase2), method="BH")
#   # Comparison
#   results@Comparison$BH.q.value.Per2_Per1 <- p.adjust(as.numeric(results@Comparison$p.value.Per2_Per1), method="BH")
#   results@Comparison$BH.q.value.Amp2_Amp1 <- p.adjust(as.numeric(results@Comparison$p.value.Amp2_Amp1), method="BH")
#   results@Comparison$BH.q.value.Ph2_Ph1 <- p.adjust(as.numeric(results@Comparison$p.value.Ph2_Ph1), method="BH")
#   
#   
#   # Calculate moderated values
#   # amp_p90 <- quantile(as.numeric(results@Comparison$std.error.Amp2_Amp1), probs=0.9, na.rm=TRUE)
#   # per_p90 <- quantile(as.numeric(results@Comparison$std.error.Per2_Per1), probs=0.9, na.rm=TRUE)
#   # phase_p90 <- quantile(as.numeric(results@Comparison$std.error.Phase2_Pase1), probs=0.9, na.rm=TRUE)
#   # # glht calculates an adjusted p value of the statistic, however, as we want to add the p90
#   # # we need to calculate the the statistic ourselves:
#   # 
#   # moderated_t_amp <- results@Comparison$estimate.Amp2_Amp1/(results@Comparison$std.error.Amp2_Amp1+amp_p90)
#   # normal_t <- results@Comparison$estimate.Amp2_Amp1/results@Comparison$std.error.Amp2_Amp1
#   # p.adjust(moderated_t_amp,method="")
#   
#   
#   # re
#   
#   # if(shiny) incProgress(total_genes/total_genes)
#   return(results)
# }