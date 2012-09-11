###IBI Calculator
Wetland_IBI <- function(data, DistinctCode=F){
  ###Pull in outside data###
  source("IBIname_match_wetlands2.r")
  load("ibiv3.RData")
  ibi <- ibiv3

  data <- IBIname_match(data)
  colnames(data)[which(colnames(data) == "FinalID")] <- "Taxa"
  colnames(data)[which(colnames(data) == "BAResult")] <- "Result"
  data$SampleID <- as.factor(data$SampleID)
  ###Calculate total count###
  total_count <- daply(data, "SampleID", function(df)sum(df$Result))
  ###Create sample count flag###
  sample_count_flag <- rep(NA, length(total_count))
  names(sample_count_flag) <- names(total_count)
  sample_count_flag[(which(total_count>=500))] <- "Adequate"
  sample_count_flag[(which(total_count < 500 & total_count >= 450))] <- "Within specifications"
  sample_count_flag[(which(total_count < 450))] <- "Inadequate"
  ###Subsample down to 500###
  datalength <- length(data)
  rarifydown <- function(data){unlist(sapply(unique(data$SampleID), function(sample){
    v <- data[data$SampleID==sample, "Result"]
    
    if(sum(v)>=500){rrarefy(v, 500)} else
    {v}
  }
  )
  )
  }
  
  require(doParallel)
  require(vegan)
  registerDoParallel()
  rarificationresult <- foreach(i=1:20, .combine=cbind, .packages="vegan") %dopar% {
    rarifydown(data)
  }
  data <- cbind(data, rarificationresult)
  colnames(data)[(datalength + 1):(datalength + 20)]<- paste("Replicate", 1:20)
  
  ###Metrics set up###
  metrics <- as.data.frame(matrix(NA, nrow = length(unique(data$SampleID)), ncol = 140))
  samplenames <- names(tapply(data$STE, data$SampleID, length))
  data$STE <- as.character(data$STE)
  ###Merge revelant ibi data into data table###
  data$MaxTol <- ibi$MaxTol[match(data$Taxa, ibi$FinalID)]
  data$MaxTol <- as.numeric(data$MaxTol)
  data$Class <- ibi$Class[match(data$Taxa, ibi$FinalID)]
  data$Order <- as.character(ibi$Order[match(data$Taxa, ibi$FinalID)])
  data$FunctionalFeedingGroup <- as.character(ibi$FunctionalFeedingGroup[match(data$Taxa, ibi$FinalID)])
  data$Family <- as.factor(ibi$Family[match(data$Taxa, ibi$FinalID)])   
  data$subfamily <- ibi$subfamily[match(data$Taxa, ibi$FinalID)]
  metrics <<- data
  for(i in 1:20){
    ###Number of Oligochaeta taxa###
    metrics[[i]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                          function(d)length(unique(d$STE[d$Class=="Oligochaeta"])))[, 2]			   
    
    ###Numer of EOT taxa
    metrics[[i+20]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                             function(d)length(unique(d$STE[d$Order %in% c("Ephemeroptera", "Odonata", "Trichoptera")])))[, 2]

    ###Number of Scraper taxa###
    metrics[[i+40]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                             function(d)length(unique(d$STE[which(d$FunctionalFeedingGroup == "SC")])))[, 2]
    
    ###Number of Predator taxa###
    metrics[[i+60]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                             function(d)length(unique(d$STE[which(d$FunctionalFeedingGroup == "P")])))[, 2]	
    
    ###Percent Coleoptera###
    metrics[[i+80]] <- ddply(data[data[[datalength + i]]>0, ], "SampleID",
                             function(d){
                               100*sum(d$Result[which(d$Order == "Coleoptera")], na.rm=T)/sum(d$Result, na.rm=T)
                             })[, 2]
    
    ###Percent EOT###
    metrics[[i+100]] <- ddply(data[data[[datalength + i]]>0, ], "SampleID",
                              function(d){
                                100*sum(d$Result[d$Order %in% c("Ephemeroptera", "Odonata", "Trichoptera")], na.rm=T)/sum(d$Result, na.rm=T)
                              })[, 2]

    ###Percent Tanypodinae/Chironomidae###
    metrics[[i+120]] <- ddply(data[data[[datalength + i]]>0, ], "SampleID",
                              function(d){
                                ifelse(length(which(d$Family == "Chironomidae"))>0,
                                       100*sum(d$Result[d$subfamily == "Tanypodinae"], na.rm=T)/sum(d$Result[d$Family == "Chironomidae"], na.rm=T),
                                0)
                              })[, 2]
    
    ###Percent Dominant###
    metrics[[i+140]] <- ddply(data[data[[datalength + i]]>0, ], "SampleID",
                              function(d){
                                100*sum(d$Result[order(d$Result, decreasing=T)][1:3])/sum(d$Result, na.rm=T)
                              })[, 2]
  }
  ###Convert metrics to scores###
  scores <- metrics
  for(i in 1:20){
    ###Oligochaeta###
    scores[, i] <- cut(metrics[, i], breaks=c(-1:6), labels=c(10, 9, 7, 5, 3, 2, 0))
    ###Numer of EOT taxa
    scores[, i+20] <- cut(metrics[, i+20], breaks=c(-1, 0, 2, 3, 4, 5, 100), labels=as.factor(c(0, 2, 4, 6, 8, 10)))
    ###Number of Scraper taxa###
    scores[, i+40] <- cut(metrics[, i+40], breaks=c(-1, 0, 1, 100), labels=c(0, 5, 10))
    ###Number of Predator taxa###
    scores[, i+60] <- cut(metrics[, i+60], breaks=c(-1, 2, 4, 3, 5, 6, 7, 9, 100), labels=c(0, 1, 3, 4, 6, 7, 9, 10))
    ###Percent Coleoptera###
    metrics[, i+80] <- round(metrics[, i+80], digits=1)
    scores[, i+80] <- cut(metrics[, i+80], breaks=c(-1, .1, .3, .4, .6, .8, 1, 1.2, 1.4, 1.5, 1.7, 100), labels=0:10)
    ###Percent EOT###
    metrics[, i+100] <- round(metrics[, i+100], digits=1)
    scores[, i+100] <- cut(metrics[, i+100], breaks=c(-1, 3, 6.1, 9.2, 12.3, 15.4, 18.5, 21.6, 24.7, 27.8, 30.9, 100), labels=0:10)
    ###Percent Tanypodinae/Chironomidae###
    metrics[, i+120] <- round(metrics[, i+120], digits=1)
    scores[, i+120] <- cut(metrics[, i+120], breaks=c(-1, 6.2, 12.4, 18.7, 24.9, 31.2, 37.4, 43.7, 50, 56.2, 62.5, 100), labels=0:10)
    ###Percent Dominant###
    metrics[, i+140] <- round(metrics[, i+140], digits=1)
    scores[, i+140] <- cut(metrics[, i+140], breaks=c(-1, 65.3, 67.5, 69.8, 72, 74.2, 76.4, 78.6, 80.9, 83.1, 87.5, 100), labels=0:10)
  }
  ###IBI###
  for(i in 0:19){
    scores[[i+161]] <- sapply(1:length(unique(data$SampleID)), function(j)(10/7)*(sum(c(scores[j, 1+i], scores[j, 20+i], 
                                                                                        scores[j, 40+i], scores[j, 60+i], scores[j, 80+i], scores[j, 100+i], scores[j, 120+i], scores[j, 140+i]))))
  }
  for(i in 1:ncol(scores)){
    scores[, i] <- as.numeric(as.character(scores[, i]))
  }
  ###Calculate means for metrics and scores###
  means <- as.data.frame(matrix(NA, nrow=length(unique(data$SampleID)), ncol = 17))
  for(i in 1:8){
    means[[i]] <- apply(metrics[, (((i-1)*20)+1):(20*i)], 1, function(d)sum(d)/20)
  }  
  for(i in 1:9){
    means[[i+8]] <- apply(scores[, (((i-1)*20)+1):(20*i)], 1, mean)
  }
  ###Construct output frame###
  results <- as.data.frame(matrix(NA, nrow=length(unique(data$SampleID)), ncol = 22))
  results[[1]] <- unique(data[, c("StationCode", "SampleID")])$StationCode
  results[[2]] <- unique(data$SampleID)
  results <- results[match(samplenames, as.character(results[[2]])),]
  results[[3]] <- total_count[!is.na(total_count)]
  results[[4]] <- sample_count_flag[!is.na(sample_count_flag)]
  results[[5]] <- rep(20, times=length(unique(data$SampleID)))
  results[which(results[,3] <= 500), 5] <- 1
  results[, 6:22] <- means
  results[, 22] <- round(results[, 22], digits=2)
  colnames(results) <- c("StationCode", "SampleID", "Total Count", "Count Flag", "Number of Iteration", 
                         "Number of Oligochaeta taxa", "Number of EOT Taxa", "Number of Scraper taxa", "Number of Predator Taxa", 
                         "Percent Coleoptera", "Percent EOT", "Percent Tanypodinae to Chironomidae", "Percent Dominant",
                         "Oligochaeta Score", "EOT Taxa Score", "Scraper Score", "Predator Taxa Score", 
                         "Coleoptera Score", "Percent EOT Score", "Percent Tanypodinae to Chironomidae Score", "Percent Dominant Score", "IBI Score")
  ###Return results###
  return(results)
}

###Compute IBI score###
dat <- read.csv("Query2.csv")
results <- Wetland_IBI(dat, DistinctCode=T)

###Write results to csv###
write.csv(results, file="Wetlands_IBI_metrics.csv")
