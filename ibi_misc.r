IBIname_match <- function(data, DistinctCode=F){
  colnames(data)[colnames(data)=="FinalID"] <- "Taxa"
  load("ibiv2.RData")
  ibi <- IBIv2
  ibi$FinalID <- as.character(ibi$FinalID)
  colnames(ibi)[1]<-"Taxa"
  ibi$SAFIT2 <- as.character(ibi$SAFIT2)
  ibi$SAFIT2 <- as.character(ibi$SAFIT2)
  load("taxonomy.RData")

  ###Fix extra spaces
  data$Taxa <- as.character(data$Taxa)
  complex <- grep("Group", data$Taxa)
  extraspace <- data$Taxa[intersect(which(data$SAFIT == "Missing"), which(!(data$Taxa %in% (data$Taxa[complex]))))]
  data$Taxa[intersect(which(data$SAFIT == "Missing"), which(!(data$Taxa %in% (data$Taxa[complex]))))] <- 
    gsub("(\\w+)\\s+$", "\\1", extraspace)
  ###Fix extra caps###
#   cap1 <- gsub("(^\\w)[[:alnum:][:space:]]+", "\\1", data$Taxa)
#   cap2 <- gsub("(^\\w)(\\w+)", "\\2", data$Taxa)
#   data$Taxa <- paste0(cap1, tolower(cap2))
  ###Convert FinalID to SAFIT2###
  data$SAFIT2 <- ibi[match(data$Taxa, ibi$Taxa), "SAFIT2"]
  data$SAFIT2[is.na(data$SAFIT2)] <- "Missing"
  ###Exclude missing from distinct coding###
  data$distinct <- rep(NA, length(data$Taxa))
  data$distinct[which(data$SAFIT2 == "Missing")] <- "Missing"
  ###Determine whether the rest of the FinalIDs are distinct###
  todetermine <- as.character(data$SAFIT2[which(is.na(data$distinct))])
  todetermine <- as.data.frame(cbind(as.character(data$SampleID[which(is.na(data$distinct))]), todetermine))
  todeterminebystation <- tapply(as.character(todetermine$todetermine), todetermine$V1, list)
  for(j in 1:length(todeterminebystation)){
    data$distinct[intersect(which(data$SampleID %in% names(todeterminebystation[j])), which(data$SAFIT2 %in% todeterminebystation[[j]]))] <- sapply(1:length(todeterminebystation[[j]]), function(i, determine, tax){
      index <- which(tax$FinalID == todeterminebystation[[j]][i])
      level <- as.numeric(tax[index, "TaxonomicLevelCode"])
      levelname <- tax[index, "TaxonomicLevelName"]
      if(level >= 60){"Distinct"} else{
        criteron1 <- which(taxonomy$FinalID %in% todeterminebystation[[j]][which(!(unlist(todeterminebystation[j]) %in% unlist(todeterminebystation[[j]][i])))])
        criteron2 <- which(tax[criteron1, levelname] == tax[index, levelname])
        criteron3 <- which(tax[criteron2, "TaxonomicLevelCode"] > tax[index, "TaxonomicLevelCode"])
        if(length(criteron3) > 0){"Non-distinct "} else
        {"Distinct"}}}, determine=todeterminebystation[[j]], tax=taxonomy)}
  if(DistinctCode == T){
    data[data$distinct == "Not Distinct" & data$DistinctCode == "Yes", "distinct"] <- "Distinct"
  }
  colnames(data)[colnames(data)=="Taxa"] <- "FinalID"
  return(data)
}

rarify<-function(inbug, sample.ID, abund, subsiz){
  start.time=proc.time();
  outbug<-inbug;
  sampid<-unique(inbug[,sample.ID]);
  nsamp<-length(sampid);
  #parameters are set up;
  #zero out all abundances in output data set;
  outbug[,abund]<-0;
  #loop over samples, rarify each one in turn;
  
  for(i in 1:nsamp) { ;
                     #extract current sample;
                     isamp<-sampid[i];
                     flush.console();
                     onesamp<-inbug[inbug[,sample.ID]==isamp,];
                     onesamp<-data.frame(onesamp,row.id=seq(1,dim(onesamp)[[1]])); #add sequence numbers as a new column;
                     #expand the sample into a vector of individuals;
                     samp.expand<-rep(x=onesamp$row.id,times=onesamp[,abund]);
                     nbug<-length(samp.expand); #number of bugs in sample;
                     #vector of uniform random numbers;
                     ranvec<-runif(n=nbug);
                     #sort the expanded sample randomly;
                     samp.ex2<-samp.expand[order(ranvec)];
                     #keep only the first piece of ranvec, of the desired fised count size;
                     #if there are fewer bugs than the fixed count size, keep them all;
                     if(nbug>subsiz){subsamp<-samp.ex2[1:subsiz]} else{subsamp<-samp.ex2};
                     #tabulate bugs in subsample;
                     subcnt<-table(subsamp);
                     #define new subsample frame and fill it with new reduced counts;
                     newsamp<-onesamp;
                     newsamp[,abund]<-0;
                     newsamp[match(newsamp$row.id,names(subcnt),nomatch=0)>0,abund]<-as.vector(subcnt);
                     outbug[outbug[,sample.ID]==isamp,abund]<-newsamp[,abund];
  }; #end of sample loop;
  
  outbug; #return subsampled data set as function value;
}; #end of function;