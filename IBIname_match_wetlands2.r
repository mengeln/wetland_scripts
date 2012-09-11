IBIname_match <- function(data){
  colnames(data)[which(colnames(data) == "FinalID")] <- "Taxa"
  colnames(data)[which(colnames(data) == "BAResult")] <- "Result"
  data <- data[which(!is.na(data$Result)), ]
  load("ibiv3.RData")
  require(plyr)
  ibi <- idata.frame(ibiv3)
  load("taxonomy_v4.RData")
  taxonomy <- idata.frame(taxonomy_v4)
  
  ###Aggregate taxa###
  data <- ddply(data, "SampleID", function(df){
    ddply(df, "Taxa", function(sdf){
      id <- unique(sdf[, !(colnames(sdf) %in% "Result")])
      Result <- sum(sdf$Result)
      cbind(id, Result)
    })
  })
  
  ###Match to STE###
  data$STE <- rep(NA, length(data$Taxa))
  data$STE <- ibi$CustomSTE[match(data$Taxa, ibi$FinalID)]
  data$STE <- as.character(data$STE)
  data$STE[which(is.na(data$STE))] <- "Missing"
  
  ###Determine Distinctiveness###
  distinctsorter <- function(taxon, data){
      level <- taxonomy$TaxonomicLevelCode[match(taxon, taxonomy$FinalID)] 
      levelname <- as.character(taxonomy$TaxonomicLevelName[match(taxon, taxonomy$FinalID)])
      samelevel <- taxonomy$FinalID[which(taxonomy[, levelname] == taxon)]
      matchedlevel <- data$STE %in% ibi$CustomSTE[match(samelevel, ibi$FinalID)]
      result <- taxonomy$TaxonomicLevelCode[match(data$STE[matchedlevel], taxonomy$FinalID)] > level
      length(which(result)) != 0
    }
  
  distinctlist <- dlply(data, "SampleID", function(df){
    sapply(1:nrow(df), function(i){
         ifelse(distinctsorter(df$STE[i], df), "Non-Distinct", "Distinct")
  })})
  data$distinct <- unlist(distinctlist)
  ###Override##
  data$DistinctCode <- as.character(data$DistinctCode)
  data$distinct[which(data$distinct == "Non-distinct" & data$DistinctCode == "Yes")] <- "Distinct" 
  ###Return##
  data
}
                         