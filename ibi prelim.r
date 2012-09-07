###scripts###
source("ibi_misc.r")

###Import Data###
data <- read.csv("BMI.csv")
colnames(data)[1] <- "SampleID"
data$SampleID <- as.factor(data$SampleID)
data <- data[!is.na(data$BAResult), ]
data <- data[!is.na(data$FinalID), ]

###Fix non-matches###
data$FinalID <- as.character(data$FinalID)
data[data$FinalID == "Simocephalus sp.", "FinalID"] <- "Simocephalus"
data[data$FinalID == "Alona sp.", "FinalID"] <- "Alona"
data[data$FinalID == "Ceriodaphnia sp.", "FinalID"] <- "Ceriodaphnia"
data[data$FinalID == "Ilyocryptus sp.", "FinalID"] <- "Ilyocryptus"
data[data$FinalID == "Anax sp.", "FinalID"] <- "Anax"
data[data$FinalID == "Daphnia sp. (DISTINCT)", "FinalID"] <- "Daphnia"
data[data$FinalID == "Ishnura", "FinalID"] <- "Ischnura"
data <- data[data$FinalID != "", ]

###Match names###
data <- IBIname_match(data)
data <- data[!duplicated(data[, c("FinalID", "SampleID", "LifeStageName", "Distinct")]),]

###print non-matches###
unique(data[data$SAFIT2=="Missing", "FinalID"])

###Subsample down to 500##

data$subsample <- rarify(data, sample.ID="SampleID", abund="BAResult", subsiz=500)$BAResult

###Scraper Richness###


data$feed <- rep(NA, length(data$BAResult))

load("ibiv2.RData")
ibi <- IBIv2
ibi$FunctionalFeedingGroup <- as.character(ibi$FunctionalFeedingGroup)
ibi$FunctionalFeedingGroup[which(ibi$FunctionalFeedingGroup == "")] <- "None reported"
ibi$FunctionalFeedingGroup <- toupper(ibi$FunctionalFeedingGroup)

data$feed <- ibi[match(data$FinalID, ibi$FinalID), "FunctionalFeedingGroup"]

scraper_richness <- function(i)length(unique(data$SAFIT2[which(i & data$feed[i]=="SC")]))

#Subsampled to 500; SAFIT Level 1; quantitative data only#
samples1 <- data$SampleID[data[["subsample"]]>0 & data$distinct=="Distinct" & data$BenthicResultType == "Quantitative"]
scraper1 <- tapply(1:length(samples1), samples1, scraper_richness)
scraper1
#Unsubsampled; SAFIT Level 1; quantitative data only#
samples2 <- data$SampleID[data[["BAResult"]]>0 & data$distinct=="Distinct" & data$BenthicResultType == "Quantitative"]
scraper2 <- tapply(1:length(samples2), samples2, scraper_richness)
scraper2
#Unsubsampled; SAFIT Level 1; quantitative and rare data#
samples3 <- data$SampleID[data[["BAResult"]]>0 & data$distinct=="Distinct"]
scraper3 <- tapply(1:length(samples3), samples3, scraper_richness)
scraper3
###EOT Richness###
EOTrichness <- function(i)sum(unique(data$SAFIT2[i]) %in% ibi$SAFIT2[ibi$Order %in% c("Ephemeroptera", "Odonata", "Trichoptera")])
#Subsampled to 500; SAFIT Level 1; quantitative data only#
EOT1 <- tapply(1:length(samples1), samples1, EOTrichness)
EOT1
#Unsubsampled; SAFIT Level 1; quantitative data only#
EOT2 <- tapply(1:length(samples2), samples2, EOTrichness)
EOT2
#Unsubsampled; SAFIT Level 1; quantitative and rare data#
EOT3 <- tapply(1:length(samples3), samples3, EOTrichness)
EOT3

###Oligochaeta richness###
Oligochaeta_richness <- function(i)sum(unique(data$SAFIT2[i]) %in% ibi$SAFIT2[ibi$Class == "Oligochaeta"])
#Subsampled to 500; SAFIT Level 1; quantitative data only#
Oligo1 <- tapply(1:length(samples1), samples1, Oligochaeta_richness)
Oligo1
#Unsubsampled; SAFIT Level 1; quantitative data only#
Oligo2 <- tapply(1:length(samples2), samples2, Oligochaeta_richness)
Oligo2
#Unsubsampled; SAFIT Level 1; quantitative and rare data#
Oligo3 <- tapply(1:length(samples3), samples3, Oligochaeta_richness)
Oligo3

###Predator Richness###
predator_richness <- function(i)length(unique(data$SAFIT2[which(i & data$feed[i]=="P")]))
#Subsampled to 500; SAFIT Level 1; quantitative data only#
Pred1 <- tapply(1:length(samples1), samples1, predator_richness)
Pred1
#Unsubsampled; SAFIT Level 1; quantitative data only#
Pred2 <- tapply(1:length(samples2), samples2, predator_richness)
Pred2
#Unsubsampled; SAFIT Level 1; quantitative and rare data#
Pred3 <- tapply(1:length(samples3), samples3, predator_richness)
Pred3

###Percent Coleoptera###
Coleoptera_abund <- function(i, count)sum(data[which(data$SAFIT2[i] %in% ibi$SAFIT2[ibi$Order == "Coleoptera"]), count])
totalabund <- function(s, count)tapply(data[which(data$SampleID %in% s), count], data$SampleID[which(data$SampleID %in% s)], sum)
#Subsampled to 500; quantitative data only#
samples4 <- data$SampleID[data[["subsample"]]>0 & data$BenthicResultType == "Quantitative"]
coleoptera1 <- 100*(tapply(1:length(samples4), samples4, Coleoptera_abund, count="subsample")/totalabund(samples4, "subsample"))
coleoptera1
#Unsubsampled; quantitative data only#
samples5 <- data$SampleID[data[["BAResult"]]>0 & data$BenthicResultType == "Quantitative"]
coleoptera2 <- 100*(tapply(1:length(samples5), samples5, Coleoptera_abund, count="BAResult")/totalabund(samples5, "BAResult"))
coleoptera2
#Unsubsampled; quantitative data and rare data#
samples6 <- data$SampleID[data[["BAResult"]]>0]
coleoptera3 <- 100*(tapply(1:length(samples6), samples6, Coleoptera_abund, count="BAResult")/totalabund(samples6, "BAResult"))
coleoptera3

###Percent EOT###
EOTabund <- function(i, count)sum(data[which(data$SAFIT2[i] %in% ibi$SAFIT2[ibi$Order %in% c("Ephemeroptera", "Odonata", "Trichoptera")]), count])
#Subsampled to 500; quantitative data only#
EOTpercent1 <- 100*(tapply(1:length(samples4), samples4, EOTabund, count="subsample")/totalabund(samples4, "subsample"))
EOTpercent1
#Subsampled to 500; quantitative data only#
EOTpercent2 <- 100*(tapply(1:length(samples5), samples5, EOTabund, count="BAResult")/totalabund(samples5, "BAResult"))
EOTpercent2
#Subsampled to 500; quantitative data and rare data#
EOTpercent3 <- 100*(tapply(1:length(samples6), samples6, EOTabund, count="BAResult")/totalabund(samples6, "BAResult"))
EOTpercent3

###Percent Dominant###
dominate <- function(i)sum(sort(data[i, "BAResult"], decreasing=T)[1:3])
percentdominant <- 100*(tapply(1:length(data$SampleID), data$SampleID, dominate)/totalabund(data$SampleID, "BAResult"))

###Percent Tanypodinae/Chironomidae###
load("taxonomy.RData")
data <- merge(data, taxonomy[,c("FinalID", "Family", "Subfamily")], all.x=T)
TanyChiro <- 100*(tapply(1:length(data$FinalID), data$SampleID,
                    function(i)sum(data$BAResult[which(data$Subfamily[i]=="Tanypodinae")]))/
                      tapply(1:length(data$FinalID), data$SampleID,
                             function(i)sum(data$BAResult[which(data$Family[i]=="Chironomidae")])))
TanyChiro[is.na(TanyChiro)] <- 0 


results <- cbind(percentdominant, TanyChiro, coleoptera1, EOTpercent1, scraper1,
                 EOT1, Oligo1, Pred1)
colnames(results) <- c("Percent three dominant taxa", "Percent Tanypodinae/ Chironomidae",
                       "Percent Coleoptera", "Percent EOT", "Scraper richness", 
                       "EOT richness", "Oligochaeta richness", "Predator richness")
results




