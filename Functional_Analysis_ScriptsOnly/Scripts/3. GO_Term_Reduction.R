library(org.Hs.eg.db)
library(GO.db)

## load enrichment output ## 

## Note: replace diversity order, e.g. Q1.00, with desired diversity order ##
## throughout this R code                                                  ##

data <- read.csv("Output/GSEA_analysis/GSEA_Q1.00.csv")

rownames(data) <- data$ID


immune_trees <- c("GO:0002376", "GO:0002682",
                  "GO:0002764", "GO:0001816", "GO:0019221", "GO:0006952",
                  "GO:0031347", "GO:0070661")

children <- as.list(GOBPOFFSPRING[immune_trees])

immune_system_processes <- c(children$`GO:0002376`, "GO:0002376")
regulation_of_immune_system_processes <- c(children$`GO:0002682`, "GO:0002682")
immune_response_regulating_signaling_pathway <- c(children$`GO:0002764`, "GO:0002764")
cytokine_production <- c(children$`GO:0001816`, "GO:0001816")
cytokine_mediated_signaling_pathway <- c(children$`GO:0019221`, "GO:0019221")
defense_response <- c(children$`GO:0006952`, "GO:0006952")
regulation_of_defense_response <- c(children$`GO:0031347`, "GO:0031347")
leukocyte_proliferation <- c(children$`GO:0070661`, "GO:0070661")

## Generating new columns ##

data$immune_system_processes <- data$X
data$regulation_of_immune_system_processes <- data$X
data$immune_response_regulating_signaling_pathway <- data$X
data$cytokine_production <- data$X
data$cytokine_mediated_signaling_pathway <- data$X
data$defense_response <- data$X
data$regulation_of_defense_response <- data$X
data$leukocyte_proliferation <- data$X

for (ID in data$ID){
  data[ID,]$immune_system_processes <- ID%in% immune_system_processes 
}


for (ID in data$ID){
  data[ID,]$regulation_of_immune_system_processes <- ID%in% regulation_of_immune_system_processes
}

for (ID in data$ID){
  data[ID,]$immune_response_regulating_signaling_pathway <- ID%in% immune_response_regulating_signaling_pathway
}


for (ID in data$ID){
  data[ID,]$cytokine_production <- ID%in% cytokine_production
}

for (ID in data$ID){
  data[ID,]$cytokine_mediated_signaling_pathway <- ID%in% cytokine_mediated_signaling_pathway
}

for (ID in data$ID){
  data[ID,]$defense_response <- ID%in% defense_response
}

for (ID in data$ID){
  data[ID,]$regulation_of_defense_response <- ID%in% regulation_of_defense_response
}

for (ID in data$ID){
  data[ID,]$leukocyte_proliferation <- ID%in% leukocyte_proliferation
}

write.csv(data, "Output/GSEA_analysis/gsea_simplified_ImmuneTerms_Q1.00.csv")

