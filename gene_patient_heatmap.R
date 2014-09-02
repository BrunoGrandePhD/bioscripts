# gene_patient_heatmap.R
# 2014-08-11
# Rebecca Lea Johnston and Bruno Grande
# Usage:
#     R --no-save --args ~/Desktop/gene_patient_heatmap/patients.txt
# Arguments:
#     1) File containing two tab-separated columns: patient name and file path
#        to MAF file

# Required Libraries ------------------------------------------------------

library(ggplot2)
library(plyr)
library(RColorBrewer)

# Argument Parsing --------------------------------------------------------

#args <- commandArgs(TRUE)

# patients_file <- args[1]

cohort_file <- file.path("/Users/bgrande/Desktop/gene_patient_heatmap/patients.txt")

n_plotted_genes <- 20

# Building the cohort_mutations data frame --------------------------------

cohort <- read.delim(cohort_file, header = FALSE, 
                     col.names = c("patient_name", "path") )

# str(cohort)

cohort_mutations <- data.frame()

# patient_name <- cohort[3, "patient_name"]
# print(patient_name)
# patient_maf <- read.delim(file.path(cohort$path[3]))
# head(patient_maf)
# wanted_columns <- patient_maf[, c("Hugo_Symbol", "Entrez_Gene_Id", 
#                                       "Variant_Classification")]
# wanted_columns["patient_name"] <- patient_name
# head(wanted_columns)
# cohort_mutations <- rbind(cohort_mutations, wanted_columns)
# head(cohort_mutations)

for (i in 1:nrow(cohort)) {
  patient_name <- cohort[i, "patient_name"]
  patient_maf <- read.delim(file.path(cohort$path[i]))
  wanted_columns <- patient_maf[, c("Hugo_Symbol", "Entrez_Gene_Id", 
                                        "Variant_Classification")]
  wanted_columns["patient_name"] <- patient_name
  cohort_mutations <- rbind(cohort_mutations, wanted_columns)
}


# str(cohort_mutations)
# head(cohort_mutations)
# tail(cohort_mutations)


# Define function for encoding variant_classification ---------------------

encode_variant_class <- function(variant_class){
  switch(as.character(variant_class),
         Silent = 0,
         Missense_Mutation = 1,
         In_Frame_Del = 2,
         In_Frame_Ins = 3,
         Splice_Site = 4,
         Nonsense_Mutation = 5,
         Del_Other = 6,
         Frame_Shift_Del = 7,
         Ins_Other = 8,
         Frame_Shift_Ins = 9,
         Nonstop_Mutation = 10)
}


# Encode variant_classification in cohort_mutations data frame ------------

cohort_mutations <- within(cohort_mutations, {
  encoded_variant_class <- vapply(Variant_Classification, encode_variant_class,
                                  FUN.VALUE = 0)
})

cohort_mutations <- droplevels(subset(cohort_mutations, 
                                      Variant_Classification != "Silent"))

genes_per_patient <- ddply(cohort_mutations, 
                           .(patient_name, Entrez_Gene_Id, Hugo_Symbol), 
                           summarize, 
                           most_deleterious = max(encoded_variant_class))

affected_patients <- ddply(genes_per_patient, 
                           .(Entrez_Gene_Id), 
                           summarize, 
                           n_affected_patients = length(patient_name))

affected_patients <- within(affected_patients, 
               Entrez_Gene_Id <- reorder(Entrez_Gene_Id, 
                                         desc(n_affected_patients)))

plotted_genes <- droplevels(arrange(affected_patients, desc(n_affected_patients), 
                         Entrez_Gene_Id)[1:n_plotted_genes,][,1])

plotted_genes_per_patient <- genes_per_patient[genes_per_patient$Entrez_Gene_Id 
                                              %in% plotted_genes,]

ggplot_data <- within(plotted_genes_per_patient, {
  gene_name = ifelse(as.character(Hugo_Symbol) == "", 
                     as.character(Entrez_Gene_Id), 
                     as.character(Hugo_Symbol))
})

my_palette <- brewer.pal(n = 10, "Spectral")

(p <- ggplot(ggplot_data, aes(patient_name, gene_name)) + 
   geom_tile(aes(fill = as.factor(most_deleterious)), colour = "white") + 
   scale_fill_brewer(palette = "Spectral")
)

base_size <- 14

(p + theme_grey(base_size = base_size) + 
   labs(x = "Patient", y = "Gene") +
   scale_x_discrete(expand = c(0, 0)) +
   scale_y_discrete(expand = c(0, 0)) + 
   scale_fill_discrete(name = "Mutation Type")
)







