# gene_patient_heatmap.R
# 2014-08-11
# Rebecca Lea Johnston and Bruno Grande
# Usage:
#     R --no-save --args ~/Desktop/gene_patient_heatmap/patients.txt
# Arguments:
#     1) File containing two tab-separated columns: patient name and file path
#        to MAF file


# Loading Libraries -------------------------------------------------------

library(ggplot2)
library(plyr)


# Argument Parsing --------------------------------------------------------

# args <- commandArgs(TRUE)

# cohort_file <- args[1]
# selection <- args[2]

cohort_file <- file.path(
    "/Users/bgrande/Desktop/gene_patient_heatmap/patients.txt"
)

selection <- as.character(50)

if (file.exists(selection)) {  
    # The specified argument is a file, so expecting a list of 
    # Ensembl gene IDs
    interesting_genes <- readLines(selection, n = -1)
} else if (is.integer(as.integer(selection))) {
    # The specified argument is an integer (N), so displaying 
    # the top N most mutated genes
    n_plotted_genes <- as.integer(selection)
} else {
    stop(paste("Specified argument is not a file with a list of Ensembl gene ",
               "IDs, nor an integer for the number of top ranking genes to ",
               "display."))
}


# Building the cohort_mutations data frame --------------------------------

cohort <- read.delim(
    cohort_file, 
    header = FALSE, 
    col.names = c("patient_name", "path") ,
    comment.char = "#"
)

cohort_mutations <- data.frame()

for (i in 1:nrow(cohort)) {
    patient_name <- cohort[i, "patient_name"]
    patient_maf <- read.delim(file.path(cohort$path[i]))
    wanted_columns <- patient_maf[
        , c("Hugo_Symbol", "Entrez_Gene_Id", "Variant_Classification")
    ]
    wanted_columns["patient_name"] <- patient_name
    cohort_mutations <- rbind(cohort_mutations, wanted_columns)
}


# Define function for encoding variant_classification ---------------------

encode_variant_class <- function(variant_class){
    switch(
        as.character(variant_class),
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
        Nonstop_Mutation = 10
    )
}


# Restrict to list of interesting genes, if specified ---------------------

if (exists("interesting_genes") > 0) {
    cohort_mutations <- cohort_mutations[
        cohort_mutations$Entrez_Gene_Id %in% interesting_genes, 
    ]
}


# Encode variant_classification in cohort_mutations data frame ------------

cohort_mutations <- within(
    cohort_mutations, 
    {encoded_variant_class <- vapply(
        Variant_Classification, 
        encode_variant_class, 
        FUN.VALUE = 0
    )}
)


# Remove silent variants --------------------------------------------------

cohort_mutations <- droplevels(
    subset(cohort_mutations, Variant_Classification != "Silent")
)


# Select the most deleterious mutation for each patient and gene ----------

genes_per_patient <- ddply(
    cohort_mutations, 
    .(patient_name, Entrez_Gene_Id, Hugo_Symbol, Variant_Classification), 
    summarize, 
    most_deleterious = max(encoded_variant_class)
)


# Create the ranked list of genes based on the number of patients  --------

affected_patients <- ddply(
    genes_per_patient, 
    .(Entrez_Gene_Id), 
    summarize, 
    n_affected_patients = length(patient_name)
)

affected_patients <- within(
    affected_patients, 
    Entrez_Gene_Id <- reorder(Entrez_Gene_Id, desc(n_affected_patients))
)


# Selecting the top N genes based, if specified -------------------

if (exists("n_plotted_genes")) {
    plotted_genes <- droplevels(
        arrange(
            affected_patients, 
            desc(n_affected_patients), 
            Entrez_Gene_Id
        )[1:n_plotted_genes,][,1]
    )
}

plotted_genes_per_patient <- genes_per_patient[
    genes_per_patient$Entrez_Gene_Id %in% plotted_genes, 
]

ggplot_data <- within(plotted_genes_per_patient, {
  gene_name = ifelse(as.character(Hugo_Symbol) == "", 
                     as.character(Entrez_Gene_Id), 
                     as.character(Hugo_Symbol))
})

ggplot_data$gene_name <- factor(ggplot_data$gene_name)

ggplot_data <- within(ggplot_data, 
    Variant_Classification <- reorder(Variant_Classification, -most_deleterious)
)

ggplot_data <- merge(ggplot_data, affected_patients)

ggplot_data <- within(ggplot_data, 
                      gene_name <- reorder(gene_name, n_affected_patients)
)

affected_genes <- ddply(ggplot_data, .(patient_name), summarize, n_affected_genes = length(gene_name))

ggplot_data <- merge(ggplot_data, affected_genes)

ggplot_data <- within(ggplot_data, 
                      patient_name <- reorder(patient_name, -n_affected_genes)
)


# Creating NA matrix for filling missing values in geom_tile plot  --------

all_combinations = expand.grid(patient_name = unique(ggplot_data$patient_name), 
                               gene_name = unique(ggplot_data$gene_name))

ggplot_data_fixed <- merge(ggplot_data, all_combinations, 
                           by = c("patient_name", "gene_name"),
                           all = TRUE)


# Plotting ----------------------------------------------------------------

(heatmap_plot <- ggplot(ggplot_data_fixed, aes(patient_name, gene_name)) +
    geom_tile(aes(fill = Variant_Classification, height = 1, width = 1), 
                hjust = 0, vjust = 1, hpad = 0, vpad = 0, colour = "#818181") +
    theme_grey(base_size = 12) +
    coord_fixed(ratio = 1) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Patient", y = "Gene") +
    theme(axis.ticks = element_blank(),
          axis.title = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(vjust = -1),
          axis.title.y = element_text(vjust = 1),
          axis.text.x = element_text(hjust = 0, vjust = 0.5, angle = 270),
          axis.text.y = element_text(vjust = 0.5)) +
    scale_fill_brewer(type = "seq", palette = "Spectral", drop =TRUE, 
                      na.value = "#E5E5E5")
)

# Rough calculation for the size of the output file
calculated_width = 0.173 * nlevels(ggplot_data_fixed$patient_name) + 3.7
calculated_height = 0.173 * nlevels(ggplot_data_fixed$gene_name) + 2

ggsave(file="~/Desktop/gene_patient_heatmap/figure.pdf",
       width = calculated_width, height = calculated_height)




















