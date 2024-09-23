# Title: Analysis of bacterial genomes evolved in the absence of phage
# Author: Reena Debray
# Date: Feb 1, 2022
# Modified: Amairani Cancino Bello (04/Sep/2024)


### Functions
## This function uses the output from the program breseq to construct a matrix of genes that were mutated in each population during experimental evolution
### It takes data in the following form: "breseq", a data frame with an entry for each mutation in each population at each time and its corresponding allele frequency
### "gene_list", a list of all genes observed in the population or subpopulation of interest
### "sample_list", a list of all samples in the population or subpopulation of interest 
### AF_min, a minimum allele frequency for consideration. I set AF_min=0 (no minimum), though note that breseq only returns polymorphisms with a population frequency of 0.05 or greater.



# --- packages ----
#' @import readxl : This package is used to read Excel files (.xls and .xlsx).
#' @import ggplot2: This is one of the most popular packages in R for data visualization. You can create basic plots (such as histograms and bar charts) as well as more complex visualizations (like heatmaps and violin plots).
#' @import viridis : This package provides color palettes for plots in R. The viridis palettes are particularly useful in visualizations like heatmaps or scatter plots where color differentiation is important.

install.packages("viridis")
library(readxl)
library(ggplot2)
library(viridis)


#--- Specify the full directory path ----
#' @title Specify the Directory Path
#' @description Sets the directory path containing the mutation tables.
#' @param breseq A string specifying the directory path.
#' @usage breseq <- "C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables"
#Directory containing the mutation tables.
breseq <- "C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables"

#--- Generate Annotared Gene Matrix---- 
#' @title Generate Annotated Gene Matrix
#' @description Creates a matrix annotated with genes and samples using the breseq mutation data.
#' @param breseq Data frame containing breseq mutation data.
#' @param gene_list List of genes to be included in the matrix.
#' @param sample_list List of samples to be included in the matrix.
#' @param AF_min Minimum allele frequency for filtering.
#' @return A matrix annotated with the number of mutation sites per gene per sample.
#' @export
gene_matrix_annotated<-function(breseq,gene_list,sample_list,AF_min){
  LOR_gene_matrix<-matrix(ncol=length(gene_list),nrow=length(sample_list))
  colnames(LOR_gene_matrix)=gene_list
  rownames(LOR_gene_matrix)=sample_list
  
  for (gene in gene_list){
    for (sample in sample_list){
      if (gene%in%breseq[breseq$Sample==sample & breseq$freq>=AF_min,"gene_2"]){num_sites<-length(unique(breseq[breseq$Sample==sample & breseq$gene_2==gene & !is.na(breseq$gene_2),"position"])); LOR_gene_matrix[sample,gene]=num_sites}
      else{LOR_gene_matrix[sample,gene]=0}
    }
  }
  return(LOR_gene_matrix)
}


#---Read and filter breseq output----
#' @title Read and Filter breseq Output
#' @description Reads the breseq mutation data files, formats the data frame, and filters the mutation data based on various criteria.
#' @return A data frame containing filtered mutation data.
#' @export

breseq_annotated<-data.frame()
filenames=list.files("C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables")
for (file in filenames){
  ann <- data.frame(read_excel(paste("C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Mutation tables", file, sep="/")))
  ann$Line<-unlist(strsplit(file,split="_"))[1]
  ann$Passage<-unlist(strsplit(file,split="_"))[2]
  ann$Sample<-paste(ann$Line,ann$Passage,sep="_")
  if (ncol(ann)>10){ann<-ann[,-c(8:12)]} # some files have an extra column of notes for annotations without evidence (*) resulting in blanks for cols 8-12
  breseq_annotated<-rbind(breseq_annotated,ann)
}
dim(breseq_annotated)

#--- Remove unassigned evidence----
#' @title Remove Unassigned Evidence
#' @description Filters out unassigned evidence from the breseq mutation data.
#' @param breseq_annotated Data frame containing mutation data.
#' @return Filtered data frame without unassigned evidence.
#' @export

breseq_annotated<-breseq_annotated[!breseq_annotated$evidence%in%c(NA,"Unassigned missing coverage evidence","*","Unassigned new junction evidence","?"),]
dim(breseq_annotated)

# --- Function: Remove Sites Differing from Reference ----
#' @title Remove Sites Differing from Reference
#' @description Removes any mutation sites that differ from the reference in the ANCDC3000 line.
#' @param breseq_annotated Data frame. Contains the annotated breseq data.
#' @return Data frame without sites differing from the reference in ANCDC3000.

## Remove any sites that differ from reference in ANCDC3000
ANCDC3000_sites<-breseq_annotated[breseq_annotated$Line=="ANCDC3000","position"]
breseq_annotated<-breseq_annotated[!breseq_annotated$position%in%ANCDC3000_sites,]
dim(breseq_annotated)


# --- Function: Remove Fixed Sites ----
#' @title Remove Fixed Sites
#' @description Removes sites that are fixed in every single line.
#' @param breseq_annotated Data frame. Contains the annotated breseq data.
#' @return Data frame without fixed mutation sites.

## Remove sites that are fixed in every single line
tab<-data.frame(table(breseq_annotated[breseq_annotated$freq==1,"Line"],breseq_annotated[breseq_annotated$freq==1,"position"])) #table of fixed sites
tab<-tab[tab$Freq>0,] # don't consider line-mutation combinations that never occurred
mut_tab<-data.frame(table(tab$Var2)) # total up number of lines in which mutation ever occurred
sites_to_remove<-as.character(mut_tab[mut_tab$Freq==length(unique(breseq_annotated$Line)),"Var1"])
breseq_annotated<-breseq_annotated[!breseq_annotated$position%in%sites_to_remove,]
dim(breseq_annotated)

# --- Function: Reformat Frameshift Mutations ----
#' @title Reformat Frameshift Mutations
#' @description Reformats frameshift mutations to be consistent with point mutations.
#' @param breseq_annotated Data frame. Contains the annotated breseq data.
#' @return Data frame with reformatted frameshift mutations.

## Reformat frameshift mutations to be consistent with point mutations
for (i in seq(1,nrow(breseq_annotated))){
  pos<-breseq_annotated[i,"position"]
  # find mutations with a ":"
  if (grepl(":",pos,fixed=TRUE)){
    pos<-unlist(strsplit(pos,split=":"))[1]
    pos<-paste(unlist(strsplit(pos,split=",")),collapse="")
    breseq_annotated[i,"position"]<-pos
  }
}
breseq_annotated$position<-as.numeric(breseq_annotated$position)


# --- Function: Reformat Gene Column ----
#' @title Reformat Gene Column
#' @description Reformats the gene column to remove direction.
#' @param breseq_annotated Data frame. Contains the annotated breseq data.
#' @return Data frame with reformatted gene column.

## Reformat gene column to remove direction
for (i in seq(1,nrow(breseq_annotated))){
  breseq_annotated[i,"gene_2"]<-substr(breseq_annotated[i,"gene"],1,nchar(breseq_annotated[i,"gene"])-2)
}

# --- Function: Filter Mutations with Sliding Window ----
#' @title Filter Mutations with Sliding Window
#' @description Filters mutations using a 50 bp sliding window within each line/timepoint.
#' @param breseq_annotated Data frame. Contains the annotated breseq data.
#' @param N Numeric. Maximum number of neighbors allowed within the sliding window.
#' @return Data frame with filtered mutations.

# Scan a 50 bp sliding window within each line/timepoint. If a mutation has more than N neighbors (including repeat calls at the same position), remove it.
breseq_annotated_filtered<-data.frame()
N=4

for (sample in unique(breseq_annotated$Sample)){
  pre_filter_POS<-sort(as.numeric(breseq_annotated[breseq_annotated$Sample==sample,"position"]))
  sites_to_remove<-c()
  for (pos in pre_filter_POS){
    # sliding window
    min=pos-50
    while (min<pos){
      max<-min+50
      if (length(pre_filter_POS[pre_filter_POS>=min & pre_filter_POS<max])>N)   {sites_to_remove<-unique(c(sites_to_remove,pre_filter_POS[pre_filter_POS>=min & pre_filter_POS<max]))} # if a window contains more than 3 mutations, remove all of them
      min<-min+1
    }
  }
 
  
  #' @title Filter and Annotate Genetic Data
  #'
  #' @Description This function filters out sites based on a provided sample and a list of sites to remove, and then categorizes the annotations in the filtered data.
  #'
  #' @param breseq_annotated A data frame containing the annotated genetic data.
  #' @param sample A string indicating the sample to filter by.
  #' @param sites_to_remove A vector of positions to remove from the data.
  #'
  #' @return A data frame with filtered annotations and updated types.  
  
  ### filter sites
  post_filter<-breseq_annotated[breseq_annotated$Sample==sample & !breseq_annotated$position%in%sites_to_remove,]
  breseq_annotated_filtered<-rbind(breseq_annotated_filtered,post_filter)
}
# Print the dimensions of the filtered data
dim(breseq_annotated_filtered)

## Summarize annotations as NS, S, pseudogene, or intergenic
for (i in seq(1,nrow(breseq_annotated_filtered))){
  anno<-substr(breseq_annotated_filtered[i,"annotation"],1,5)
  if (anno=="pseud"){breseq_annotated_filtered[i,"Type"]="Pseudogene"}
  else if (anno=="inter"){breseq_annotated_filtered[i,"Type"]="Intergenic"}
  else {
    AA1<-substr(anno,1,1)
    AA2<-substr(anno,5,5)
    if (AA1==AA2){breseq_annotated_filtered[i,"Type"]="S"}
    else {breseq_annotated_filtered[i,"Type"]="NS"}
  }
} 
  
## Count number of fixed sites and polymorphisms in each population at passage 12

#' @title Count Fixed Sites and Polymorphisms in Each Population at Passage 12
#'
#' @description This script calculates the number of fixed mutations and polymorphisms in each population at passage 12, excluding resistance mutations acquired before experimental evolution. It also performs statistical analyses to compare the counts between different populations and identifies mutator populations.
#'
#' @details
#' The script includes the following analyses:
#' \itemize{
#'   \item Counts the number of fixed mutations (both nonsynonymous and synonymous) in each population at passage 12, excluding resistance mutations.
#'   \item Computes the mean and standard deviation of the fixed mutations count and performs a t-test to compare the number of fixed mutations between different populations.
#'   \item Counts the number of polymorphisms (both nonsynonymous and synonymous) in each population at passage 12, excluding resistance mutations.
#'   \item Computes the mean and standard deviation of the polymorphism count and performs a t-test to compare the number of polymorphisms between different populations.
#'   \item Identifies mutator populations based on specific criteria and performs a t-test to compare the number of fixed mutations between mutator and non-mutator populations.
#' }
#'
#' @importFrom readxl read_excel
#' @param Costs_of_res A dataframe containing positions of resistance mutations to exclude from the count.
#' @param breseq_annotated_filtered A dataframe containing mutation data, including passage, type, position, frequency, and line information.
#' @return Prints the mean, standard deviation, and results of t-tests for the number of fixed mutations and polymorphisms in each population, and provides comparisons between mutator and non-mutator populations.
#' @examples
#' \dontrun{
#' # Load necessary data and run the script to calculate counts and perform analyses
#' }
#' @export

### Fixed mutations (all types)
### Exclude resistance mutations from count (because they were acquired before experimental evolution)
Costs_of_res <- read_excel("C:/Users/amair/Desktop/R ladies/Reprohack/Sesión 2/loss-of-resistance-main/loss-of-resistance-main/Costs_of_res.xlsx")
res_mutations<-unique(Costs_of_res$Position)
a<-aggregate(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"freq"],by=list(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"Line"]),FUN=function(x){length(x[x==1])})
mean(a$x)
sd(a$x)
t.test(x~(substr(a$Group.1,1,3)=="ANC"),a,var.equal=T)

### Polymorphisms (all types)
p<-aggregate(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"freq"],by=list(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"Line"]),FUN=function(x){length(x[x<1])})
mean(p$x)
sd(p$x) 
t.test(x~(substr(p$Group.1,1,3)=="ANC"),p,var.equal=T)

### Fixed mutations in mutator populations (within genes)
mut<-c("FMS13","MS1","MS10","MS15","QAC4")
p[p$Group.1%in%mut,"Mutator"]<-"Y"
p[!p$Group.1%in%mut,"Mutator"]<-"N"
t.test(x~Mutator,p,alternative="less",var.equal=T)


#--- Figure 4: Genomic parallelism across populations---- 

#' @title Genomic Parallelism Across Populations: Sørensen-Dice Similarity Analysis
#'
#' @description This script analyzes genomic parallelism across populations initiated with phage-resistant bacteria. It identifies non-synonymous mutations that are not fixed at passage 0, constructs a gene matrix, calculates Sørensen-Dice similarity coefficients for all pairs of samples, and annotates resistance genes of the populations.
#'
#' @details
#' The script performs the following steps:
#' \itemize{
#'   \item Identifies populations initiated with phage-resistant bacteria, excluding ancestral controls.
#'   \item Constructs a data frame of non-synonymous mutations that are not fixed at passage 0 across these populations.
#'   \item Creates a gene matrix using only mutations not fixed at passage 0.
#'   \item Calculates Sørensen-Dice similarity coefficients for each pair of samples to assess the extent of overlap in gene mutations between samples.
#'   \item Annotates the resistance genes of the populations and determines whether pairs of samples have the same or different resistance genes.
#' }
#'
#' @param breseq_annotated_filtered A dataframe containing annotated mutation data, including line, passage, position, type, frequency, gene, and sample information.
#' @param Costs_of_res A dataframe containing information about resistance genes associated with different populations.
#' @return A dataframe (`sim_coef_NS_noP0`) containing Sørensen-Dice similarity coefficients for each pair of samples, annotated with resistance gene information and whether the samples share the same or different resistance genes.
#' @examples
#' \dontrun{
#' # Run the script to compute similarity coefficients and annotate resistance genes across populations
#' }
#' @export


### Generate names of populations initiated with phage-resistant bacteria only (exclude ancestral controls)
n<-unique(breseq_annotated_filtered[substr(breseq_annotated_filtered$Line,1,3)!="ANC","Line"])

### Form a data frame of non-synonymous mutations not fixed at Psg 0 only
breseq_NS_noP0<-data.frame(matrix(nrow=0,ncol=ncol(breseq_annotated_filtered)))
for (line in n){
  P0_mutations<-breseq_annotated_filtered[breseq_annotated_filtered$Line==line & breseq_annotated_filtered$Passage=="P0" & breseq_annotated_filtered$freq==1,"position"]
  ###rows corresponding to P2/P12 mutations not fixed in P0 of this sample
  sample_data<-breseq_annotated_filtered[breseq_annotated_filtered$Line==line & breseq_annotated_filtered$Passage!="P0" & !breseq_annotated_filtered$position%in%P0_mutations & breseq_annotated_filtered$Type=="NS",]
  
  breseq_NS_noP0<-rbind(breseq_NS_noP0,sample_data)
}

### Make a gene matrix using only mutations not fixed in P0
gene_list<-sort(unique(breseq_NS_noP0$"gene_2"))
sample_list<-unique(breseq_NS_noP0$Sample)
LOR_gene_matrix_NS_noP0<-gene_matrix_annotated(breseq_NS_noP0,gene_list,sample_list,0)

### Calculate Sørensen-Dice similarity coefficients for every pair of samples
all_pairs<-combn(sample_list,2)
sim_coef_NS_noP0<-data.frame(matrix(nrow=0,ncol=3))

for (i in seq(1,ncol(all_pairs))){
  sample1<-all_pairs[1,i]
  sample2<-all_pairs[2,i]
  
  ### Names of genes with any mutations from ancestral genotype (any number of independent mutations per gene)
  sample1_genes<-names(LOR_gene_matrix_NS_noP0[sample1,apply(data.frame(LOR_gene_matrix_NS_noP0[sample1,]),1,function(x){x>0})])
  sample2_genes<-names(LOR_gene_matrix_NS_noP0[sample2,apply(data.frame(LOR_gene_matrix_NS_noP0[sample2,]),1,function(x){x>0})])
  
  ### Calculate the extent of overlap
  sim_coef<-(2*length(intersect(sample1_genes,sample2_genes))) / (length(sample1_genes) + length(sample2_genes))
  sim_coef_NS_noP0<-rbind(sim_coef_NS_noP0,c(sample1,sample2,sim_coef))
}
colnames(sim_coef_NS_noP0)=c("sample1","sample2","sim_coef")
sim_coef_NS_noP0$sim_coef<-as.numeric(sim_coef_NS_noP0$sim_coef)

### Annotate the resistance genes of the populations
for (i in seq(1,nrow(sim_coef_NS_noP0))){
  sim_coef_NS_noP0[i,"line1"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample1"],split="_"))[1]
  sim_coef_NS_noP0[i,"passage1"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample1"],split="_"))[2]
  sim_coef_NS_noP0[i,"line2"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample2"],split="_"))[1]
  sim_coef_NS_noP0[i,"passage2"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample2"],split="_"))[2]
  
  sim_coef_NS_noP0[i,"gene1"]=as.character(unique(Costs_of_res[Costs_of_res$Population==sim_coef_NS_noP0[i,"line1"],"Gene"]))
  sim_coef_NS_noP0[i,"gene2"]=as.character(unique(Costs_of_res[Costs_of_res$Population==sim_coef_NS_noP0[i,"line2"],"Gene"]))
}

# Annotate whether pairs have the same or different resistance genes
sim_coef_NS_noP0[!is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2) & sim_coef_NS_noP0$gene1==sim_coef_NS_noP0$gene2,"same_gene"]="same"
sim_coef_NS_noP0[!is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2) & sim_coef_NS_noP0$gene1!=sim_coef_NS_noP0$gene2,"same_gene"]="different"

##--- Randomization test (Day 6)---- 
#' @title Randomization Test and Visualization for Genomic Parallelism Analysis (Day 6)
#'
#' @description This script performs a randomization test to assess the significance of the overlap in acquired mutations between populations, based on the Sørensen-Dice similarity coefficients. It constructs a test statistic, performs randomization and resampling, calculates the p-value, and visualizes the results using a boxplot.
#'
#' @details 
#' The script performs the following steps:
#' \itemize{
#'   \item Constructs a test statistic representing the difference in mean similarity coefficients between pairs of samples with the same resistance gene and those with different resistance genes.
#'   \item Conducts a randomization test by permuting the labels ("same" or "different" resistance gene) and recalculating the test statistic 10,000 times to generate a distribution of the test statistic under the null hypothesis.
#'   \item Calculates the p-value as the proportion of permuted test statistics that are more extreme than the observed test statistic.
#'   \item Prepares the data for visualization and creates a boxplot to illustrate the overlap in acquired mutations between samples with the same or different resistance genes.
#'   \item Annotates the plot with the calculated p-value.
#' }
#'
#' @param sim_coef_NS_noP0 A dataframe containing Sørensen-Dice similarity coefficients for each pair of samples, annotated with resistance gene information.
#' @return A plot displaying the Sørensen-Dice similarity coefficients for pairs of samples with the same or different resistance genes, and the p-value for the randomization test.
#' @examples
#' \dontrun{
#' # Run the script to perform the randomization test and generate the plot
#' }
#' @export
  

### Construct test statistic
tmp<-sim_coef_NS_noP0[sim_coef_NS_noP0$passage1=="P2" & sim_coef_NS_noP0$passage2=="P2" & !is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2),]
test_same_P2<-tmp[tmp$same_gene=="same","sim_coef"]
test_diff_P2<-tmp[tmp$same_gene=="different","sim_coef"]
test_gene_P2<-mean(test_same_P2)-mean(test_diff_P2)

### Randomization & resampling
perm_gene_P2<-c()
set.seed(123)
for (i in seq(1,10000)){
  tmp$same_gene<-sample(tmp$same_gene,replace=F)
  perm_same_P2<-mean(tmp[tmp$same_gene=="same","sim_coef"])
  perm_diff_P2<-mean(tmp[tmp$same_gene=="different","sim_coef"])
  perm_gene_P2<-c(perm_gene_P2,perm_same_P2-perm_diff_P2)
}

### p-value (proportion of permutations more extreme than observed value)
print(length(perm_gene_P2[perm_gene_P2>=test_gene_P2])/length(perm_gene_P2))

#--- Modified to add value p----

# Prepare data
test_gene_df <- data.frame(
  sim_coef = c(test_same_P2, test_diff_P2),
  label = c(rep("Same\nresistance\ngene", length(test_same_P2)), rep("Different\nresistance\ngene", length(test_diff_P2)))
)
colnames(test_gene_df) <- c("sim_coef", "label")
test_gene_df$label <- factor(test_gene_df$label, levels = c("Same\nresistance\ngene", "Different\nresistance\ngene"))

# Calculate day 6 p-value
p_value_P2 <- length(perm_gene_P2[perm_gene_P2 >= test_gene_P2]) / length(perm_gene_P2)
print(p_value_P2)

# Set up p-value annotationn
if (p_value_P2 < 0.001) {
  annotation <- "p < 0.001***"
} else {
  annotation <- paste0("p = ", format(p_value_P2, digits = 3))
}

#--- Fig. 4A ---- 
# Creation and display of the graph
p <- ggplot(test_gene_df, aes(label, sim_coef)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.8) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.2) +
  theme_classic(base_size = 18) +
  xlab("") +
  ylab("Overlap in acquired mutations\n(Pairwise Sørensen-Dice similarity)") +
  theme(panel.spacing = unit(3, "lines"),
        strip.text = element_text(size = 18, face = "bold")) +
  ylim(0.1, 0.55) +
  annotate("text", x = 1.5, y = 0.55, label = annotation, size = 6, hjust = 0.5, vjust = 1.5, color = "black")

print(p)

#--- Fig.4B ---- 
#' @title Fig. 4B: Heatmap of Mutation Percentages by Resistance Gene
#'
#' @description This script generates a heatmap displaying the percentage of populations with each resistance gene that have mutations in specific top genes by Day 6 of experimental evolution. The heatmap visualizes the proportions of populations with mutations in the top 20 genes, categorized by their resistance gene.
#'
#' @details 
#' The script performs the following steps:
#' \itemize{
#'   \item Identifies the top 20 genes mutated across populations at Day 6 of experimental evolution.
#'   \item Constructs a data frame of mutations that appeared by Day 6, excluding mutations fixed before the start of the experiment and phage resistance mutations.
#'   \item Annotates each population according to its resistance gene: "rfbA", "glycosyl", or "glycoside".
#'   \item Calculates the percentage of populations with each resistance gene that have mutations in each of the top genes.
#'   \item Prepares and plots a heatmap showing the percentage of populations with mutations in the top genes, using a color scale to represent the percentage of populations.
#' }
#'
#' @param breseq_NS_noP0 A data frame containing information on non-synonymous mutations not fixed at passage 0, including gene descriptions, lines, and mutation counts.
#' @param Costs_of_res A data frame containing information on resistance genes associated with specific populations.
#' @return A heatmap displaying the percentage of populations with each resistance gene that have mutations in the top genes by Day 6.
#' @examples
#' \dontrun{
#' # Run the script to generate the heatmap for the top 20 genes
#' }
#' @export

## Identify which genes were mutated in which populations
top_genes_in_study<-data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P2","description"],breseq_NS_noP0[breseq_NS_noP0$Passage=="P2","Line"]))
colnames(top_genes_in_study)=c("description","Line","hits")
top_genes<-names(tail(sort(table(top_genes_in_study[top_genes_in_study$hits>0,"description"])),20))

### Construct a data frame of mutations that appeared by Day 6, excluding phage resistance mutations or any other mutations fixed before the start of the experiment
breseq_gene_analysis<-data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P2","description"],breseq_NS_noP0[breseq_NS_noP0$Passage=="P2","Line"]))
colnames(breseq_gene_analysis)=c("gene","Line","hits")

# Define the vectors or lists of values corresponding to each resistance gene
rfbA <- c("FMS6", "VCM4", "MS2", "FMS4", "MS10", "QAC4")  
glycosyl <- c("VCM19", "MS12", "FMS13", "SNK12", "QAC3", "SNK11", "VCM17", "SNK6", "MS1")  
glycoside <- c("SNK7", "FMS11", "FMS12", "FMS9", "FMS16")  

### Annotate each population by its resistance gene
breseq_gene_analysis[breseq_gene_analysis$Line %in% rfbA, "res_gene"] <- "rfbA"
breseq_gene_analysis[breseq_gene_analysis$Line %in% glycosyl, "res_gene"] <- "glycosyl"
breseq_gene_analysis[breseq_gene_analysis$Line %in% glycoside, "res_gene"] <- "glycoside"

# Remove rows where res_gene is NA
breseq_gene_analysis <- breseq_gene_analysis[!is.na(breseq_gene_analysis$res_gene), ]

### Calculate the percentage of populations with each resistance gene that mutated in each of the other genes
breseq_gene_props <- data.frame(matrix(nrow = 0, ncol = 3))
l <- list(rfbA, glycosyl, glycoside)
names(l) <- c("rfbA mutants", "PSPTO_4988 mutants", "PSPTO_4991 mutants")

# Loop through genes and calculate mutation proportions
for (gene in top_genes) {
  for (j in seq_along(l)) {
    res_gene <- l[[j]]
    prop <- nrow(breseq_gene_analysis[breseq_gene_analysis$Line %in% res_gene & 
                                        breseq_gene_analysis$gene == gene & 
                                        breseq_gene_analysis$hits > 0, ]) / 
      nrow(breseq_gene_analysis[breseq_gene_analysis$Line %in% res_gene & 
                                  breseq_gene_analysis$gene == gene, ])
    breseq_gene_props <- rbind(breseq_gene_props, c(gene, names(l)[j], prop))
  }
}

# Sort the 'res_gene' column so that "rfbA mutants" appears first, followed by "PSPTO_4988 mutants" and then "PSPTO_4991 mutants"
colnames(breseq_gene_props)
colnames(breseq_gene_props) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props$res_gene <- factor(breseq_gene_props$res_gene, levels = c("PSPTO_4991 mutants","PSPTO_4988 mutants","rfbA mutants"))
breseq_gene_props$percent_pops <- as.numeric(breseq_gene_props$percent_pops)
breseq_gene_props[breseq_gene_props$percent_pops == 0, "percent_pops"] <- 10^-4

# Replace problematic hyphens and ensure proper encoding
breseq_gene_props$top_gene <- gsub("‑", "-", breseq_gene_props$top_gene)  

# Maintain the order of factors on the X-axis
breseq_gene_props$top_gene <- factor(breseq_gene_props$top_gene, levels = rev(unique(breseq_gene_props$top_gene)))

#--- Plot Fig 4B---- 
ggplot(breseq_gene_props, aes(top_gene, res_gene, fill = percent_pops * 100)) +
  geom_tile(color = "grey80") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, family = "Arial")  # Define all adjustments in one call
  ) +
  scale_fill_viridis(
    name = "Populations\nwith mutation (%)",
    limits = c(0, 100),  # Define scale limits
    breaks = c(0, 25, 50, 75, 100)  # Define specific breaks to display
  ) +
  ylab("") +
  xlab("") +
  ggtitle("Day 6 of experimental evolution")


##--- Randomization Test (Day 36)----
#' @title Randomization Test and Visualization for Genomic Parallelism Analysis (Day 36)
#'
#' @description This script performs a randomization test to assess the significance of the overlap in acquired mutations between populations, based on the Sørensen-Dice similarity coefficients. It constructs a test statistic, performs randomization and resampling, calculates the p-value, and visualizes the results using a boxplot.
#'
#' @details 
#' The script performs the following steps:
#' \itemize{
#'   \item Constructs a test statistic representing the difference in mean Sørensen-Dice similarity coefficients between pairs of samples with the same resistance gene and those with different resistance genes for Day 36 (P12).
#'   \item Conducts a randomization test by permuting the labels ("same" or "different" resistance gene) and recalculating the test statistic 10,000 times to generate a distribution of the test statistic under the null hypothesis.
#'   \item Calculates the p-value as the proportion of permuted test statistics that are more extreme than the observed test statistic.
#'   \item Prepares the data for visualization and creates a boxplot to illustrate the overlap in acquired mutations between samples with the same or different resistance genes.
#' }
#'
#' @param sim_coef_NS_noP0 A dataframe containing Sørensen-Dice similarity coefficients for each pair of samples, annotated with passage and resistance gene information.
#' @return A ggplot object displaying the Sørensen-Dice similarity coefficients for pairs of samples with the same or different resistance genes, without p-value annotation.
#' @examples
#' \dontrun{
#' # Run the script to perform the randomization test and generate the plot
#' randomization_test_day36()
#' }
#' @import ggplot2
#' @importFrom stats sd
#' @export
#' 

# Construction of the Test Statistic
  tmp2 <- sim_coef_NS_noP0[sim_coef_NS_noP0$passage1 == "P12" & 
                             sim_coef_NS_noP0$passage2 == "P12" & 
                             !is.na(sim_coef_NS_noP0$gene1) & 
                             !is.na(sim_coef_NS_noP0$gene2), ]
  
  test_same_P12 <- tmp2[tmp2$same_gene == "same", "sim_coef"]
  test_diff_P12 <- tmp2[tmp2$same_gene == "different", "sim_coef"]
  test_gene_P12 <- mean(test_same_P12) - mean(test_diff_P12)
  
  # Randomization and Resampling
  perm_gene_P12 <- c()
  set.seed(123)
  for (i in seq(1, 10000)) {
    tmp2$same_gene <- sample(tmp2$same_gene, replace = FALSE)
    perm_same_P12 <- mean(tmp2[tmp2$same_gene == "same", "sim_coef"])
    perm_diff_P12 <- mean(tmp2[tmp2$same_gene == "different", "sim_coef"])
    perm_gene_P12 <- c(perm_gene_P12, perm_same_P12 - perm_diff_P12)
  }
  
  # P-value (proportion of permutations more extreme than the observed value)
  p_value_P12 <- length(perm_gene_P12[perm_gene_P12 >= test_gene_P12]) / length(perm_gene_P12)
  print(p_value_P12)
  
  # Preparing Data for Visualization
  test_gene_df <- data.frame(
    sim_coef = c(test_same_P12, test_diff_P12),
    label = c(rep("Same\nresistance\ngene", length(test_same_P12)), rep("Different\nresistance\ngene", length(test_diff_P12)))
  )
  colnames(test_gene_df) <- c("sim_coef", "label")
  test_gene_df$label <- factor(test_gene_df$label, levels = c("Same\nresistance\ngene", "Different\nresistance\ngene"))
  
  #--- Fig 5A ---- 
  #Creating and visualizing the plot without p-value annotation
  p <- ggplot(test_gene_df, aes(label, sim_coef)) +
    geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.8) +
    geom_jitter(width = 0.05, size = 2, alpha = 0.2) +
    theme_classic(base_size = 18) +
    xlab("") +
    ylab("Overlap in acquired mutations\n(Paired Sørensen-Dice similarity)") +
    theme(panel.spacing = unit(3, "lines"),
          strip.text = element_text(size = 18, face = "bold")) +
    ylim(0.1, 0.55)
  
  print(p)


## Identifying which genes mutated in which populations
  #' Identifying and Analyzing Gene Mutations in Different Populations
  #'
  #' This script identifies which genes are mutated in different populations and analyzes the proportions of mutations for each gene within various resistance gene groups. It also constructs a plot to visualize the percentage of populations with mutations for each gene.
  #'
  #' @details
  #' The script performs the following steps:
  #' \itemize{
  #'   \item Constructs a data frame of top genes with mutations in the P12 passage.
  #'   \item Creates a gene analysis table for the P12 passage.
  #'   \item Defines vectors for resistance genes and annotates each population based on these vectors.
  #'   \item Calculates the percentage of populations with each resistance gene that mutated in each of the other genes.
  #'   \item Sorts and processes the data frame for plotting.
  #'   \item Generates a heatmap plot (Fig 5B) to visualize the mutation proportions across different top genes and resistance gene groups.
  #' }
  #'
  #' @return A ggplot object showing the percentage of populations with mutations for each gene across different resistance gene groups.
  #'
  #' @import ggplot2
  #' @import viridis
  #' @export
  #'
  #' @examples
  #' # Load necessary libraries
  #' library(ggplot2)
  #' library(viridis)
  #'
  #' # Assuming 'breseq_NS_noP0' is loaded in the environment
  
top_genes_in_study2 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "Line"]))
colnames(top_genes_in_study2) <- c("description", "Line", "hits")
top_genes <- names(tail(sort(table(top_genes_in_study2[top_genes_in_study2$hits > 0, "description"])), 20))

### Constructing a gene analysis table for P12
breseq_gene_analysis2 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "Line"]))
colnames(breseq_gene_analysis2) <- c("gene", "Line", "hits")

# Define the resistance gene vectors
rfbA <- c("FMS6", "VCM4", "MS2", "FMS4", "MS10", "QAC4")  
glycosyl <- c("VCM19", "MS12", "FMS13", "SNK12", "QAC3", "SNK11", "VCM17", "SNK6", "MS1")  
glycoside <- c("SNK7", "FMS11", "FMS12", "FMS9", "FMS16")  

# Annotate each population according to its resistance gene
breseq_gene_analysis2[breseq_gene_analysis2$Line %in% rfbA, "res_gene"] <- "rfbA"
breseq_gene_analysis2[breseq_gene_analysis2$Line %in% glycosyl, "res_gene"] <- "glycosyl"
breseq_gene_analysis2[breseq_gene_analysis2$Line %in% glycoside, "res_gene"] <- "glycoside"

# Remove rows where res_gene is NA
breseq_gene_analysis2 <- breseq_gene_analysis2[!is.na(breseq_gene_analysis2$res_gene), ]

### Calculate the percentage of populations with each resistance gene that mutated in each of the other genes
breseq_gene_props <- data.frame(matrix(nrow = 0, ncol = 3))
l <- list(rfbA = rfbA, glycosyl = glycosyl, glycoside = glycoside)
names(l) <- c("rfbA mutants", "PSPTO_4988 mutants", "PSPTO_4991 mutants")

# Loop through genes and calculate mutation proportions
for (gene in top_genes) {
  for (j in seq_along(l)) {
    res_gene <- l[[j]]
    prop <- nrow(breseq_gene_analysis2[breseq_gene_analysis2$Line %in% res_gene & 
                                         breseq_gene_analysis2$gene == gene & 
                                         breseq_gene_analysis2$hits > 0, ]) / 
      nrow(breseq_gene_analysis2[breseq_gene_analysis2$Line %in% res_gene & 
                                   breseq_gene_analysis2$gene == gene, ])
    breseq_gene_props <- rbind(breseq_gene_props, c(gene, names(l)[j], prop))
  }
}

# Sort the 'res_gene' column so that "rfbA mutants" appears first, followed by "PSPTO_4988 mutants" and then "PSPTO_4991 mutants"
colnames(breseq_gene_props)
colnames(breseq_gene_props) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props$res_gene <- factor(breseq_gene_props$res_gene, levels = c("PSPTO_4991 mutants", "PSPTO_4988 mutants", "rfbA mutants"))

# Assign column names
colnames(breseq_gene_props) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props$percent_pops <- as.numeric(breseq_gene_props$percent_pops)
breseq_gene_props[breseq_gene_props$percent_pops == 0, "percent_pops"] <- 10^-4

# Replace problematic dashes and ensure proper encoding
breseq_gene_props$top_gene <- gsub("‑", "-", breseq_gene_props$top_gene)

# Maintain the order of factors on the X-axis
breseq_gene_props$top_gene <- factor(breseq_gene_props$top_gene, levels = rev(unique(breseq_gene_props$top_gene)))

#--- Plot Fig 5B----
ggplot(breseq_gene_props, aes(top_gene, res_gene, fill = percent_pops * 100)) +
  geom_tile(color = "grey80") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, family = "Arial")) +
  scale_fill_viridis(name = "Populations\nwith mutation (%)", limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  ylab("") +
  xlab("") +
  ggtitle("Day 36 of experimental evolution")



sessionInfo()
#R version 4.4.0 (2024-04-24 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 22631)

#Matrix products: default


#locale:
 # [1] LC_COLLATE=Spanish_Mexico.utf8  LC_CTYPE=Spanish_Mexico.utf8    LC_MONETARY=Spanish_Mexico.utf8
#[4] LC_NUMERIC=C                    LC_TIME=Spanish_Mexico.utf8    

#time zone: America/Mexico_City
#tzcode source: internal

#attached base packages:
 # [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
 # [1] viridis_0.6.5     viridisLite_0.4.2 ggplot2_3.5.1     readxl_1.4.3     

#loaded via a namespace (and not attached):
#  [1] vctrs_0.6.5       cli_3.6.3         knitr_1.48        rlang_1.1.4       xfun_0.47        
#[6] generics_0.1.3    labeling_0.4.3    glue_1.7.0        colorspace_2.1-0  htmltools_0.5.8.1
#[11] gridExtra_2.3     scales_1.3.0      fansi_1.0.6       rmarkdown_2.28    grid_4.4.0       
#[16] cellranger_1.1.0  evaluate_0.24.0   munsell_0.5.1     tibble_3.2.1      fastmap_1.2.0    
#[21] lifecycle_1.0.4   compiler_4.4.0    dplyr_1.1.4       pkgconfig_2.0.3   rstudioapi_0.16.0
#[26] farver_2.1.2      digest_0.6.37     R6_2.5.1          tidyselect_1.2.1  utf8_1.2.4       
#[31] pillar_1.9.0      magrittr_2.0.3    withr_3.0.1       tools_4.4.0       gtable_0.3.5     
