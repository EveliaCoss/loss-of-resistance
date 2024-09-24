# Analysis of bacterial genomes evolved in the absence of phage
# Reena Debray
# Feb 1, 2022
# Modified: Amairani Cancino Bello (04/Sep/2024)


### Functions
## This function uses the output from the program breseq to construct a matrix of genes that were mutated in each population during experimental evolution
### It takes data in the following form: "breseq", a data frame with an entry for each mutation in each population at each time and its corresponding allele frequency
### "gene_list", a list of all genes observed in the population or subpopulation of interest
### "sample_list", a list of all samples in the population or subpopulation of interest 
### AF_min, a minimum allele frequency for consideration. I set AF_min=0 (no minimum), though note that breseq only returns polymorphisms with a population frequency of 0.05 or greater.



# --- packages ----
library(readxl)
library(ggplot2)
install.packages("viridis")
library(viridis)


#--- Specify the full directory path ----
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
breseq_annotated<-breseq_annotated[!breseq_annotated$evidence%in%c(NA,"Unassigned missing coverage evidence","*","Unassigned new junction evidence","?"),]
dim(breseq_annotated)

# --- Remove Sites Differing from Reference ----
## Remove any sites that differ from reference in ANCDC3000
ANCDC3000_sites<-breseq_annotated[breseq_annotated$Line=="ANCDC3000","position"]
breseq_annotated<-breseq_annotated[!breseq_annotated$position%in%ANCDC3000_sites,]
dim(breseq_annotated)


# --- Remove Fixed Sites ----
## Remove sites that are fixed in every single line
tab<-data.frame(table(breseq_annotated[breseq_annotated$freq==1,"Line"],breseq_annotated[breseq_annotated$freq==1,"position"])) #table of fixed sites
tab<-tab[tab$Freq>0,] # don't consider line-mutation combinations that never occurred
mut_tab<-data.frame(table(tab$Var2)) # total up number of lines in which mutation ever occurred
sites_to_remove<-as.character(mut_tab[mut_tab$Freq==length(unique(breseq_annotated$Line)),"Var1"])
breseq_annotated<-breseq_annotated[!breseq_annotated$position%in%sites_to_remove,]
dim(breseq_annotated)

# --- Reformat Frameshift Mutations ----
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


# --- Reformat Gene Column ----
## Reformat gene column to remove direction
for (i in seq(1,nrow(breseq_annotated))){
  breseq_annotated[i,"gene_2"]<-substr(breseq_annotated[i,"gene"],1,nchar(breseq_annotated[i,"gene"])-2)
}

# --- Function: Filter Mutations with Sliding Window ----
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
### Construct test statistic
tmp_day6 <- sim_coef_NS_noP0[sim_coef_NS_noP0$passage1 == "P2" & sim_coef_NS_noP0$passage2 == "P2" & !is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2),]
test_same_P2_day6 <- tmp_day6[tmp_day6$same_gene == "same", "sim_coef"]
test_diff_P2_day6 <- tmp_day6[tmp_day6$same_gene == "different", "sim_coef"]
test_gene_P2_day6 <- mean(test_same_P2_day6) - mean(test_diff_P2_day6)

### Randomization & resampling
perm_gene_P2_day6 <- c()
set.seed(123)
for (i in seq(1, 10000)) {
  tmp_day6$same_gene <- sample(tmp_day6$same_gene, replace = F)
  perm_same_P2 <- mean(tmp_day6[tmp_day6$same_gene == "same", "sim_coef"])
  perm_diff_P2 <- mean(tmp_day6[tmp_day6$same_gene == "different", "sim_coef"])
  perm_gene_P2_day6 <- c(perm_gene_P2_day6, perm_same_P2 - perm_diff_P2)
}

### p-value (proportion of permutations more extreme than observed value)
print(length(perm_gene_P2_day6[perm_gene_P2_day6 >= test_gene_P2_day6]) / length(perm_gene_P2_day6))

#--- Modified for added value p----

# Prepare data
test_gene_df_day6 <- data.frame(
  sim_coef = c(test_same_P2_day6, test_diff_P2_day6),
  label = c(rep("Same\nresistance\ngene", length(test_same_P2_day6)), rep("Different\nresistance\ngene", length(test_diff_P2_day6)))
)
colnames(test_gene_df_day6) <- c("sim_coef", "label")
test_gene_df_day6$label <- factor(test_gene_df_day6$label, levels = c("Same\nresistance\ngene", "Different\nresistance\ngene"))

# Calculate day 6 p-value
p_value_P2_day6 <- length(perm_gene_P2_day6[perm_gene_P2_day6 >= test_gene_P2_day6]) / length(perm_gene_P2_day6)
print(p_value_P2_day6)

# Set up p-value annotation
if (p_value_P2_day6 < 0.001) {
  annotation_day6 <- "p < 0.001***"
} else {
  annotation_day6 <- paste0("p = ", format(p_value_P2_day6, digits = 3))
}

#--- Fig. 4A ---- 
# Creation and display of the graph
p_day6 <- ggplot(test_gene_df_day6, aes(label, sim_coef)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.8) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.2) +
  theme_classic(base_size = 18) +
  xlab("") +
  ylab("Overlap in acquired mutations\n(Pairwise Sørensen-Dice similarity)") +
  theme(panel.spacing = unit(3, "lines"),
        strip.text = element_text(size = 18, face = "bold")) +
  ylim(0.1, 0.55) +
  annotate("text", x = 1.5, y = 0.55, label = annotation_day6, size = 6, hjust = 0.5, vjust = 1.5, color = "black")

print(p_day6)

#--- Fig. 4B ---- 
## Identify which genes were mutated in which populations
top_genes_in_study_day6 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage == "P2", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage == "P2", "Line"]))
colnames(top_genes_in_study_day6) <- c("description", "Line", "hits")
top_genes_day6 <- names(tail(sort(table(top_genes_in_study_day6[top_genes_in_study_day6$hits > 0, "description"])), 20))

### Construct a data frame of mutations that appeared by Day 6, excluding phage resistance mutations or any other mutations fixed before the start of the experiment
breseq_gene_analysis_day6 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage == "P2", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage == "P2", "Line"]))
colnames(breseq_gene_analysis_day6) <- c("gene", "Line", "hits")

# Define the vectors or lists of values corresponding to each resistance gene
rfbA <- c("FMS6", "VCM4", "MS2", "FMS4", "MS10", "QAC4")  
glycosyl <- c("VCM19", "MS12", "FMS13", "SNK12", "QAC3", "SNK11", "VCM17", "SNK6", "MS1")  
glycoside <- c("SNK7", "FMS11", "FMS12", "FMS9", "FMS16")  

### Annotate each population by its resistance gene
breseq_gene_analysis_day6[breseq_gene_analysis_day6$Line %in% rfbA, "res_gene"] <- "rfbA"
breseq_gene_analysis_day6[breseq_gene_analysis_day6$Line %in% glycosyl, "res_gene"] <- "glycosyl"
breseq_gene_analysis_day6[breseq_gene_analysis_day6$Line %in% glycoside, "res_gene"] <- "glycoside"

# Remove rows where res_gene is NA
breseq_gene_analysis_day6 <- breseq_gene_analysis_day6[!is.na(breseq_gene_analysis_day6$res_gene), ]

### Calculate the percentage of populations with each resistance gene that mutated in each of the other genes
breseq_gene_props_day6 <- data.frame(matrix(nrow = 0, ncol = 3))
l <- list(rfbA, glycosyl, glycoside)
names(l) <- c("rfbA mutants", "PSPTO_4988 mutants", "PSPTO_4991 mutants")

# Loop through genes and calculate mutation proportions
for (gene in top_genes_day6) {
  for (j in seq_along(l)) {
    res_gene <- l[[j]]
    prop <- nrow(breseq_gene_analysis_day6[breseq_gene_analysis_day6$Line %in% res_gene & 
                                             breseq_gene_analysis_day6$gene == gene & 
                                             breseq_gene_analysis_day6$hits > 0, ]) / 
      nrow(breseq_gene_analysis_day6[breseq_gene_analysis_day6$Line %in% res_gene & 
                                       breseq_gene_analysis_day6$gene == gene, ])
    breseq_gene_props_day6 <- rbind(breseq_gene_props_day6, c(gene, names(l)[j], prop))
  }
}

# Sort the 'res_gene' column so that "rfbA mutants" appears first, followed by "PSPTO_4988 mutants" and then "PSPTO_4991 mutants"
colnames(breseq_gene_props_day6)
colnames(breseq_gene_props_day6) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props_day6$res_gene <- factor(breseq_gene_props_day6$res_gene, levels = c("PSPTO_4991 mutants", "PSPTO_4988 mutants", "rfbA mutants"))

# Assign column names
colnames(breseq_gene_props_day6) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props_day6$percent_pops <- as.numeric(breseq_gene_props_day6$percent_pops)
breseq_gene_props_day6[breseq_gene_props_day6$percent_pops == 0, "percent_pops"] <- 10^-4

# Replace problematic hyphens and ensure proper encoding
breseq_gene_props_day6$top_gene <- gsub("‑", "-", breseq_gene_props_day6$top_gene)  

# Maintain the order of factors on the X-axis
breseq_gene_props_day6$top_gene <- factor(breseq_gene_props_day6$top_gene, levels = rev(unique(breseq_gene_props_day6$top_gene)))

#--- Plot Fig 4B---- 
ggplot(breseq_gene_props_day6, aes(top_gene, res_gene, fill = percent_pops * 100)) +
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
# Construction of the Test Statistic
tmp2_day36 <- sim_coef_NS_noP0[sim_coef_NS_noP0$passage1 == "P12" & sim_coef_NS_noP0$passage2 == "P12" & !is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2), ]
    
test_same_P12_day36 <- tmp2_day36[tmp2_day36$same_gene == "same", "sim_coef"]
test_diff_P12_day36 <- tmp2_day36[tmp2_day36$same_gene == "different", "sim_coef"]
test_gene_P12_day36 <- mean(test_same_P12_day36) - mean(test_diff_P12_day36)
    
# Randomization and Resampling
perm_gene_P12_day36 <- c()
set.seed(123)
for (i in seq(1, 10000)) {
tmp2_day36$same_gene <- sample(tmp2_day36$same_gene, replace = FALSE)
perm_same_P12_day36 <- mean(tmp2_day36[tmp2_day36$same_gene == "same", "sim_coef"])
perm_diff_P12_day36 <- mean(tmp2_day36[tmp2_day36$same_gene == "different", "sim_coef"])
perm_gene_P12_day36 <- c(perm_gene_P12_day36, perm_same_P12_day36 - perm_diff_P12_day36)
    }
    
# P-value (proportion of permutations more extreme than the observed value)
p_value_P12_day36 <- length(perm_gene_P12_day36[perm_gene_P12_day36 >= test_gene_P12_day36]) / length(perm_gene_P12_day36)
print(p_value_P12_day36)
    
# Preparing Data for Visualization
test_gene_df_day36 <- data.frame(
sim_coef = c(test_same_P12_day36, test_diff_P12_day36),
label = c(rep("Same\nresistance\ngene", length(test_same_P12_day36)), rep("Different\nresistance\ngene", length(test_diff_P12_day36)))
)
colnames(test_gene_df_day36) <- c("sim_coef", "label")
test_gene_df_day36$label <- factor(test_gene_df_day36$label, levels = c("Same\nresistance\ngene", "Different\nresistance\ngene"))
  
#--- Fig 5A ---- 
#Creating and visualizing the plot without p-value annotation
p_day36 <- ggplot(test_gene_df_day36, aes(label, sim_coef)) +
geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.8) +
geom_jitter(width = 0.05, size = 2, alpha = 0.2) +
theme_classic(base_size = 18) +
xlab("") +
ylab("Overlap in acquired mutations\n(Paired Sørensen-Dice similarity)") +
theme(panel.spacing = unit(3, "lines"),
strip.text = element_text(size = 18, face = "bold")) +
ylim(0.1, 0.55)
    
print(p_day36)
    
## Identifying which genes mutated in which populations
# Assuming 'breseq_NS_noP0' is loaded in the environment
    
top_genes_in_study2_day36 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "Line"]))
colnames(top_genes_in_study2_day36) <- c("description", "Line", "hits")
top_genes_day36 <- names(tail(sort(table(top_genes_in_study2_day36[top_genes_in_study2_day36$hits > 0, "description"])), 20))
    
### Constructing a gene analysis table for P12
breseq_gene_analysis2_day36 <- data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "description"], breseq_NS_noP0[breseq_NS_noP0$Passage=="P12", "Line"]))
colnames(breseq_gene_analysis2_day36) <- c("gene", "Line", "hits")
    
# Define the resistance gene vectors
rfbA <- c("FMS6", "VCM4", "MS2", "FMS4", "MS10", "QAC4")  
glycosyl <- c("VCM19", "MS12", "FMS13", "SNK12", "QAC3", "SNK11", "VCM17", "SNK6", "MS1")  
glycoside <- c("SNK7", "FMS11", "FMS12", "FMS9", "FMS16")  
    
# Annotate each population according to its resistance gene
breseq_gene_analysis2_day36[breseq_gene_analysis2_day36$Line %in% rfbA, "res_gene"] <- "rfbA"
breseq_gene_analysis2_day36[breseq_gene_analysis2_day36$Line %in% glycosyl, "res_gene"] <- "glycosyl"
breseq_gene_analysis2_day36[breseq_gene_analysis2_day36$Line %in% glycoside, "res_gene"] <- "glycoside"
    
# Remove rows where res_gene is NA
breseq_gene_analysis2_day36 <- breseq_gene_analysis2_day36[!is.na(breseq_gene_analysis2_day36$res_gene), ]
    
### Calculate the percentage of populations with each resistance gene that mutated in each of the other genes
breseq_gene_props_day36 <- data.frame(matrix(nrow = 0, ncol = 3))
l_day36 <- list(rfbA = rfbA, glycosyl = glycosyl, glycoside = glycoside)
names(l_day36) <- c("rfbA mutants", "PSPTO_4988 mutants", "PSPTO_4991 mutants")
    
# Loop through genes and calculate mutation proportions
for (gene in top_genes_day36) {
 for (j in seq_along(l_day36)) {
res_gene <- l_day36[[j]]
        prop <- nrow(breseq_gene_analysis2_day36[breseq_gene_analysis2_day36$Line %in% res_gene & 
                                                   breseq_gene_analysis2_day36$gene == gene & 
                                                   breseq_gene_analysis2_day36$hits > 0, ]) / 
          nrow(breseq_gene_analysis2_day36[breseq_gene_analysis2_day36$Line %in% res_gene & 
                                             breseq_gene_analysis2_day36$gene == gene, ])
        breseq_gene_props_day36 <- rbind(breseq_gene_props_day36, c(gene, names(l_day36)[j], prop))
      }
    }
    
# Sort the 'res_gene' column so that "rfbA mutants" appears first, followed by "PSPTO_4988 mutants" and then "PSPTO_4991 mutants"
colnames(breseq_gene_props_day36)
colnames(breseq_gene_props_day36) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props_day36$res_gene <- factor(breseq_gene_props_day36$res_gene, levels = c("PSPTO_4991 mutants", "PSPTO_4988 mutants", "rfbA mutants"))
    
# Assign column names
colnames(breseq_gene_props_day36) <- c("top_gene", "res_gene", "percent_pops")
breseq_gene_props_day36$percent_pops <- as.numeric(breseq_gene_props_day36$percent_pops)
breseq_gene_props_day36[breseq_gene_props_day36$percent_pops == 0, "percent_pops"] <- 10^-4
    
# Replace problematic dashes and ensure proper encoding
breseq_gene_props_day36$top_gene <- gsub("‑", "-", breseq_gene_props_day36$top_gene)
    
# Maintain the order of factors on the X-axis
breseq_gene_props_day36$top_gene <- factor(breseq_gene_props_day36$top_gene, levels = rev(unique(breseq_gene_props_day36$top_gene)))
    
#--- Plot Fig 5B----
    ggplot(breseq_gene_props_day36, aes(top_gene, res_gene, fill = percent_pops * 100)) +
      geom_tile(color = "grey80") +
      theme_classic(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, family = "Arial")) +
      scale_fill_viridis(name = "Populations\nwith mutation (%)", limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
      ylab("") +
      xlab("") +
      ggtitle("Day 36 of experimental evolution")
    
# Saving objects needed for the figures
save(p_day6, p_day36, annotation_day6, test_gene_df_day6, annotation_day6, breseq_gene_props_day6, test_gene_df_day36, breseq_gene_props_day36, 
     file = "ReproHack_data_figures.RData")



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
