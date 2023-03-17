#title: "Tempeh-ITS"
#author: "Jonathan Sogin"
#date: "2023"


#Importing libraries
#######################################################
#pre-processing and data handling packages
library("phyloseq"); packageVersion("phyloseq")
library("PERFect"); packageVersion("PERFect")
library("decontam"); packageVersion("decontam")

#visualization packages
library("ggpubr"); packageVersion("ggpubr")
library("ggtext"); packageVersion("ggtext")
library("ggVennDiagram"); packageVersion("ggVennDiagram")

#data analysis packages
library("microbiome"); packageVersion("microbiome")
library("vegan"); packageVersion("vegan")
library("GUniFrac"); packageVersion("GUniFrac")

#setting seed
addTaskCallback(function(...) {set.seed(02221997);TRUE})
#######################################################


#Custom functions
#######################################################

#Function to create joined table of OTU taxonomy to sample abundances
summary_table <- function(physeq){
  arg_name <- deparse(substitute(physeq))
  if(!taxa_are_rows(physeq)){otu_table(physeq)<-t(otu_table(physeq))}
  table <- as.data.frame(cbind(tax_table(physeq), get_taxa(physeq, sample_names(physeq))))
  var_name <- paste("summary_table", arg_name, sep="_")
  assign(var_name, table, envir = globalenv())
}

#Function to agglomerate taxa, making sure that unreseolved taxa (NA) are not agglomerated.
#This is supposed to be builtin to phyloseq but it was not behaving correctly.
glom_tax <- function(physeq, rank){
  taxa_ranks <- matrix(seq(1,7), ncol=1)
    rownames(taxa_ranks) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rank_number <- taxa_ranks[rank, ]
  tmp_physeq <- physeq
  test <- is.na(tax_table(tmp_physeq)[,rank_number])
  if(!is.na(table(test)["TRUE"])){
    unresolved <- rownames(subset(tax_table(tmp_physeq), test))
    uniqueid_check <- c()
    for(otu in unresolved){
      name_rank=rank_number-1
      if(is.na(tax_table(tmp_physeq)[otu, name_rank])){name_rank=rank_number-2}
      tax_table(tmp_physeq)[otu, rank_number] <- paste0("UR ", tax_table(tmp_physeq)[otu, name_rank], " ", gsub("(.{1})(.{7})", "\\1", otu)) #renaming the taxa with UR higher_tax_rank .{4} from rownames, which are hash ids of the sequences; the chances of the four character code being the same are low
      if(length(setdiff(tax_table(tmp_physeq)[otu, rank_number], uniqueid_check))==0){
        print("some tax hashes were not unique enough to produce a unique unresolved species ID; check hashes")
        stop()
      uniqueid_check <- c(uniqueid_check, tax_table(tmp_physeq)[otu, rank_number])
      }
    }
  }
  tax_glom(tmp_physeq, taxrank=rank, NArm=FALSE)
}
#######################################################


#Importing Data
#######################################################
#file names as environmental variables
biom_file <- "./Qiime_Data/ITS.biom"
sequence_file <- "./Qiime_Data/ITS-dna-sequences.fasta"
metadata_file <- read.csv("./Qiime_Data/metadata.csv", fileEncoding="UTF-8-BOM")

#removing ITS negative controls from metadata_file
metadata <- subset(metadata_file, Sample_ID_ITS!="nd")
rownames(metadata) <- metadata$Sample_ID_ITS

#converting data to phyloseq objects
fungdata <- import_biom(biom_file, NULL, sequence_file)
fungdata <- merge_phyloseq(fungdata, sample_data(metadata))

#Dereplicating Product and Time_Point Variables
product_variable = as.character(get_variable(fungdata, "Product"))
time_point_variable = as.character(get_variable(fungdata, "Time_point"))
sample_data(fungdata)$dereplicated <- mapply(paste, product_variable, time_point_variable)

#applying correct order to phyloseq object
sample_data(fungdata)$Time_point=factor(sample_data(fungdata)$Time_point, levels = c("pre-inoculation", "post-inoculation", "early", "middle", "end", "packaged"))
sample_data(fungdata)$dereplicated=factor(sample_data(fungdata)$dereplicated, levels = c("soy pre-inoculation", "soy post-inoculation", "soy early", "soy middle", "soy end", "soy packaged", "multigrain pre-inoculation", "multigrain post-inoculation", "multigrain early", "multigrain middle", "multigrain end", "multigrain packaged"))
sample_data(fungdata)$Product=factor(sample_data(fungdata)$Product, levels = c("soy", "multigrain"))
#######################################################


#Renaming taxonomy 
#######################################################
#naming taxonomic ranks
colnames(tax_table(fungdata)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Removing extraneous characters from UNITE taxonomy
tax_table(fungdata) = gsub("[[:lower:]]\\_\\_", "", tax_table(fungdata)) 

#Removing underscores formatted as spaces
tax_table(fungdata) = gsub("_{1}", " ", tax_table(fungdata)) 

write.csv(summary_table(fungdata), "Fungdata_OTU_table.csv", quote=F)
#######################################################


#Filtering for fungi
#######################################################
#filtering features
fungdataFilt <- subset_taxa(fungdata,
    (Kingdom=="Fungi") & #retaining features classified as fungi
    (!is.na(Phylum)) #excluding features unclassified at Phylum level
    )
#######################################################




#summarizing number of reads after filtering
#######################################################
#table summarizing number of filtered reads
sum_table_filt <- cbind(
    sample_sums(fungdata),
    sample_sums(subset_taxa(fungdata, is.na(Phylum))),
    sample_sums(fungdataFilt)
    )
colnames(sum_table_filt) <- c("total unfiltered","unidentified at phylum", "total passing filters")
#######################################################





#Filtering based on decontamination with fix for taxa believed to be inaccurately classified
#######################################################
#separating extraction and PCR controls for processing
Samples <- subset_samples(fungdataFilt, TRUE_or_CONTROL=="TRUE SAMPLE")
Ext_controls <- subset_samples(fungdataFilt, Sample_ID_ITS=="Ext_neg_1_1_ITS" | Sample_ID_ITS=="Ext_neg_1_2_ITS")
PCR_controls <- subset_samples(fungdataFilt, Sample_ID_ITS=="PCR_neg_1_ITS" | Sample_ID_ITS=="PCR__neg_2_ITS")

#running decontam with extraction controls
#13 is hardcoded here for formatting reasons of the summary table
decontam_Ext <- merge_phyloseq(Samples, Ext_controls)
sample_data(decontam_Ext)$is.neg <- sample_data(decontam_Ext)$TRUE_or_CONTROL == "CONTROL SAMPLE"
contam_Ext <- isContaminant(decontam_Ext, batch="ITS_Amplification_Batch", method="prevalence", neg="is.neg")
contam_table_Ext <- as.data.frame(cbind(tax_table(decontam_Ext), contam_Ext, get_taxa(decontam_Ext, sample_names(decontam_Ext))))
#exporting a summary table
contam_table_Ext <- contam_table_Ext[with(contam_table_Ext, order(-contaminant, p)),]
#adding some metadata to the table 
  sample_sums_Ext <- t(matrix(c(rep(NA, 13), sample_sums(Samples), sample_sums(Ext_controls))))
    colnames(sample_sums_Ext) <- c(colnames(contam_table_Ext)[1:13], names(sample_sums(Samples)), names(sample_sums(Ext_controls)))
    rownames(sample_sums_Ext) <- "Sample_Sums"
  samples_data <- t(sample_data(Samples)[, c("Product", "Time_point", "ITS_Amplification_Batch")])
  Ext_data <- t(sample_data(Ext_controls)[, c("Product", "Time_point", "ITS_Amplification_Batch")])
  blank <- matrix(rep(NA, 13*3), ncol=13)
    colnames(blank) <- colnames(contam_table_Ext)[1:13]
    rownames(blank) <- rownames(samples_data)
  samples_data <- cbind(blank, samples_data, Ext_data)
contam_table_Ext <- rbind(samples_data, sample_sums_Ext, contam_table_Ext)
contam_table_Ext['e6d3abfbd1efd8747548b863efdb4743','contaminant'] <- T #after looking at data this appeared to be a contaminant
contam_table_Ext['ec303f053b6e878db7ee7cf2e5270c64','contaminant'] <- F #after looking at data this did not appear to be a contaminant
write.csv(contam_table_Ext, "ITS_contams_Ext.csv", quote = F)



#running decontam with PCR controls
decontam_PCR <- merge_phyloseq(Samples, PCR_controls)
sample_data(decontam_PCR)$is.neg <- sample_data(decontam_PCR)$TRUE_or_CONTROL == "CONTROL SAMPLE"
contam_PCR <- isContaminant(decontam_PCR, batch="ITS_Amplification_Batch", method="prevalence", neg="is.neg")
contam_table_PCR <- as.data.frame(cbind(tax_table(decontam_PCR), contam_PCR, get_taxa(decontam_PCR, sample_names(decontam_PCR))))
#exporting a summary table
contam_table_PCR <- contam_table_PCR[with(contam_table_PCR, order(-contaminant, p)),]
#contam_table_PCR <- subset(contam_table_PCR, contam_table_PCR[,"PCR_neg_1_ITS"]!=0 | contam_table_PCR[,"PCR_neg_2_ITS"]!=0)
  #adding some metadata to the table  
  sample_sums_PCR <- t(matrix(c(rep(NA, 13), sample_sums(Samples), sample_sums(PCR_controls))))
    colnames(sample_sums_PCR) <- c(colnames(contam_table_PCR)[1:13], names(sample_sums(Samples)), names(sample_sums(PCR_controls)))
    rownames(sample_sums_PCR) <- "Sample_Sums"
  samples_data <- t(sample_data(Samples)[, c("Product", "Time_point", "ITS_Amplification_Batch")])
  PCR_data <- t(sample_data(PCR_controls)[, c("Product", "Time_point", "ITS_Amplification_Batch")])
  blank <- matrix(rep(NA, 13*3), ncol=13)
    colnames(blank) <- colnames(contam_table_PCR)[1:13]
    rownames(blank) <- rownames(samples_data)
  samples_data <- cbind(blank, samples_data, PCR_data)
contam_table_PCR <- rbind(samples_data, sample_sums_PCR, contam_table_PCR)
write.csv(contam_table_PCR, "ITS_contams_PCR.csv", quote = F)


#pooling contaminants from Extraction and PCR and removing OTUs identified as contaminants
contams <- union(row.names(subset(contam_table_Ext, contaminant==T)), row.names(subset(contam_table_PCR, contaminant==T)))
fungdata_decontam <- prune_taxa(setdiff(taxa_names(Samples), contams), fungdataFilt)

#summary table describing reads retained through pipeline
decontam_counts <- as.data.frame(sample_sums(fungdata_decontam))
sum_table <- cbind(sum_table_filt, as.data.frame(sample_sums(fungdata_decontam)))
colnames(sum_table) <- c("total unfiltered", "unidentified at phylum", "passing filters", "contaminants removed/clean")
sum_table
summary(sum_table[, "contaminants removed/clean"])

fungdata_decontam <- subset_samples(fungdata_decontam, TRUE_or_CONTROL=="TRUE SAMPLE")
summary(sample_sums(fungdata_decontam))

#removing taxa with zero OTUs because negative controls are not included
fungdata_decontam <- prune_taxa(taxa_sums(fungdata_decontam)!=0, fungdata_decontam)
summary_table(fungdata_decontam)

write.csv(summary_table(fungdata_decontam), "Fungdata_decontamed_table.csv", quote=F)
write.csv(sum_table, "Fungdata_filtering_counts.csv", quote=F)
#######################################################


#Relative Abundance Plotting
#######################################################
samples_plotting <- fungdata_decontam

#agglomerating data to Family level for plotting
fungdata_decontamGlomFamily <- glom_tax(samples_plotting, "Family")
fungdata_decontamGlomFamilyNorm <- transform_sample_counts(fungdata_decontamGlomFamily, function(x) x / sum(x))

plotting <- fungdata_decontamGlomFamilyNorm

sample_data(plotting)$Product <- as.factor(sample_data(plotting)$Product)
sample_data(plotting)$Time_point <- as.factor(sample_data(plotting)$Time_point)
sample_data(plotting)$dereplicated <- as.factor(sample_data(plotting)$dereplicated)
sample_data(plotting) <- sample_data(plotting)[,c("Product", "Time_point", "dereplicated", "Replicate")]
pruned = filter_taxa(plotting, function(x) mean(x) > 0.025, T)

otu_other = matrix(1-sample_sums(pruned), nrow=1)
  colnames(otu_other) = sample_names(pruned)
  rownames(otu_other) = "other"
  otu_plot = rbind(otu_table(pruned), otu_other)
tax_other = matrix(rep("other",7), nrow=1)
  rownames(tax_other)="other"
  colnames(tax_other) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_plot = rbind(tax_table(pruned), tax_other)
otu_object = otu_table(otu_plot, taxa_are_row=T)
tax_object = tax_table(tax_plot)
plot_object= phyloseq(otu_object, tax_object, sample_data(pruned))

plottable = psmelt(plot_object)
  plottable$Sample =factor(plottable$Sample)
  plottable$Time_point=factor(plottable$Time_point, levels = c("pre-inoculation", "post-inoculation", "early", "middle", "end", "packaged"))
  plottable$Product=factor(plottable$Product, levels = c("soy", "multigrain"))
  plottable$Replicate=factor(plottable$Replicate, levels = c(1, 2, 3, 4, 5, 6))
  plottable$Family[plottable$Family!="other"]=gsub("^", "<i>", plottable$Family[plottable$Family!="other"])
  plottable$Family[plottable$Family!="other"]=gsub("$", "</i>", plottable$Family[plottable$Family!="other"])
  plottable$Family=factor(plottable$Family, levels=c(setdiff(plottable$Family, "other"), "other"))
retail <- subset(plottable, Product=="soy")
multigrain <- subset(plottable, Product=="multigrain")

retail_rel_abundance_plot <- ggbarplot(retail, x="Replicate", y="Abundance", fill="Family", width=1, xlab=F, ylab="soy", palette = "simpsons", )+ 
  facet_grid(~Time_point, switch="x", scales="fixed", space="free")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  font("ylab", face="bold", size=15)+
  font("legend.title")+
  theme(text=element_text(family="serif"))+
  theme(legend.text=element_markdown())
multigrain_rel_abundance_plot <- ggbarplot(multigrain, x="Replicate", y="Abundance", fill="Family", width=1, xlab=F, ylab="multigrain",palette = "simpsons")+
  facet_grid(~Time_point, switch="x", scales="fixed", space="free")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  font("ylab", face="bold", size=15)+
  font("legend.title")+
  theme(text=element_text(family="serif"))+
  theme(legend.text=element_markdown())
combined_rel_abundance_plot <- ggarrange(retail_rel_abundance_plot, multigrain_rel_abundance_plot, ncol=1, common.legend = T, legend = "bottom", align="v")
#combined_rel_abundance_plot <- annotate_figure(combined_rel_abundance_plot, top = text_grob("Fungal Family Relative Abundance", face = "bold", family="serif", size=18))
#exported as tiff Fungal_Relative_Abundance
ggsave(plot=combined_rel_abundance_plot, filename="Fungal_Relative_Abundance.tiff", width=8, height=6, units="in", dpi="print")
#######################################################


#Data analysis
#subsetting and denoising data
#######################################################
#subsetting data
#not including beginning samples due to the presence of several zeros across the data
fungdata_decontam_nostart <- subset_samples(fungdata_decontam, Time_point!="pre-inoculation" & Time_point!="post-inoculation")
  fungdata_decontam_nostart <- filter_taxa(fungdata_decontam_nostart, function(x) sum(x)!=0, prune=T)

#Denoising data
#Using PERFect, which is a statistical means of eliminating taxa that do not contribute to covariance
#this will inherently affect alpha diversity measurements, but will do so with the benefit of greater confidence in preventing artifically high alpha diversity measurements due to sequencing artifacts

#splitting up fungdata_decontam_nostart object, as the otu table will be modified and the separate parts will need to be merged into a new objectafter denoising
fungdatadecontam_otu <- otu_table(fungdata_decontam_nostart)
fungdatadecontam_tax <- tax_table(fungdata_decontam_nostart)
fungdatadecontam_sample <- sample_data(fungdata_decontam_nostart)
fungdatadecontam_seqs <- refseq(fungdata_decontam_nostart)

#running PERFseq
#transposing data to use in PERFseq
Counts <- t(fungdatadecontam_otu)
dim(Counts)

res_sim <- PERFect_sim(X = Counts)
  dim(res_sim$filtX)  
simultaneous_filtering_plot <- pvals_Plots(PERFect = res_sim, X = Counts, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)
  simultaneous_filtering_plot <- simultaneous_filtering_plot$plot + ggtitle("Simultanenous Filtering")

res_perm <- PERFect_perm(X = Counts, Order = "pvals", pvals_sim = res_sim, algorithm = "full")
  dim(res_perm$filtX)
permutational_filtering_plot <- pvals_Plots(res_perm, Counts)
  permutational_filtering_plot <- permutational_filtering_plot$plot + ggtitle("Full Algorithm")

new_otu <- t(res_perm$filtX)

#remerging data to phyloseq object
fungdata_denoised <- phyloseq(new_otu, fungdatadecontam_tax, fungdatadecontam_sample, fungdatadecontam_seqs)

#table summarizing number of filtered reads
sum_table_denoised <- cbind(sum_table[sample_names(fungdata_denoised),], data.frame(unlist(sample_sums(fungdata_denoised))))
colnames(sum_table_denoised) <- c("total unfiltered", "unidentified at phylum", "passing filters", "contaminants removed/clean", "denoised")
sum_table_denoised

write.csv(summary_table(fungdata_denoised), "Fungdata_denoised_table.csv", quote=F)
write.csv(sum_table_denoised, "Fungdata_filtering_counts_denoised.csv", quote=F)
#######################################################


#######################################################
analysis <- fungdata_denoised

#Rarefying data for analysis
#rarefying data for alpha and beta diversity analyses
min_analysis <- summary(sample_sums(analysis))['Min.']
analysis_rare <- rarefy_even_depth(analysis, rngseed=02221997, sample.size=min_analysis, replace=F)

#normalizing data and agglomerating to species for differential abundance analysis
analysis_GlomSpecies <- glom_tax(analysis, "Species")
analysis_GlomSpecies_norm <-transform_sample_counts(analysis_GlomSpecies, function(x) x / sum(x))
#######################################################


#Alpha diversity analysis
#######################################################
#calculating alpha diversity indices
phyloseq4alpha <- analysis_rare
fungdata_diversity <- alpha(phyloseq4alpha, index=c("inverse_simpson", "shannon", "coverage", "pielou", "observed", "chao1", "relative", "core_abundance", "low_abundance"))
indices <- names(fungdata_diversity)
index_names <- as.list(c("Observed", "Chao1", "Inverse Simpson", "Shannon-Wiener", "Coverage (50%)", "Pielou", "Relative", "Core Abundance", "Low Abundance"))
names(index_names)<-indices
fungdata_diversity <- cbind(sample_data(phyloseq4alpha), fungdata_diversity)
fungdata_diversity$Time_point=factor(fungdata_diversity$Time_point, levels = c("pre-inoculation", "post-inoculation", "early", "middle", "end", "packaged"))
fungdata_diversity$dereplicated=factor(fungdata_diversity$dereplicated, levels = c("soy pre-inoculation", "soy post-inoculation", "soy early", "soy middle", "soy end", "soy packaged", "multigrain pre-inoculation", "multigrain post-inoculation", "multigrain early", "multigrain middle", "multigrain end", "multigrain packaged"))
fungdata_diversity$Product=factor(fungdata_diversity$Product, levels = c("soy", "multigrain"))
#fungdata_diversity is a data frame holding all the alpha diversity calculations

measures <- c("Product", "Time_point", "dereplicated")
plots <- c()
stats <- c()
for(i in indices){
  plots[[i]] <- ggboxplot(fungdata_diversity, x="Product", y=i, fill="Product", palette="simpsons", ylab=index_names[i], xlab=F, add="point", add.params=list(color="darkgrey"))+rremove("x.text")+theme(text=element_text(family="serif"))
  pvals <- c()
for(m in measures){
    stats[[i]][[m]] <- kruskal.test(unlist(fungdata_diversity[i]), unlist(fungdata_diversity[m]))
    pvals[m] <- stats[[i]][[m]]$p.value
  }
  stats[[i]][["p.adjust"]] <- p.adjust(pvals, method="fdr")
  }

combined_diversity_plot <- ggarrange(plots$diversity_shannon, plots$diversity_inverse_simpson, plots$diversity_coverage, nrow=1, legend="none")
  combined_diversity_plot <- annotate_figure(combined_diversity_plot, left=text_grob("Diversity", face="bold", rot=90, size=14, family="serif"))
combined_richness_plot <- ggarrange(plots$observed, plots$chao1, nrow=1, legend="none")
  combined_richness_plot <- annotate_figure(combined_richness_plot, left=text_grob("Richness", face="bold", rot=90, size=14, family="serif"))
combined_evenness_plot <- ggarrange(plots$evenness_pielou, nrow=1, legend="none")
  combined_evenness_plot <- annotate_figure(combined_evenness_plot, left=text_grob("Evenness", face="bold", rot=90, size=14, family="serif"))
combined_richness_evenness_plot <- ggarrange(combined_richness_plot, combined_evenness_plot, nrow=1, legend="none",widths =c(2, 1))
combined_dominance_plot <- ggarrange(plots$dominance_relative, plots$dominance_core_abundance, nrow=1, legend="none")
  combined_dominance_plot <- annotate_figure(combined_dominance_plot, left=text_grob("Dominance", face="bold", rot=90, size=14, family="serif"))
combined_rarity_plot <- ggarrange(plots$rarity_low_abundance, nrow=1, legend="none")
  combined_rarity_plot <- annotate_figure(combined_rarity_plot, left=text_grob("Rarity", face="bold", rot=90, size=14, family="serif"))
combined_dominance_rarity_plot <- ggarrange(combined_dominance_plot, combined_rarity_plot, nrow=1, legend="none", widths=c(2, 1))
combined_alpha_diversity_plot <- ggarrange(combined_diversity_plot, combined_richness_evenness_plot, combined_dominance_rarity_plot, common.legend=T, legend="bottom", legend.grob=get_legend(plots$diversity_shannon), ncol=1)
#  combined_alpha_diversity_plot <- annotate_figure(combined_alpha_diversity_plot, top=text_grob("Fungal Alpha Diversity Measures", face="bold", size=18, family="serif"))
#exported as tiff Fungal_Alpha_Diversity_Product
ggsave(plot=combined_alpha_diversity_plot, filename="Fungal_Alpha_Diversity.tiff", width=8, height=6, units="in", dpi="print")
#######################################################


#Beta diversity ordination
#######################################################
#Plotting beta diversity ordination
phyloseq4ordination <- analysis_GlomSpecies_norm
ordination_NMDS <- ordinate(phyloseq4ordination, "NMDS", "bray")
ordination_plot_data <- plot_ordination(phyloseq4ordination, ordination_NMDS, type="biplot", justDF=T)
ordination_plot_data$Time_point=factor(ordination_plot_data$Time_point, levels = c("pre-inoculation", "post-inoculation", "early", "middle", "end", "packaged"))
ordination_plot_data$Product=factor(ordination_plot_data$Product, levels = c("soy", "multigrain"))
ordination_plot_samples <- subset(ordination_plot_data, id.type=="Samples")
ordination_plot_taxa <- subset(ordination_plot_data, id.type=="Taxa")
ordination_plot_product <- ggscatter(data=ordination_plot_samples, x="NMDS1", y="NMDS2", color="Product", palette="simpsons", legend="right", size=3, ellipse=T, ellipse.alpha=0.25, star.plot=T, ellipse.type="convex")+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank())+
  theme(legend.title=element_text(face="bold"))+
  theme(text=element_text(family="serif"))
ordination_plot_time_point <- ggscatter(data=ordination_plot_samples, x="NMDS1", y="NMDS2", color="Time_point", palette="futurama", legend="right", size=3, ellipse=T, ellipse.alpha=0.25, star.plot=T, ellipse.type="convex")+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank())+
  theme(legend.title=element_text(face="bold"))+
  theme(text=element_text(family="serif"))
combined_ordination_plot <- ggarrange(ordination_plot_product, ordination_plot_time_point, ncol=1, align="v")
#  labs(title="Fungal Species Bray-Curtis Ordination")+
#  theme(plot.title=element_text(family="serif", face="bold", size=18))
#exported as tiff Fungal_Bray_Ordination
ggsave(plot=combined_ordination_plot, filename="Fungal_Bray_Ordination.tiff", width=8, height=6, units="in", dpi="print")
#######################################################

#Beta diversity PERMANOVA
#######################################################
#Conducting PERMANOVA analysis to determine effect of Product and Time_point variables on overall community structure

#Agglomerating data to Species level to reduce the dimensionality of the data
phyloseq4permanova <- analysis_GlomSpecies_norm

otu <- t(abundances(phyloseq4permanova))
dist <- vegdist(otu, method="bray")
meta <- meta(phyloseq4permanova)

#comparing reduced models for the interaction of Product and Time_point and the main effects
dbrda.main <- dbrda(dist ~ Product + Time_point, data=meta)
dbrda.main.aov <- anova(dbrda.main, by = 'margin')
dbrda.main.aov

dbrda.interaction <- dbrda(dist ~ Product*Time_point, data=meta)
dbrda.interaction.aov <- anova(dbrda.interaction, by = 'margin')
dbrda.interaction.aov

#running permanova
#this one is the influence of time with product controlled
permanova_Product <- adonis(otu ~ Time_point + Product, data=meta, method="bray", permutations=10000)
permanova_Product[c("call", "aov.tab")]

#this one is the influence of time with product controlled
permanova_Time_point <- adonis(otu ~ Product + Time_point, data=meta, method="bray", permutations=10000)
permanova_Time_point[c("call", "aov.tab")]

#checking condition for equal variance around centroid
permdisp2_Product <- betadisper(dist, meta$Product)
permdisp2_Product.aov <- anova(permdisp2_Product)
permdisp2_Product$call
permdisp2_Product.aov

permdisp2_Time_point <- betadisper(dist, meta$Time_point)
permdisp2_Time_point.aov <- anova(permdisp2_Time_point)
permdisp2_Time_point$call
permdisp2_Time_point.aov

coef_Product <- coefficients(permanova_Product)["Product1",]
top_coef_Product <- as.data.frame(coef_Product[rev(order(abs(coef_Product)))])
names(top_coef_Product) <- "coefficient"
top_coef_Product$otu_id <- row.names(top_coef_Product)
for(i in 1:nrow(top_coef_Product)){
  otuid = top_coef_Product$otu_id[i]
  species = tax_table(phyloseq4permanova)[otuid][[7]]
  species_name = paste0("<i>", species, "</i>")
    species_name = gsub("^<i>UR ", "UR <i>", species_name)
    species_name = gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", species_name)
  top_coef_Product$species_name[i] <- species_name
}
top_coef_Product <- top_coef_Product[order(top_coef_Product$coefficient),]

#plotting the top product species coefficients
plot_coef_Product <- ggbarplot(data=top_coef_Product, x="species_name", y="coefficient", fill="coefficient", orientation="horiz", legend="none")+
  scale_fill_gradient(low="#FEDD61", high="#9BB8EA", limits=c(-0.06, 0.065))+
  theme(axis.text.y=element_markdown())+
  scale_y_continuous(breaks=seq(-0.6, 0.06, 0.02), limits=c(-0.06, 0.065))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(text=element_text(family="serif"))

#exported as tiff Fungal_Product_PERMANOVA_Coefficients.tiff
ggsave(plot=plot_coef_Product, filename="Fungal_Product_PERMANOVA_Coefficients.tiff", width=8, height=3, units="in", dpi="print")
#######################################################


#Using ZicoSeq to do differential abundance testing
#######################################################
phyloseq4ZicoSeq <- analysis_GlomSpecies
  phyloseq4ZicoSeq_norm <- analysis_GlomSpecies_norm

zico_meta_data <-data.frame(sample_data(phyloseq4ZicoSeq))
colnames(zico_meta_data) <- colnames(sample_data(phyloseq4ZicoSeq))

feature_vals <- as.integer(otu_table(phyloseq4ZicoSeq))
  dims <- dim(otu_table(phyloseq4ZicoSeq))
zico_otu <- matrix(feature_vals, nrow=dims[1])
  rownames(zico_otu) <- as.character(rownames(otu_table(phyloseq4ZicoSeq)))
  colnames(zico_otu) <- as.character(colnames(otu_table(phyloseq4ZicoSeq)))

zico_meta_data <- data.frame(sample_data(phyloseq4ZicoSeq))

#ZicoSeq filtering parameters consistent between Time_point and Product analyses
  meta.dat=zico_meta_data
  feature.dat=zico_otu
  prev.filter=0.33
  mean.abund.filter=0.0005
  max.abund.filter=0
  min.prop=0
  outlier.pct=0.03
  perm.no=10000

#running ZicoSeq to distinguish taxa different by Product
ZicoSeq_Product <- ZicoSeq(meta.dat=meta.dat, feature.dat=feature.dat, grp.name="Product", adj.name="Time_point", feature.dat.type="count", prev.filter=prev.filter, mean.abund.filter=mean.abund.filter, max.abund.filter=max.abund.filter, min.prop=min.prop, winsor.end="top", perm.no=perm.no, return.feature.dat=T)

  otus_tested_Product <- names(ZicoSeq_Product$p.adj.fdr)
  rel_proportion_tested <- sample_sums(prune_taxa(otus_tested_Product, phyloseq4ZicoSeq_norm))

res_product <- data.frame(cbind(ZicoSeq_Product$p.raw, ZicoSeq_Product$p.adj.fdr))
  colnames(res_product) <- c("p.raw", "p.adj.fdr")
  res_product <- res_product[order(res_product[,"p.adj.fdr"]),]
  res_product_tax <- tax_table(phyloseq4ZicoSeq)[rownames(res_product),]
  res_product_mean_soy <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Product=="soy"))[rownames(res_product)])
  res_product_cv_soy <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Product=="soy"))[rownames(res_product)], 1, function(x) sd(x)/mean(x))
  res_product_mean_multigrain <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Product=="multigrain"))[rownames(res_product)])
  res_product_cv_multigrain <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Product=="multigrain"))[rownames(res_product)], 1, function(x) sd(x)/mean(x))
  res_product <- cbind(res_product, res_product_tax, res_product_mean_soy, res_product_cv_soy, res_product_mean_multigrain, res_product_cv_multigrain)
  res_product <- data.frame(res_product)
write.csv(res_product, "ITS_ZicoSeq_Product.csv", quote=F)

#running ZicoSeq to distinguish taxa different by Time_point
ZicoSeq_Time_point <- ZicoSeq(meta.dat=meta.dat, feature.dat=feature.dat, grp.name="Time_point", adj.name="Product", feature.dat.type="count", prev.filter=prev.filter, mean.abund.filter=mean.abund.filter, max.abund.filter=max.abund.filter, min.prop=min.prop, winsor.end="top", perm.no=perm.no, return.feature.dat=T)

  otus_tested_Time_point <- names(ZicoSeq_Time_point$p.adj.fdr)
  rel_proportion_tested <- sample_sums(prune_taxa(otus_tested_Time_point, phyloseq4ZicoSeq_norm))

res_Time_point <- data.frame(cbind(ZicoSeq_Time_point$p.raw, ZicoSeq_Time_point$p.adj.fdr))
  colnames(res_Time_point) <- c("p.raw", "p.adj.fdr")
  res_Time_point <- res_Time_point[order(res_Time_point[,"p.adj.fdr"]),]
  res_Time_point_tax <- tax_table(phyloseq4ZicoSeq)[rownames(res_Time_point),]
  res_Time_point_mean_early <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Time_point=="early"))[rownames(res_Time_point)])
  res_Time_point_cv_early <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Time_point=="early"))[rownames(res_Time_point)], 1, function(x) sd(x)/mean(x))
  res_Time_point_mean_middle <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Time_point=="middle"))[rownames(res_Time_point)])
  res_Time_point_cv_middle <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Time_point=="middle"))[rownames(res_Time_point)], 1, function(x) sd(x)/mean(x))
  res_Time_point_mean_end <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Time_point=="end"))[rownames(res_Time_point)])
  res_Time_point_cv_end <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Time_point=="end"))[rownames(res_Time_point)], 1, function(x) sd(x)/mean(x))
  res_Time_point_mean_packaged <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Time_point=="packaged"))[rownames(res_Time_point)])
  res_Time_point_cv_packaged <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Time_point=="packaged"))[rownames(res_Time_point)], 1, function(x) sd(x)/mean(x))
  res_Time_point <- cbind(res_Time_point, res_Time_point_tax, res_Time_point_mean_early, res_Time_point_cv_early, res_Time_point_mean_middle, res_Time_point_cv_middle, res_Time_point_mean_end, res_Time_point_cv_end, res_Time_point_mean_packaged, res_Time_point_cv_packaged)
  res_Time_point <- data.frame(res_Time_point)
write.csv(res_Time_point, "ITS_ZicoSeq_Time_point.csv", quote=F)
#######################################################


#Core microbiome plotting
#######################################################
phyloseq4core <- analysis_GlomSpecies_norm

#venn diagram core microbiome https://microbiome.github.io/tutorials/core_venn.html
#time based core analysis
time_states <- sort(factor(unique(as.character(sample_data(phyloseq4core)$Time_point)), levels = c("pre-inoculation", "post-inoculation", "early", "middle", "end", "packaged")))

list_core_time <- c() # an empty object to store information
for (n in time_states){ # for each variable n in product_states
    ps.sub <- subset_samples(phyloseq4core, Time_point == n) # Choose sample from product_states by n
    core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                           detection = 0.005, # 0.001 in atleast 90% samples 
                           prevalence = 0.66) #detected in at least 75% of samples
    print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each product_states
    list_core_time[[n]] <- core_m # add to a list core taxa for each group.
}

variables_time <- names(list_core_time)
tax_core_time <- c()
for (n in variables_time){ # for each variable n in variables
    list <- unlist(list_core_time[n])
    tax_list <- c()
    for (i in list){
        tax_list <- rbind(tax_list, tax_table(phyloseq4core)[i])
    }
    tax_core_time[[n]] <- tax_list
}

tax_core <- unique(unlist(list_core_time))

phyloseq_core_time <- prune_taxa(intersect(tax_core, taxa_names(phyloseq4core)), phyloseq4core)
######manually coloring numbers and species########
time_intersect <- tax_table(phyloseq_core_time)[intersect(rownames(tax_core_time$`early`), intersect(rownames(tax_core_time$`middle`), intersect(rownames(tax_core_time$`end`), rownames(tax_core_time$packaged))))]
unique_middleend <- tax_table(phyloseq_core_time)[setdiff(rownames(tax_core_time$`middle`), rownames(time_intersect))]
time_coloring <- rev(c("red", "red", "blue", "black"))


###################################################
otu_other = matrix(1-sample_sums(phyloseq_core_time), nrow=1)
  colnames(otu_other) = sample_names(phyloseq_core_time)
  rownames(otu_other) = "other"
  otu_plot = rbind(otu_table(phyloseq_core_time), otu_other)
tax_other = matrix(rep("other",7), nrow=1)
  rownames(tax_other)="other"
  colnames(tax_other) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_plot = rbind(tax_table(phyloseq_core_time), tax_other)
otu_object = otu_table(otu_plot, taxa_are_row=T)
tax_object = tax_table(tax_plot)
phyloseq_core_time_plotobject= phyloseq(otu_object, tax_object, sample_data(phyloseq_core_time))
heat_data_time <- merge_samples(phyloseq_core_time_plotobject, "Time_point", fun="mean")
heat_data_time_norm <- transform_sample_counts(heat_data_time, function(x) -log(x / sum(x)))
otu_tab_time <- psmelt(heat_data_time_norm)
for(i in 1:nrow(otu_tab_time)){
  otuid = otu_tab_time$OTU[i]
  species = tax_table(heat_data_time_norm)[otuid][[7]]
  species_name = paste0("<i>", species, "</i>")
      species_name = gsub("^<i>UR ", "UR <i>", species_name)
      species_name = gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", species_name)
  if(otuid=="other"){species_name="other"}
  otu_tab_time$species_name[i] <- species_name
}
otu_tab_time$Sample=factor(otu_tab_time$Sample, levels = c("pre-inoculation", "post-inoculation", "early", "middle", "end", "packaged"))
otu_tab_time$OTU=factor(otu_tab_time$OTU, levels=c(names(rev(sort(taxa_sums(phyloseq_core_time)))), "other"))
otu_tab_time=otu_tab_time[order(otu_tab_time$OTU, otu_tab_time$Sample),]
otu_tab_time$species_name=factor(otu_tab_time$species_name, levels=rev(unique(otu_tab_time$species_name)))
heat_plot_time <- ggplot(data=otu_tab_time, aes(x=Sample, y=species_name, fill=Abundance))+
  theme_classic()+
  scale_fill_gradient(high="#FED439", low="#709AE1", limits=c(0, 10))+
  geom_tile()+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.line=element_blank())+
  theme(text=element_text(family="serif"),axis.text.y=element_markdown(color=time_coloring, size=rel(1.2)), axis.text.x=element_text(size=rel(1.2)))+
  theme(legend.position="right")+
  labs(fill = "ln(abund)")+
  labs(title="Time Point")+
  theme(plot.title=element_text(family="serif", size=15, hjust=0.5))

venn_time <- Venn(list_core_time)
data_time <- process_data(venn_time)
venn_plot_time <- ggplot()+
  theme_classic()+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.line=element_blank(), axis.text=element_blank())+
  geom_sf(aes(fill=count), data=venn_region(data_time), show.legend=FALSE)+  # region coloring
  geom_sf(color="darkgrey", lwd=0.25, data=venn_setedge(data_time), show.legend=FALSE)+  # region edges
  geom_sf_text(aes(label=name), fontface="plain", family="serif", data=venn_setlabel(data_time), vjust=c(0,-0.02,-0.1,0))+  # group labels
  geom_sf_label(aes(label=count), data=subset(venn_region(data_time), count>0), fontface="bold", family="serif", alpha=1, color=c("blue", "red"))+ # region labels
  scale_fill_gradient(low="#f8fafe", high="#7EA4E4", limits=c(0,10))+
  scale_x_continuous(expand = c(0.1,0.1))

#product based core analysis
product_states <- sort(factor(unique(as.character(sample_data(phyloseq4core)$Product)), levels = c("soy", "multigrain")))

list_core_product <- c() # an empty object to store information
for (n in product_states){ # for each variable n in product_states
    ps.sub <- subset_samples(phyloseq4core, Product == n) # Choose sample from product_states by n
    core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                           detection = 0.005, # 0.001 in atleast 90% samples 
                           prevalence = 0.66) #detected in at least 75% of samples
    print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each product_states
    list_core_product[[n]] <- core_m # add to a list core taxa for each group.
}

variables_product <- names(list_core_product)
tax_core_product <- c()
for (n in variables_product){ # for each variable n in variables
    list <- unlist(list_core_product[n])
    tax_list <- c()
    for (i in list){
        tax_list <- rbind(tax_list, tax_table(phyloseq4core)[i])
    }
    tax_core_product[[n]] <- tax_list
}

tax_core <- unique(unlist(list_core_product))

phyloseq_core_product <- prune_taxa(intersect(tax_core, taxa_names(phyloseq4core)), phyloseq4core)
######manually coloring numbers and species########
product_intersect <- tax_table(phyloseq_core_product)[intersect(rownames(tax_core_product$`soy`), rownames(tax_core_product$`multigrain`))]
product_coloring <- rev(c("red", "red", "red", "black"))


###################################################
otu_other = matrix(1-sample_sums(phyloseq_core_product), nrow=1)
  colnames(otu_other) = sample_names(phyloseq_core_product)
  rownames(otu_other) = "other"
  otu_plot = rbind(otu_table(phyloseq_core_product), otu_other)
tax_other = matrix(rep("other",7), nrow=1)
  rownames(tax_other)="other"
  colnames(tax_other) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_plot = rbind(tax_table(phyloseq_core_product), tax_other)
otu_object = otu_table(otu_plot, taxa_are_row=T)
tax_object = tax_table(tax_plot)
phyloseq_core_product_plotobject= phyloseq(otu_object, tax_object, sample_data(phyloseq_core_product))
heat_data_product <- merge_samples(phyloseq_core_product_plotobject, "Product", fun="mean")
heat_data_product_norm <- transform_sample_counts(heat_data_product, function(x) -log(x / sum(x)))
otu_tab_product <- psmelt(heat_data_product_norm)
for(i in 1:nrow(otu_tab_product)){
  otuid = otu_tab_product$OTU[i]
  species = tax_table(heat_data_product_norm)[otuid][[7]]
  species_name = paste0("<i>", species, "</i>")
      species_name = gsub("^<i>UR ", "UR <i>", species_name)
      species_name = gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", species_name)
  if(otuid=="other"){species_name="other"}
  otu_tab_product$species_name[i] <- species_name
}
otu_tab_product$Sample=factor(otu_tab_product$Sample, levels = c("soy", "multigrain"))
otu_tab_product$OTU=factor(otu_tab_product$OTU, levels=c(names(rev(sort(taxa_sums(phyloseq_core_product)))), "other"))
otu_tab_product=otu_tab_product[order(otu_tab_product$OTU, otu_tab_product$Sample),]
otu_tab_product$species_name=factor(otu_tab_product$species_name, levels=rev(unique(otu_tab_product$species_name)))
heat_plot_product <- ggplot(data=otu_tab_product, aes(x=Sample, y=species_name, fill=Abundance))+
  theme_classic()+
  scale_fill_gradient(high="#FED439", low="#709AE1", limits=c(0, 10))+
  geom_tile()+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.line=element_blank())+
  theme(text=element_text(family="serif"), axis.text.y=element_markdown(color=product_coloring, size=rel(1.2)), axis.text.x=element_text(size=rel(1.2)))+
  theme(legend.position="right")+
  labs(fill = "-ln(abund)")+
  labs(title="Product")+
  theme(plot.title=element_text(family="serif", size=15, hjust=0.5))

venn_product <- Venn(list_core_product)
data_product <- process_data(venn_product)
venn_plot_product <- ggplot()+ 
theme_classic()+
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.line=element_blank(), axis.text=element_blank())+
  geom_sf(aes(fill=count), data=venn_region(data_product), show.legend=FALSE)+  # region coloring
  geom_sf(color="darkgrey", lwd=0.25, data=venn_setedge(data_product), show.legend=FALSE)+  # region edges
  geom_sf_text(aes(label=name), fontface="plain", family="serif", data=venn_setlabel(data_product))+  # group labels
  geom_sf_label(aes(label=count), data=subset(venn_region(data_product), count>0), fontface="bold", family="serif", alpha=1, color="red")+ # region labels
  scale_fill_gradient(low="#f8fafe", high="#7EA4E4", limits=c(0,10))+
  scale_x_continuous(expand=c(0.2,0.2))

combined_heat_plots <- ggarrange(heat_plot_product, heat_plot_time, align="h", widths=c(2,3), common.legend=T, legend="bottom")
combined_venn_plots <- ggarrange(venn_plot_product, venn_plot_time)
combined_core_plot <- ggarrange(combined_venn_plots, combined_heat_plots, ncol=1, heights=c(2,3), align="v")
#  labs(title="Fungal Core Species")+

#Exported as tiff Fungal_Core
ggsave(plot=combined_core_plot, filename="Fungal_Core.tiff", width=8, height=6, units="in", dpi="print")
#######################################################


#Minorly looking at taxa from starting samples
#######################################################
fungdata_start <- subset_samples(fungdata_decontam, Time_point=="pre-inoculation" | Time_point=="post-inoculation")
  fungdata__start <- filter_taxa(fungdata_start, function(x) sum(x)!=0, prune=T)

#aglommerating data to Family level for plotting
fungdata_startGlomSpecies <- glom_tax(fungdata__start, "Species")
fungdata_startGlomSpeciesNorm <- transform_sample_counts(fungdata_startGlomSpecies, function(x) x / sum(x))

summary_table(fungdata_startGlomSpeciesNorm)
df <- as.data.frame(summary_table_fungdata_startGlomSpeciesNorm[8:19])
df <- sapply(df, function(x) as.numeric(as.character(x)))

mean <- rowMeans(df)
max <- apply(df, 1, function(x) max(x))
prevalence <- colSums(apply(df, 1, function(x) x>0))

summary_table_fungdata_startGlomSpeciesNorm <- cbind(summary_table_fungdata_startGlomSpeciesNorm, mean, max, prevalence)

write.csv(summary_table_fungdata_startGlomSpeciesNorm, "fungdata_start.csv", quote=F)

#######################################################
