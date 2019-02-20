#' Further improved solution to summarize OpenBis.se gDDA workflow results on peptide and protein level per MS injection.
#' Purpose: easy accessibility and visualization of DDA results before library creation.
#' Major change: Filter proteins and peptides to ca. 1 % FDR as assessed by Mayu each
#' MH 20.02.2019

## Set parameters
mayu_protein_fdr_target = 0.01 # Modify to override the FDR target and corresponding i-probability cutoff used based on the Mayu table
mayu_peptide_fdr_target = 0.01 # Modify to override the FDR target and corresponding i-probability cutoff used based on the Mayu table

## install + load required packages
if(is.null(require(data.table))){
  install.packages("data.table")
}
if(is.null(require(ggplot2))){
  install.packages("ggplot2")
}
if(is.null(require(pheatmap))){
  install.packages("pheatmap")
}
if(is.null(require(corrplot))){
  install.packages("corrplot")
}
library("data.table")
library("ggplot2")
library("pheatmap")
library("corrplot")

## load data
peptides = fread("peptides.tsv")
mayu = fread("mayuout_main_1.07.csv")

# Comment: NO USAGE OF protein prophet results, the Mayu FDR estimation is used

## report mayu fdr's of the current analyte set
# reporting
min_iprobability = min(peptides$probability_ip)
message("Minimum iProbability: ", min_iprobability)

# Mayu tests a few probability cutoffs in a fixed spacing (parameter of spacing set when run)
# Thus the cutoffs tested do not precisely match that of the i-prophet fdr model
# Therefore, report FDRs and analyte ids at the closest cutoff
# --> get closest value in mayu table and report FDRs estimated by Mayu
mayu_ip_cutoffs = mayu$`IP/PPs`
mayu_ip_cutoffs_delta = abs(mayu_ip_cutoffs - min_iprobability)
closest_cutoff = mayu_ip_cutoffs[which(mayu_ip_cutoffs_delta == min(mayu_ip_cutoffs_delta))]
message("Mayu-estimated FDR of the unfiltered peptides.tsv results at the closest evaluated iProbability cutoff ", closest_cutoff, ":")
message(paste(c("iProbability", "pepFDR", "protFDR", "target_pepID", "target_protID"), mayu[`IP/PPs` == closest_cutoff, .(`IP/PPs`, pepFDR, protFDR, target_pepID, target_protID)], sep = "="))
message("Summary of the unfiltered peptides.tsv results with minimal iProbability ", min_iprobability)
decoy_idx = grep("DECOY", peptides$protein)
target_idx = grep("DECOY", peptides$protein, invert = T)
message("Peptides up to target pepFDR for all proteins passing target protFDR are maintained")
message("N Target proteins: ", length(unique(peptides[target_idx]$protein)))
message("N Decoy proteins: ", length(unique(peptides[decoy_idx]$protein)))
message("N Target peptides: ", length(unique(peptides[target_idx]$peptide)))
message("N Decoy peptides: ", length(unique(peptides[decoy_idx]$peptide)))

# Now extract the iProbability cutoffs to achieve peptide and protein FDR target values as defined above (default = 0.01)
mayu[, dpepFDR:=abs(pepFDR-mayu_peptide_fdr_target)]
mayu[, dprotFDR:=abs(protFDR-mayu_protein_fdr_target)]
mayu_peptide_fdr_iprob_cutoff= mayu[dpepFDR == min(dpepFDR), unique(`IP/PPs`)]
mayu_protein_fdr_iprob_cutoff= mayu[dprotFDR == min(dprotFDR), unique(`IP/PPs`)]

# Filter peptides to FDR target
peptides_unfiltered = copy(peptides)
peptides_1pcpepfdr = peptides[probability_ip >= mayu_peptide_fdr_iprob_cutoff]

message("Summary metrics of the Mayu-pepFDR-filtered peptides.tsv table, filtered with iProbability cutoff ", mayu_peptide_fdr_iprob_cutoff, " are:")
message(paste(c("iProbability", "pepFDR", "protFDR", "target_pepID", "target_protID"), mayu[`IP/PPs` == mayu_peptide_fdr_iprob_cutoff, .(`IP/PPs`, pepFDR, protFDR, target_pepID, target_protID)], sep = "="))

# Filter the FDR-controlled peptide table to proteins passing the protein FDR target
mayu_protein_fdr_conforming_proteins = peptides[probability_ip >= mayu_protein_fdr_iprob_cutoff, unique(protein)]
peptides_1pcpepfdr_1pcprotfdr = peptides_1pcpepfdr[protein %in% mayu_protein_fdr_conforming_proteins]

message("Summary of the Mayu-pepFDR-protFDR filtered results, peptides filtered with iProbability cutoff ",
        mayu_peptide_fdr_iprob_cutoff, " and protein identifiers filtered with cutoff ", mayu_protein_fdr_iprob_cutoff, ":")

decoy_idx = grep("DECOY", peptides_1pcpepfdr_1pcprotfdr$protein)
target_idx = grep("DECOY", peptides_1pcpepfdr_1pcprotfdr$protein, invert = T)
message("Peptides were filtered to fdr target ", mayu_peptide_fdr_target, "based on Mayu result table")
message("Proteins were filtered to fdr target ", mayu_protein_fdr_target, "based on Mayu result table")
message("Peptides up to target pepFDR for all proteins passing target protFDR are maintained")
message("N Target proteins: ", length(unique(peptides_1pcpepfdr_1pcprotfdr[target_idx]$protein)))
message("N Decoy proteins: ", length(unique(peptides_1pcpepfdr_1pcprotfdr[decoy_idx]$protein)))
message("N Target peptides: ", length(unique(peptides_1pcpepfdr_1pcprotfdr[target_idx]$peptide)))
message("N Decoy peptides: ", length(unique(peptides_1pcpepfdr_1pcprotfdr[decoy_idx]$peptide)))

# Convert the filtered peptide level info to wide format & report
long <- copy(peptides_1pcpepfdr_1pcprotfdr)
long[, filename:= sapply(long$spectrum, function(x){strsplit(x, split = "\\.")[[1]][1]})]
# long = long[filename %in% unique(filename)[1:2]] # testing plotting settings for two-run-only datasets
long[, max_iprob:= max(probability_ip), by = peptide]
long[, n_spectra:= .N, by = c("filename","peptide")]
long[, spectralcounter:=1]
long[, n_parent_proteins:= length(unique(protein)), peptide]
long[, parent_proteins:= paste(unique(protein), collapse = ","), peptide]
long[, n_peptides:=length(unique(peptide)), parent_proteins]
long[, peptides:=paste(unique(peptide), collapse = ","), parent_proteins]

# summarize per peptide in wide format
peptides_wide <- dcast(data = long, formula = peptide + protein + n_parent_proteins + parent_proteins + max_iprob ~ filename, value.var = "spectralcounter", fun.aggregate = sum)
setorder(peptides_wide, peptide)
setorder(peptides_wide, parent_proteins)

# summarize spectral counts per protein in wide format
proteins_wide <- dcast(data = long, formula = parent_proteins + n_parent_proteins + n_peptides + peptides ~ filename, value.var = "spectralcounter", fun.aggregate = sum, fill = 0)
setorder(proteins_wide, parent_proteins)

# write out tables
write.table(peptides_wide, file = paste0("peptides_wide_pepFDR",mayu_peptide_fdr_target, "_protFDR", mayu_protein_fdr_target, ".tsv"), sep = "\t", row.names = F, quote = F)
write.table(proteins_wide, file = paste0("proteins_wide_speccount_pepFDR",mayu_peptide_fdr_target, "_protFDR", mayu_protein_fdr_target, ".tsv"), sep = "\t", row.names = F, quote = F)

#' Visualization:
# consider decoys separately/label for plotting
long[, is_decoy:=0]
long$is_decoy[grep("DECOY", long$parent_proteins)] = 1
long$is_decoy = factor(long$is_decoy)

# MZ maps
long[, mz:=(precursor_neutral_mass+assumed_charge)/assumed_charge]
ggplot(long[is_decoy == 0], aes(retention_time_sec, mz)) +
  geom_point(pch = ".", aes(color = is_decoy)) +
  facet_wrap(~filename) +
  ggtitle("MZ-RT maps of target peptide - spectrum matches")
ggsave("visualizeObisDDAresults_mzmaps.pdf", device = "pdf", height = 4+length(unique(long$filename))/2, width = 4+length(unique(long$filename)))

# peptides per protein global
ggplot(unique(long[, .(parent_proteins, n_peptides, is_decoy)]), aes(n_peptides, fill = is_decoy)) +
  geom_bar(position = position_dodge(), alpha = 0.6) +
  scale_fill_manual(values = c("black", "red")) +
  coord_cartesian(xlim = c(0,10)) +
  scale_x_continuous(breaks = seq(1:10))+
  ylab("Number of proteins (groups)") + xlab ("n_peptides per protein (unique naked sequences)") +
  ggtitle("Number of peptides (unique naked sequences) per protein (group)")
ggsave("visualizeObisDDAresults_pepsperprot.pdf", device = "pdf", height = 7, width = 4+length(unique(long$filename))/2)


# Peptides per run
long[, peps_per_run:=length(unique(peptide)), .(filename, is_decoy)]
ggplot(unique(long[, .(filename, peps_per_run, is_decoy)]), aes(filename, y = peps_per_run, fill = is_decoy)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.6) +
  scale_fill_manual(values = c("black", "red")) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of peptides identified per run")
ggsave("visualizeObisDDAresults_pepsperrun.pdf", device = "pdf", height = 7, width = 4+length(unique(long$filename))/2) 

# Proteins per run
long[, prots_per_run:=length(unique(parent_proteins)), .(filename, is_decoy)]

ggplot(unique(long[, .(filename, prots_per_run, is_decoy)]), aes(filename, y = prots_per_run, fill = is_decoy)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.6) +
  scale_fill_manual(values = c("black", "red")) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of protein groups identified per run, >= 1 pep per prot")
ggsave("visualizeObisDDAresults_protsperrun_min1pep.pdf", device = "pdf", height = 7, width = 4+length(unique(long$filename))/2) 

# Proteins per run with 2 or more peptides
long_min2 = long[n_peptides >= 2]
long_min2[, prots_per_run:=length(unique(parent_proteins)), .(filename, is_decoy)]
ggplot(unique(long_min2[, .(filename, prots_per_run, is_decoy)]), aes(filename, y = prots_per_run, fill = is_decoy)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.6) +
  scale_fill_manual(values = c("black", "red")) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of protein groups identified per run, >= 2 pep per prot")
ggsave("visualizeObisDDAresults_protsperrun_min2pep.pdf", device = "pdf", height = 7, width = 4+length(unique(long$filename))/2) 

# Spectral counts per protein distributions
proteins_long = melt(proteins_wide, id.vars = names(proteins_wide)[1:4])
proteins_long$is_decoy = 0
proteins_long[grep("DECOY", parent_proteins)]$is_decoy = 1
proteins_long[, is_decoy:=as.factor(is_decoy)]

ggplot(proteins_long, aes(x = variable, y = value, fill = is_decoy, colour = is_decoy)) +
  geom_violin(position = "dodge", alpha = 0.5) +
  geom_point(size = 1, position = position_jitterdodge(), alpha = 0.5) +
  scale_y_log10() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values = c("black", "red")) +
  ylab("Spectral counts per protein (log scaled axis)") + xlab("filename") +
  ggtitle("Distributions of spectral counts per protein (violin plot)")
ggsave("visualizeObisDDAresults_speccountsperprotperrun_violin.pdf", device = "pdf", height = 7, width = 4+length(unique(long$filename))/1.5) 

# Spectral counts per protein boxplots
ggplot(proteins_long, aes(x = variable, y = value, fill = is_decoy, colour = is_decoy)) +
  geom_boxplot(position = "dodge", alpha = 0.5) +
  geom_point(size = 1, position = position_jitterdodge(), alpha = 0.5) +
  scale_y_log10() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values = c("black", "red")) +
  ylab("Spectral counts per protein (log scaled axis)") + xlab("filename") +
  ggtitle("Distributions of spectral counts per protein (box plot)")
ggsave("visualizeObisDDAresults_speccountsperprotperrun_boxplots.pdf", device = "pdf", height = 7, width = 4+length(unique(long$filename))) 

# Correlation analysis, sample-sample
#######################
protmatrix = as.matrix(proteins_wide[, 5:ncol(proteins_wide)])
rownames(protmatrix) = proteins_wide$parent_proteins
log10_protmatrix = log(protmatrix)
log10_protmatrix[is.infinite(log10_protmatrix)] = 0

# sample-sample correlation heatmap
ss = cor(log10_protmatrix) # log transform? ss = cor(log(protmatrix)[(!is.infinite(log(protmatrix)))])
corrplot.mixed(ss, upper = "color", lower = "number", tl.pos = "lt")
pdf("visualizeObisDDAresults_corr_sample-sample.pdf", width = 2+ncol(protmatrix), height = 2+ncol(protmatrix))
corrplot.mixed(ss, upper = "color", lower = "number", tl.pos = "lt")
dev.off()

# Sample-sample and protein-protein visualized in heatmap (log10-transformed counts)
pheatmap(log10_protmatrix)
dev.off()
pdf("visualizeObisDDAresults_corr_heatmap_reordered.pdf", width = 4+ncol(protmatrix)/2, height = 4+nrow(protmatrix)/8)
pheatmap(log10_protmatrix, main = "DDA spectral counts heatmap, log10-transformed counts")
dev.off()


# For comparison, also convert and report unfiltered tables
# Convert the filtered peptide level info to wide format & report
long.u <- copy(peptides_unfiltered)
long.u[, filename:= sapply(long$spectrum, function(x){strsplit(x, split = "\\.")[[1]][1]})]
long.u[, max_iprob:= max(probability_ip), by = peptide]
long.u[, n_spectra:= .N, by = c("filename","peptide")]
long.u[, spectralcounter:=1]
long.u[, n_parent_proteins:= length(unique(protein)), peptide]
long.u[, parent_proteins:= paste(unique(protein), collapse = ","), peptide]
long.u[, n_peptides:=length(unique(peptide)), parent_proteins]
long.u[, peptides:=paste(unique(peptide), collapse = ","), parent_proteins]

# summarize per peptide in wide format
peptides_wide.u <- dcast(data = long.u, formula = peptide + protein + n_parent_proteins + parent_proteins + max_iprob ~ filename, value.var = "spectralcounter", fun.aggregate = sum)
setorder(peptides_wide.u, peptide)
setorder(peptides_wide.u, parent_proteins)

# summarize spectral counts per protein in wide format
proteins_wide.u <- dcast(data = long.u, formula = parent_proteins + n_parent_proteins + n_peptides + peptides ~ filename, value.var = "spectralcounter", fun.aggregate = sum, fill = 0)
setorder(proteins_wide.u, parent_proteins)

# write out unfiltered tables
write.table(peptides_wide.u, file = "unfiltered_peptides_wide.tsv", sep = "\t", row.names = F, quote = F)
write.table(proteins_wide.u, file = "unfiltered_proteins_wide_speccount.tsv", sep = "\t", row.names = F, quote = F)

