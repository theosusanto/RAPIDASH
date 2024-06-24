library(riboWaltz)
df=create_annotation(gtfpath = "data/gencode.v43.primary_assembly.basic.annotation.gtf")
bam_list <- bamtolist("data/", annotation = df, refseq_sep = "|")

# Try filter only read lengths that satisfy a periodicity threshold
# (i.e. 70% of the reads fall in one of the three reading frames along the CDS)
#filtered_bam_list <- length_filter(data = bam_list,
#                                   length_filter_mode = "periodicity",
#                                   periodicity_threshold = 70)
# This doesn't seem to work? It threw away 99.9% of the reads. The reads kept have much longer length.
# Apparently it is because I don't have good periodicity to begin with. RNAse or CHX might be an issue.

offsets = psite(data = bam_list)
samples = unique(offsets$sample)
for (i in 1:length(samples)) {
  sam_offs = subset(offsets, sample == samples[i], 
                    select = c("length", "corrected_offset_from_3"))
  colnames(sam_offs) = c("length", "p_offset")
  write.table(sam_offs, paste("output/Try2", samples[i], "_p-site-offsets", sep=""),sep = "\t", 
              row.names = FALSE)
}

reads_psite_list <- psite_info(bam_list, offsets)
codon_coverage_example <- codon_coverage(reads_psite_list, df, psite = FALSE)
cds_coverage_example <- cds_coverage(reads_psite_list, df)

Nterm_1_psite_region <- region_psite(reads_psite_list, df, sample = "Nterm_dep_rep1_RFP_Aligned.toTranscriptome.out")
Nterm_1_psite_region[["plot"]]
Nterm_1_frames_stratified <- frame_psite_length(reads_psite_list, sample = "Nterm_dep_rep1_RFP_Aligned.toTranscriptome.out",
                                                region = "all", cl = 90)
Nterm_1_frames_stratified[["plot"]]

comparison_list <- list()6
comparison_list[["subsample_28nt"]] <- reads_psite_list[["Nterm_dep_rep1_RFP_Aligned.toTranscriptome.out"]][length == 28]
comparison_list[["whole_sample"]] <- reads_psite_list[["Nterm_dep_rep1_RFP_Aligned.toTranscriptome.out"]]
sample_list <- list("Only_28" = c("subsample_28nt"),
                    "All" = c("whole_sample"))

Nterm_1_metaheatmap <- metaheatmap_psite(comparison_list, df, sample = sample_list,
                                         utr5l = 20, cdsl = 40, utr3l = 20, log_colour = F, 
                                         plot_title = "Comparison metaheatmap")
Nterm_1_metaheatmap[["plot"]]
Nterm_1_metaprofile <- metaprofile_psite(reads_psite_list, df, sample = "Nterm_dep_rep1_RFP_Aligned.toTranscriptome.out",
                                         utr5l = 20, cdsl = 40, utr3l = 20,
                                         plot_title = "Nterm_1.transcript")
Nterm_1_metaprofile[["plot_Nterm_dep_rep1_RFP_Aligned.toTranscriptome.out"]]

#To plot read counts against position given a particular gene in a given sample:

transcript_id = "ENST00000519479.2"
reads_psite_sample = "wt_dep_rep1_RFP_Aligned.toTranscriptome.out"

psite_counts <- reads_psite_list$reads_psite_sample[
  reads_psite_list$reads_psite_sample$transcript == transcript_id,]$psite

library(tidyr)
library(dplyr)

# Convert counts table to data frame, change th
counts <- as.data.frame(table(psite_counts))
colnames(counts) <- c("Position", "ReadCount")

# Generate a sequence of numbers from 0 or min. number to the max. number
counts$Position <- as.numeric(as.character(counts$Position))
number_range <- seq(pmin(min(counts$Position), 0), max(counts$Position))

# Complete the data frame with missing numbers and set their counts to 0
counts_complete <- counts %>%
  complete(Position = number_range) %>%
  replace_na(list(ReadCount = 0))

ggplot(counts_complete, aes(x = Position, y = ReadCount)) +
  geom_line() +
  labs(x = "Position", y = "Read Count", title = "Read Counts vs. Positions") +
  theme_minimal()