# Initial data processing for OCT1 DMS data.
# Reads the Enrich2 output, generates a dataframe with features annotated
# Outputs csvs for further interactive plotting and analysis

library(dplyr)
library(tidyverse)

source("dms_analysis_utilities.R")

# Enrich2 score files - input
oct1_scores_file ="../data/Enrich2/Oct1_full/tsv/OCT1_full_exp/main_identifiers_scores.tsv"
oct1_scores_shared_file ="../data/Enrich2/Oct1_full/tsv/OCT1_full_exp/main_identifiers_scores_shared.tsv"

# Enrich2 score files - output
oct1_combined_scores_file ="../data/oct1_combined_scores.csv"

# Read in the OCT1 scores

oct1_variantscore_colnames <- c("hgvs", "SE", "epsilon", "score")
oct1_variantscore_coltypes <- c("character", "numeric", "numeric", "numeric")

oct1_scores <- read_delim(oct1_scores_file, col_names = oct1_variantscore_colnames, skip=3)
oct1_scores <- process_hgvs_df(oct1_scores)

oct1_variantscore_colnames <- c("hgvs", "SM73_0_SE", "SM73_0_epsilon", "SM73_0_score", "SM73_1_SE", "SM73_1_epsilon", "SM73_1_score", "GFP_SE", "GFP_epsilon", "GFP_score")
oct1_variantscore_coltypes <- c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")

oct1_shared_colnames <- c("hgvs", "SM73_0_SE_R1", "SM73_0_score_R1",
                          "SM73_0_SE_R2", "SM73_0_score_R2",
                          "SM73_0_SE_R3", "SM73_0_score_R3",
                          "SM73_1_SE_R1", "SM73_1_score_R1",
                          "SM73_1_SE_R2", "SM73_1_score_R2",
                          "SM73_1_SE_R3", "SM73_1_score_R3",
                          "GFP_SE_R1", "GFP_score_R1",
                          "GFP_SE_R2", "GFP_score_R2",
                          "GFP_SE_R3", "GFP_score_R3")
oct1_shared_coltypes <- c("character", "numeric", "numeric",
                          "numeric", "numeric",
                          "numeric", "numeric",
                          "numeric", "numeric",
                          "numeric", "numeric",
                          "numeric", "numeric",
                          "numeric", "numeric",
                          "numeric", "numeric",
                          "numeric", "numeric")

oct1_scores <- read_delim(oct1_scores_file, col_names = oct1_variantscore_colnames, skip=3)
oct1_scores <- process_hgvs_df(oct1_scores)
oct1_scores$is.wt = oct1_scores$mutation_type == "S"

oct1_scores_shared <- read_delim(oct1_scores_shared_file, col_names = oct1_shared_colnames, skip=5)

# Enrich2 does not calculate scores if a variant in a replicate goes to 0 counts at some point.
# This calculates a score for those cases by averaging the scores for the other replicates.
oct1_scores_shared <- oct1_scores_shared %>%  rowwise() %>%
  mutate(SM73_0_score = mean(c_across(starts_with("SM73_0_score")), na.rm=TRUE),
         SM73_1_score = mean(c_across(starts_with("SM73_1_score")), na.rm=TRUE),
         GFP_score = mean(c_across(starts_with("GFP_score")), na.rm=TRUE))

oct1_scores_shared <- process_hgvs_df(oct1_scores_shared %>% select(c(hgvs, SM73_0_score, SM73_1_score, GFP_score)))
oct1_scores_shared$is.wt = oct1_scores_shared$mutation_type == "S"

oct1_combined_scores <-
  oct1_scores %>% left_join(oct1_scores_shared %>% 
                              select(c('SM73_0_score', 'SM73_1_score', 'GFP_score', 'hgvs')), by = 'hgvs') %>%
  mutate(SM73_0_score = coalesce(SM73_0_score.x, SM73_0_score.y),
         SM73_1_score = coalesce(SM73_1_score.x, SM73_1_score.y),
         GFP_score = coalesce(GFP_score.x, GFP_score.y)) %>%
  select(-one_of(c("SM73_0_score.x", 'SM73_0_score.y',
                   "SM73_1_score.x", 'SM73_1_score.y',
                   "GFP_score.x", 'GFP_score.y')))


# Since variant processing removes entries that are unobserved, for plotting and
# analysis purposes we add them as missing data to this dataframe, now. We take
# them from our designed variant file.
oct1_designed_variants_file ="../designed_variants/oct1_variants.csv"
oct1_variants <- read_delim(oct1_designed_variants_file)

oct1_combined_scores <- oct1_combined_scores %>% full_join(oct1_variants %>% select(pos, mutation_type,
                                                            mutation, length, hgvs),
                                   by = c("hgvs" = "hgvs",
                                          "pos" = "pos",
                                          "mutation_type" = "mutation_type",
                                          "len" = "length",
                                          "variants" = "mutation"))

# Generate a few amino acid property indices.
# Polarity: GRAR740102 Polarity (Grantham, 1974)
#polarity_index <- aa.index[111][[1]]$I
polarity_index <- aa.index$GRAR740102$I

#ZIMJ680104	Isoelectric point (Zimmerman et al., 1968)	
isoelectric_index <- aa.index$ZIMJ680104$I

#ROSM880102	Side chain hydropathy, corrected for solvation (Roseman, 1988)	
hydropathy_index <- aa.index$ROSM880102$I

#ENGD860101 Hydrophobicity index (Engelman et al., 1986)
hydropathy_index_2 <- aa.index$ENGD860101$I

#PONJ960101 Average volumes of residues (Pontius et al., 1996)
volume_index <- aa.index$PONJ960101$I

#ZHOH040103 Buriability (Zhou-Zhou, 2004)
buriability_index <- aa.index$ZHOH040103$I

#ZHOH040102 The relative stability scale extracted from mutation experiments (Zhou-Zhou,  2004)
stability_index <- aa.index$ZHOH040102$I

#GEOR030108 Linker propensity from helical (annotated by DSSP) dataset (George-Heringa,  2003)
linker_index_1 <- aa.index$GEOR030108$I

#GEOR030105 Linker propensity from small dataset (linker length is less than six  residues) (George-Heringa, 2003)
linker_index_2 <- aa.index$GEOR030105$I

#GEOR030106 Linker propensity from medium dataset (linker length is between six and 14  residues) (George-Heringa, 2003)
linker_index_3 <- aa.index$GEOR030106$I

#GEOR030107 Linker propensity from long dataset (linker length is greater than 14  residues) (George-Heringa, 2003)
linker_index_4 <- aa.index$GEOR030107$I

#FAUJ880109 Number of hydrogen bond donors (Fauchere et al., 1988)
donor_index <- aa.index$FAUJ880109$I

oct1_combined_scores <- oct1_combined_scores %>%
  left_join(oct1_combined_scores |> 
              group_by(pos) |>
              filter(pos > 0) |>
              mutate(wt_pos = substr(hgvs, 4, 4))) %>% ungroup()

oct1_combined_scores <- oct1_combined_scores %>% 
  left_join( oct1_combined_scores |> 
               filter(mutation_type %in% c("S", "M")) |> 
               mutate(polarity_distance = polarity_index[wt_pos] - polarity_index[variants],
                      hydropathy_distance = hydropathy_index[wt_pos] - hydropathy_index[variants],
                      hydropathy_2_distance = hydropathy_index_2[wt_pos] - hydropathy_index_2[variants],
                      volume_distance = volume_index[wt_pos] - volume_index[variants],
                      buriability_distance = buriability_index[wt_pos] - buriability_index[variants],
                      isoelectric_distance = isoelectric_index[wt_pos] - isoelectric_index[variants],
                      stability_distance = stability_index[wt_pos] - stability_index[variants],
                      linker_1_distance = linker_index_1[wt_pos] - linker_index_1[variants],
                      linker_2_distance = linker_index_2[wt_pos] - linker_index_2[variants],
                      linker_3_distance = linker_index_3[wt_pos] - linker_index_3[variants],
                      linker_4_distance = linker_index_4[wt_pos] - linker_index_4[variants],
                      donor_distance = donor_index[wt_pos] - donor_index[variants]))

write_csv(oct1_combined_scores, oct1_combined_scores_file)