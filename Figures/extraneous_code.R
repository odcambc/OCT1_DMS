# Extra code.

oct1_combined_scores %>% filter(mutation_type != "X") %>% 
  group_by(pos) %>% 
  mutate(mean_sm73 = mean(-1*SM73_1_score, na.rm = TRUE),
         mean_gfp = mean(GFP_score, na.rm = TRUE)) %>%
  ggplot() +
  geom_density(aes(x = mean_sm73)) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() 


oct1_combined_scores %>% filter(mutation_type != "X") %>% 
  group_by(pos) %>% 
  mutate(mean_sm73 = mean(-1*SM73_1_score, na.rm = TRUE),
         mean_gfp = mean(GFP_score, na.rm = TRUE)) %>%
  ggplot() +
  geom_density(aes(x = mean_gfp)) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() 

oct1_combined_scores %>% filter(mutation_type != "X") %>% 
  group_by(pos) %>% 
  mutate(mean_sm73 = mean(-1*SM73_1_score, na.rm = TRUE),
         mean_gfp = mean(GFP_score, na.rm = TRUE)) %>%
  ggplot() +
  geom_point(aes(x = mean_gfp, y = mean_sm73)) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() 

linreg_reg_spec <- 
  linear_reg() %>% 
  set_engine("lm")

linreg_reg_spec

linreg_reg_fit <- linreg_reg_spec %>% fit(mean_sm73 ~ mean_gfp, data = oct1_combined_scores %>%
                                            filter(mutation_type != "X" &
                                                     abs(SM73_1_SE) < 0.4 &
                                                     abs(GFP_score) < 0.4) %>% 
                                            group_by(pos) %>% 
                                            mutate(mean_sm73 = mean(-1*SM73_1_score, na.rm = TRUE),
                                                   mean_gfp = mean(GFP_score, na.rm = TRUE)))

intercept = linreg_reg_fit$fit$coefficients[[1]]
slope = linreg_reg_fit$fit$coefficients[[2]]

oct1_combined_scores %>%
  filter(mutation_type != "X" &
           abs(SM73_1_SE) < 0.4 &
           abs(GFP_score) < 0.4) %>% 
  group_by(pos) %>% 
  mutate(regressed_function = -1*SM73_1_score - (slope*GFP_score + intercept)) %>%
  ggplot() +
  geom_density(aes(x = regressed_function, color = mutation_type), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + 
  theme_classic() 


oct1_combined_scores %>% filter(mutation_type != "X") %>%
  ungroup() %>%
  mutate(regressed_function = SM73_1_score - (slope*GFP_score + intercept)) %>%
  ggplot() +
  geom_point(aes(x = regressed_function, y = GFP_score ,color = mutation_type), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + 
  theme_classic() 


# Fitting

function_em_fit <- normalmixEM(
  oct1_combined_scores %>% filter(!is.na(SM73_1_score)) %>% .$SM73_1_score
)

folding_em_fit <- normalmixEM(
  oct1_combined_scores %>% filter(!is.na(GFP_score)) %>% .$GFP_score
)


function_em_fit <- normalmixEM(
  oct1_combined_scores %>% filter(!is.na(SM73_1_score) & mutation_type == "S") %>% .$SM73_1_score
)

folding_em_fit <- normalmixEM(
  oct1_combined_scores %>% filter(!is.na(GFP_score) & mutation_type == "S") %>% .$GFP_score
)

function_em_fit$sigma
function_em_fit$mu

function_em_fit$mu[[1]] + function_em_fit$sigma[[1]]
function_em_fit$mu[[1]] - function_em_fit$sigma[[1]]

function_em_fit$mu[[2]] + function_em_fit$sigma[[2]]
function_em_fit$mu[[2]] - function_em_fit$sigma[[2]]


plot(function_em_fit, whichplots = 2)

plot(folding_em_fit)


oct1_combined_scores %>% filter(mutation_type != "X") %>% ungroup() %>% mutate(regressed_function = -1*SM73_1_score - (0.6321804*GFP_score + 0.1797040)) %>%
  ggplot() +
  geom_point(aes(y = regressed_function, x = GFP_score)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() 

oct1_combined_scores %>% filter(mutation_type != "X") %>% 
  ungroup() %>% 
  mutate(regressed_function = -1*SM73_1_score - (0.6321804*GFP_score + 0.1797040)) %>%
  group_by(pos) %>%
  summarize(mean = mean(regressed_function, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(y = mean, x = pos)) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic() 


#Define a new parameter to remove the conditional dependence of observed
#function on stability. We expect that the observed functional score (SM73_1_score)
#will depend on the abundance levels. First, just do a linear regression to correct
#to first order.

linreg_reg_spec <- 
  linear_reg() %>% 
  set_engine("lm")

linreg_reg_spec

linreg_reg_fit <- linreg_reg_spec %>% fit(-1*SM73_1_score ~ GFP_score, data = oct1_combined_scores %>% 
                                            filter(!mutation_type %in% c("X") & abs(SM73_1_SE) < 0.4 & abs(GFP_score) < 0.4) %>% 
                                            ungroup())

intercept = linreg_reg_fit$fit$coefficients[[1]]
slope = linreg_reg_fit$fit$coefficients[[2]]

oct1_combined_scores %>% filter(mutation_type != "X") %>%
  ungroup() %>%
  mutate(regressed_function = -1*SM73_1_score - (slope*GFP_score + intercept)) %>%
  ggplot() +
  geom_density(aes(x = regressed_function, color = mutation_type), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + 
  theme_classic() 




#Generate some synthetic data for the importance scores:
#  sums of absolute values of normally distributed variables? scaled somehow?

total_scores = c()
total_scores_2 = c()

for (i in 1:10000){
  sum_scores = sum(abs(rnorm(20)))
  total_scores = c(total_scores, sum_scores)
}

for (i in 1:10000){
  sum_scores_2 = sum(abs(rnorm(500)))
  total_scores_2 = c(total_scores_2, sum_scores_2)
}

hist(total_scores)
hist(total_scores_2)




#Investigate some of the statistical properties of these scores.

ggplot(oct1_pos_scores) +
  geom_point(aes(x = low_GFP, y = low_function)) +
  xlab("Stability") +
  ylab("Function") +
  theme_classic()


ggplot(oct1_combined_scores %>%
         group_by(pos) %>%
         summarize(all_func = sum(abs(SM73_1_score)),
                   all_gfp = sum(abs(GFP_score)))) +
  geom_point(aes(x = all_gfp, y = all_func)) +
  geom_point(data = oct1_pos_scores, aes(x = low_GFP, y = low_function), color = "grey") +
  xlab("Stability") +
  ylab("Function") +
  theme_classic()


ggplot() +
  geom_density(data = oct1_combined_scores,
               aes(x = GFP_score), size = 2) +
  geom_density(data = oct1_combined_scores %>% 
                 filter(GFP_score <= abundance_cutoff),
               aes(x = GFP_score), color = "red") +
  geom_density(data = oct1_combined_scores %>% 
                 filter(GFP_score > abundance_cutoff),
               aes(x = GFP_score), color = "blue") +
  theme_classic()

ggplot() +
  geom_density(data = oct1_combined_scores,
               aes(x = SM73_1_score), size = 2) +
  geom_density(data = oct1_combined_scores %>% 
                 filter(GFP_score <= abundance_cutoff),
               aes(x = SM73_1_score), color = "red", linetype = 2) +
  geom_density(data = oct1_combined_scores %>% 
                 filter(GFP_score > abundance_cutoff),
               aes(x = SM73_1_score), color = "blue", linetype = 2) +
  theme_classic()

oct1_combined_scores %>% filter(GFP_score <= abundance_cutoff) %>%
  summarize(summed_gfp = sum(GFP_score * GFP_score),
            var_gfp = var(GFP_score),
            n_gfp = n()) %>%
  ggplot() +
  geom_density(aes(x = var_gfp)) +
  geom_density(aes(x = summed_gfp/n_gfp))

oct1_combined_scores %>% filter(!is.na(GFP_score))%>%
  group_by(pos) %>% 
  summarize(summed_gfp = sum(GFP_score * GFP_score),
            var_gfp = var(GFP_score),
            n_gfp = n()) %>%
  ggplot() +
  geom_density(aes(x = var_gfp)) +
  geom_density(aes(x = summed_gfp/n_gfp))




summed_fit <- fitdistr(
  oct1_combined_scores %>% filter(!is.na(GFP_score))%>%
    group_by(pos) %>% 
    summarize(summed_gfp = sum(GFP_score * GFP_score)) %>%
    .$summed_gfp,
  densfun = "chi-squared",
  start = list(df = 2,
               ncp = 0)
)


low_gfp_fit <- fitdistr(
  oct1_pos_scores %>% filter(low_GFP > 0) %>% .$low_GFP,
  densfun = "chi-squared",
  start = list(df = 2,
               ncp = 0),
  method = "L-BFGS-B" 
)

low_func_fit <- fitdistr(
  oct1_pos_scores$low_function,
  densfun = "chi-squared",
  start = list(df = 2,
               ncp = 0)
)

plot(low_func_fit)

ggplot(oct1_pos_scores) +
  geom_histogram(aes(x = low_GFP), trim = FALSE, fill = "red") +
  geom_histogram(aes(x = low_function), trim = FALSE, fill = "blue", alpha = 0.5) +
  geom_function(fun = dexp, args = list(rate = low_func_fit$estimate[["rate"]]
  ),
  colour = "blue") +
  theme_classic()


dens2 <- ggplot(oct1_combined_scores %>% filter(mutation_type != "X"),
                aes(x = GFP_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  geom_function(fun = dnorm, args = list(mean = folding_synonymous_fit$estimate[["mean"]],
                                         sd = folding_synonymous_fit$estimate[["sd"]]
  ),
  colour = "red") +
  geom_vline(xintercept = folding_synonymous_fit$estimate[["mean"]] +
               2 * folding_synonymous_fit$estimate[["sd"]], linetype = 2) +
  geom_vline(xintercept = folding_synonymous_fit$estimate[["mean"]] -
               2 * folding_synonymous_fit$estimate[["sd"]], linetype = 2) +
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()


oct1_combined_scores %>%
  filter(SM73_1_score >= sm73_cutoff ) %>%
  group_by(pos) %>%
  summarize(abundance = sum(abs(GFP_score)),
            func = sum(abs(SM73_1_score))) %>%
  ggplot() +
  geom_density(aes(x = abundance), color = "red") +
  geom_density(aes(x = func), color = "blue") +
  theme_classic()

oct1_combined_scores %>%
  filter(SM73_1_score >= sm73_cutoff ) %>%
  ggplot() +
  geom_density(aes(x = GFP_score), color = "red") +
  geom_density(aes(x = SM73_1_score), color = "blue") +
  theme_classic()


library(fitdistrplus)

descdist(oct1_pos_scores$low_GFP)
descdist(oct1_pos_scores$low_function)

descdist(oct1_pos_scores %>% filter(low_GFP > 0) %>% .$low_GFP)
descdist(oct1_pos_scores %>% filter(low_function > 0) %>% .$low_function)


plotdist(oct1_pos_scores %>% filter(low_function > 0) %>% .$low_function, "gamma", para = list(shape = 1.24686656, rate = 0.37761061))
plotdist(oct1_pos_scores$low_function, "gamma", para = list(shape = 1.24686656, rate = 0.37761061))

ggplot(oct1_pos_scores) +
  geom_density(aes(x = low_GFP), color = "red") +
  geom_density(aes(x = low_function), color = "blue") +
  theme_classic()

low_gfp_fit <- fitdistr(
  oct1_pos_scores$low_GFP,
  densfun = "gamma"
)

low_func_fit <- fitdistr(
  oct1_pos_scores %>% filter(low_function > 0) %>% .$low_function,
  densfun = "gamma"
)

plot(low_func_fit)




# Fig 3

oct1_pos_cors <- oct1_combined_scores %>% filter(mutation_type == "M") %>%
  group_by(pos) %>%
  mutate(pocket = ifelse(pos %in% pocket_residues, 1, 0),
         pocket_A = ifelse(pos %in% pocket_residues_A, 1, 0),
         pocket_B = ifelse(pos %in% pocket_residues_B, 1, 0),
         pocket_C = ifelse(pos %in% pocket_residues_C, 1, 0),
         pocket_D = ifelse(pos %in% pocket_residues_D, 1, 0)) %>%
  summarize(correlation = cor(-1*SM73_1_score, GFP_score, use="na.or.complete", method = "spearman"),
            correlation_sm73_GFP = cor(-1*SM73_1_score, GFP_score, use="na.or.complete", method = "spearman"),
            correlation_pol_sm73 = cor(polarity_distance, -1*SM73_1_score, use="na.or.complete", method = "spearman"),
            correlation_pol_GFP = cor(polarity_distance, GFP_score, use="na.or.complete", method = "spearman"),
            correlation_h_sm73 = cor(hydropathy_distance, -1*SM73_1_score, use="na.or.complete", method = "spearman"),
            correlation_h_GFP = cor(hydropathy_distance, GFP_score, use="na.or.complete", method = "spearman"),
            correlation_h2_sm73 = cor(hydropathy_2_distance, -1*SM73_1_score, use="na.or.complete", method = "spearman"),
            correlation_h2_GFP = cor(hydropathy_2_distance, GFP_score, use="na.or.complete", method = "spearman"),
            correlation_v_sm73 = cor(volume_distance, -1*SM73_1_score, use="na.or.complete", method = "spearman"),
            correlation_v_GFP = cor(volume_distance, GFP_score, use="na.or.complete", method = "spearman"),
            correlation_b_sm73 = cor(buriability_distance, -1*SM73_1_score, use="na.or.complete", method = "spearman"),
            correlation_b_GFP = cor(buriability_distance, GFP_score, use="na.or.complete", method = "spearman"),
            correlation_i_sm73 = cor(isoelectric_distance, -1*SM73_1_score, use="na.or.complete", method = "spearman"),
            correlation_i_GFP = cor(isoelectric_distance, GFP_score, use="na.or.complete", method = "spearman"),
            correlation_s_sm73 = cor(stability_distance, -1*SM73_1_score, use="na.or.complete", method = "spearman"),
            correlation_s_GFP = cor(stability_distance, GFP_score, use="na.or.complete", method = "spearman"),
            correlation_d_sm73 = cor(donor_distance, -1*SM73_1_score, use="na.or.complete", method = "spearman"),
            correlation_d_GFP = cor(donor_distance, GFP_score, use="na.or.complete", method = "spearman"),
            pocket = first(pocket),
            pocket_A = first(pocket_A),
            pocket_B = first(pocket_B),
            pocket_C = first(pocket_C),
            pocket_D = first(pocket_D)                                                    
  ) %>% ungroup()

# Make a long version for plotting

oct1_pos_cors_long <- oct1_pos_cors %>% pivot_longer(cols = matches("correlation_.*_sm73"),
                                                     values_to = "corr",
                                                     names_pattern = "correlation_(.*)_sm73",
                                                     names_to = c("property")) %>%
  select(c(pos, pocket,
           pocket_A, pocket_B,
           pocket_C, pocket_D,
           property, corr))

#Now plot.


#polarity
plot_pocket_pol_cor <- ggplot(data = oct1_pos_cors) +
  geom_density(aes(x = correlation_pol_sm73, color = factor(pocket))) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_classic()

print(plot_pocket_pol_cor) 

# hydrophobicity
plot_pocket_h_cor <- ggplot(data = oct1_pos_cors) +
  geom_density(aes(x = correlation_h_sm73, color = factor(pocket))) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_classic()

print(plot_pocket_h_cor) 

# volume
plot_pocket_v_cor <- ggplot(data = oct1_pos_cors) +
  geom_density(aes(x = correlation_v_sm73, color = factor(pocket))) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_classic()

print(plot_pocket_v_cor) 

# buriability
plot_pocket_b_cor <- ggplot(data = oct1_pos_cors) +
  geom_density(aes(x = correlation_b_sm73, color = factor(pocket))) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_classic()

print(plot_pocket_b_cor) 

# isoelectric point
plot_pocket_i_cor <- ggplot(data = oct1_pos_cors) +
  geom_density(aes(x = correlation_i_sm73, color = factor(pocket))) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_classic()

print(plot_pocket_i_cor) 

# stability score
plot_pocket_s_cor <- ggplot(data = oct1_pos_cors) +
  geom_density(aes(x = correlation_s_sm73, color = factor(pocket))) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_classic()

print(plot_pocket_s_cor) 

# h-bond donation
plot_pocket_d_cor <- ggplot(data = oct1_pos_cors) +
  geom_density(aes(x = correlation_d_sm73, color = factor(pocket))) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_classic()

print(plot_pocket_d_cor) 



#Plot pocket vs non-pocket correlations between different properties.

plot_correlations <- ggplot(data = oct1_pos_cors_long) +
  geom_violin(aes(y = corr, x = factor(property), fill = factor(pocket))) +
  geom_boxplot(aes(y = corr, x = factor(property), color = factor(pocket)), width = 0.5) +
  scale_fill_manual(values = c("grey", "blue"),
                    name="Residue class",
                    labels=c("Nonbinding", "Binding")) +
  scale_color_manual(values = c("grey", "blue"),
                     name="Residue class",
                     labels=c("Nonbinding", "Binding")) +
  xlab("Residue property") +
  ylab("Positional correlation with function") +
  theme_classic()

print(plot_correlations)

ggsave("output/Figure_3_correlations.png", width = 5, height = 3, plot_correlations)
ggsave("output/Figure_3_correlations.pdf", width = 5, height = 3, plot_correlations)

plot_correlations <- ggplot(data = oct1_pos_cors_long) +
  geom_violin(aes(y = corr, x = factor(property), fill = factor(pocket))) +
  scale_fill_manual(values = c("grey", "blue"),
                    name="Residue class",
                    labels=c("Nonbinding", "Binding")) +
  xlab("Residue property") +
  ylab("Positional correlation with function") +
  theme_classic()

ggsave("output/Figure_3_correlations_nobox.png", width = 5, height = 3, plot_correlations)
ggsave("output/Figure_3_correlations_nobox.pdf", width = 5, height = 3, plot_correlations)



# Quantify movement between states
b_apo_xyz <- b_apo$xyz[atom.select(b_apo, "calpha")$xyz]
d_mpp_xyz <- d_mpp$xyz[atom.select(d_mpp, "calpha")$xyz]

d_mpp_xyz_aln <- fit.xyz(fixed = b_apo_xyz, mobile = d_mpp_xyz)

plot.dmat(dist.xyz(b_apo_xyz, d_mpp_xyz_aln))

bd_distance <- dist.xyz(b_apo_xyz, d_mpp_xyz_aln)

plot.dmat(bd_distance[TM2, TM7])

bd_df <- as.data.frame.table(bd_distance)

bd_df[1:2] <- lapply(bd_df[1:2], as.numeric)

bd_difference_plot <- ggplot(data = bd_df) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 50, 100)) +
  geom_segment(x = 145, xend = 145, y = 0, yend = 535, linetype = 2) + 
  geom_segment(x = 46, xend = 46, y = 0, yend = 535, linetype = 2) + 
  geom_segment(x = 172, xend = 172, y = 0, yend = 535, linetype = 2) + 
  geom_segment(x = 462, xend = 462, y = 0, yend = 535, linetype = 2) + 
  geom_segment(x = 490, xend = 490, y = 0, yend = 535, linetype = 2) + 
  geom_segment(y = 145, yend = 145, x = 0, xend = 535, linetype = 2) + 
  geom_segment(y = 46, yend = 46, x = 0, xend = 535, linetype = 2) + 
  geom_segment(y = 172, yend = 172, x = 0, xend = 535, linetype = 2) + 
  geom_segment(y = 462, yend = 462, x = 0, xend = 535, linetype = 2) + 
  geom_segment(y = 490, yend = 490, x = 0, xend = 535, linetype = 2) + 
  coord_fixed() +
  guides(fill = guide_legend(title = "Distance")) +
  theme_classic()

ggsave("output/Distance_change.png", width = 6, height = 6, bd_difference_plot)
ggsave("output/Distance_change.pdf", width = 6, height = 6, bd_difference_plot)


ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM7))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 50, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM2))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 50, 100)) +
  coord_fixed() +
  theme_classic()



#----
TM1 <- c(16:43)
TM3 <- c(175:195)
TM4 <- c(199:227)
TM5 <- c(230:257)
TM6 <- c(261:281)
TM7 <- c(341:371)
TM8 <- c(375:400)
TM9 <- c(402:422)
TM10 <- c(428:456)
TM11 <- c(461:490)
TM12 <- c(492:512)

ICD <- c(282:340)



ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM1))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM2))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM3))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM4))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM5))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM6))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM7))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM8))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM9))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM10))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM11))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM12))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()





ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM1))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM2))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM3))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM4))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM5))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM6))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM7))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM8))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM9))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM10))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM11))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()

ggplot(data = bd_df %>% filter(Var1 %in% c(TM7) & Var2 %in% c(TM12))) + 
  geom_raster(aes(x = Var1, y = Var2, fill = Freq)) +
  xlab("Residue") + 
  ylab("Residue") +
  scale_fill_fermenter(type = "seq",
                       direction = 1,
                       palette = 1,
                       breaks = c(0, 10, 25, 40, 65, 100)) +
  coord_fixed() +
  theme_classic()



sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM1)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM1)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM2)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM2)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM3)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM3)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM4)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM4)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM5)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM5)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM6)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM6)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM7)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM7)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM8)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM8)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM9)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM9)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM10)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM10)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM11)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM11)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM12)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM1) & Var2 %in% c(TM12)))[1]


#---


sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM1)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM1)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM2)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM2)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM3)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM3)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM4)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM4)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM5)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM5)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM6)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM6)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM7)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM7)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM8)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM8)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM9)) %>% .$Freq)  / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM9)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM10)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM10)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM11)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM11)))[1]

sum(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM12)) %>% .$Freq) / 
  dim(bd_df %>% filter(Var1 %in% c(TM2) & Var2 %in% c(TM12)))[1]
