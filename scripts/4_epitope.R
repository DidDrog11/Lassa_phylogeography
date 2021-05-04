source(here::here("scripts", "library_api.R"))

#Potential epitope for BP is 0.35
#Potential epitope for BP_2 is 0.5
#Potential epitope for CF is 0.978
#Potential epitope for E is 1
#Potential epitope for KS is 0.978
#Potential epitope for KT is 1.036
#Potential epitope for KT is 0.725

gp_epitope <- read_csv(here("data", "epitope_prediction", "epitope_prediction_gp.csv")) %>%
  dplyr::select(Position_BP_2, Residue_BP_2, all_of(starts_with("Score"))) %>%
  rename("Position" = Position_BP_2,
         "Residue" = Residue_BP_2) %>%
  pivot_longer(names_to = "Algorithm",
               cols = c("Score_BP_2", "Score_BP",
                        "Score_CF", "Score_E",
                        "Score_KS", "Score_KT",
                        "Score_P"),
               values_to = "Score",
               names_prefix = "Score_") %>%
  mutate(Algorithm = case_when(Algorithm == "BP_2" ~ "Bepipred-2",
                               Algorithm == "BP" ~ "Bepipred",
                               Algorithm == "CF" ~ "Chou-Fasman",
                               Algorithm == "E" ~ "Emini",
                               Algorithm == "KS" ~ "Karplus-Schulz",
                               Algorithm == "KT" ~ "Kolasker-Tongaonkar",
                               Algorithm == "P" ~ "Parker",
                               TRUE ~ "Missing"))

gp_epitope %<>%
  group_by(Algorithm) %>%
  mutate(Alpha = case_when(Algorithm == "Bepipred" & Score < 0.35 ~ 0.2,
                           Algorithm == "Bepipred-2" & Score < 0.5 ~ 0.2,
                           Algorithm == "Chou-Fasman" & Score < 0.978 ~ 0.2,
                           Algorithm == "Emini" & Score < 1 ~ 0.2,
                           Algorithm == "Karplus-Schulz" & Score < 0.978 ~ 0.2,
                           Algorithm == "Kolasker-Tongaonkar" & Score < 1.036 ~ 0.2,
                           Algorithm == "Parker" & Score < 0.725 ~ 0.2,
                           is.na(Score) ~ 0,
                           TRUE ~ 1),
         Proportional_score = case_when(Alpha == 0.2 | Alpha == 0 ~ 0,
                                        TRUE ~ Score),
         Proportional_score = ntile(Proportional_score, 100))

gp_epitope %>%
  ggplot() +
  geom_col(aes(x = Position, y = Score, alpha = Alpha, fill = Proportional_score)) +
  facet_wrap(~ Algorithm, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_c(option = "plasma", direction = 1, limits = c(49, 100)) +
  labs(title = "Epitope prediction: Glycoprotein")
  
ggsave(filename = here("outputs", "glycoprotein_epitope.png"), plot = last_plot())

np_epitope <- read_csv(here("data", "epitope_prediction", "epitope_prediction_np.csv")) %>%
  dplyr::select(Position_BP_2, Residue_BP_2, all_of(starts_with("Score"))) %>%
  rename("Position" = Position_BP_2,
         "Residue" = Residue_BP_2) %>%
  pivot_longer(names_to = "Algorithm",
               cols = c("Score_BP_2", "Score_BP",
                        "Score_CF", "Score_E",
                        "Score_KS", "Score_KT",
                        "Score_P"),
               values_to = "Score",
               names_prefix = "Score_") %>%
  mutate(Algorithm = case_when(Algorithm == "BP_2" ~ "Bepipred-2",
                               Algorithm == "BP" ~ "Bepipred",
                               Algorithm == "CF" ~ "Chou-Fasman",
                               Algorithm == "E" ~ "Emini",
                               Algorithm == "KS" ~ "Karplus-Schulz",
                               Algorithm == "KT" ~ "Kolasker-Tongaonkar",
                               Algorithm == "P" ~ "Parker",
                               TRUE ~ "Missing"))


np_epitope %<>%
  group_by(Algorithm) %>%
  mutate(Alpha = case_when(Algorithm == "Bepipred" & Score < 0.35 ~ 0.2,
                           Algorithm == "Bepipred-2" & Score < 0.5 ~ 0.2,
                           Algorithm == "Chou-Fasman" & Score < 0.978 ~ 0.2,
                           Algorithm == "Emini" & Score < 1 ~ 0.2,
                           Algorithm == "Karplus-Schulz" & Score < 0.978 ~ 0.2,
                           Algorithm == "Kolasker-Tongaonkar" & Score < 1.036 ~ 0.2,
                           Algorithm == "Parker" & Score < 0.725 ~ 0.2,
                           is.na(Score) ~ 0,
                           TRUE ~ 1),
         Proportional_score = case_when(Alpha == 0.2 | Alpha == 0 ~ 0,
                                        TRUE ~ Score),
         Proportional_score = ntile(Proportional_score, 100))

np_epitope %>%
  ggplot() +
  geom_col(aes(x = Position, y = Score, alpha = Alpha, fill = Proportional_score)) +
  facet_wrap(~ Algorithm, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_c(option = "plasma", direction = 1, limits = c(30, 100)) +
  labs(title = "Epitope prediction: Nucleoprotein")

ggsave(filename = here("outputs", "nucleoprotein_epitope.png"), plot = last_plot())
