source(here::here("scripts", "library_api.r"))

search_date <- "2021-09-24"

all_sequences <- read_rds(here("data", paste0("all_sequences_", search_date, ".rds"))) %>%
  mutate(region = snakecase::to_sentence_case(region),
         region = recode(region,
                         "Ebonyi state" = "Ebonyi",
                         "Eddo" = "Edo",
                         "Edo state" = "Edo",
                         "Enugu state" = "Enugu",
                         "Fct" = "Federal capital territory",
                         "Nasarawa state" = "Nasarawa",
                         "Ondo state" = "Ondo",
                         "Plateau state" = "Plateau",
                         "Songo hospital sebgwema" = "Sebgwema",
                         "Tararba" = "Taraba"),
         accession_lassa = as.character(accession_lassa),
         accession_lassa = str_remove(accession_lassa, "\\.1"))

no_country <- all_sequences %>%
  filter(country == "NULL")

reference_sequences <- all_sequences %>%
  filter(accession_lassa %in% c("NC_004297", "NC_004296")) %>%
  mutate(segment = c("L","S"))
# These are the genbank reference sequences

all_sequences <- all_sequences %>%
  filter(!accession_lassa %in% no_country$accession_lassa)

np_ref_length <- reference_sequences %>%
  filter(!is.na(np_length)) %>%
  .$np_length
gpc_ref_length <- reference_sequences %>%
  filter(!is.na(gpc_length)) %>%
  .$gpc_length
l_ref_length <- reference_sequences %>%
  filter(!is.na(l_length)) %>%
  .$l_length
z_ref_length <- reference_sequences %>%
  filter(!is.na(z_length)) %>%
  .$z_length
l_seq_ref_length <- reference_sequences %>%
  filter(segment == "L") %>%
  .$nuc_length
s_seq_ref_length <- reference_sequences %>%
  filter(segment == "S") %>%
  .$nuc_length

np_sequences <- all_sequences %>%
  filter(!is.na(np_length)) %>%
  select(accession_lassa, defined_lassa, host, country, region, year, gene, np_seq, np_length) %>%
  mutate(np_length_exp = ifelse(np_length > np_ref_length*0.9 & np_length < np_ref_length*1.1, T, F)) %>%
  filter(country != "NULL") %>% # remove entries without attributed countries
  filter(np_length_exp == T) # remove entries that are >10% different in length from the reference

gpc_sequences <- all_sequences %>%
  filter(!is.na(gpc_length)) %>%
  select(accession_lassa, defined_lassa, host, country, region, year, gene, gpc_seq, gpc_length) %>%
  mutate(gpc_length_exp = ifelse(gpc_length > gpc_ref_length*0.9 & gpc_length < gpc_ref_length*1.1, T, F)) %>%
  filter(country != "NULL")%>%
  filter(gpc_length_exp == T)

l_sequences <- all_sequences %>%
  filter(!is.na(l_length)) %>%
  select(accession_lassa, defined_lassa, host, country, region, year, gene, l_seq, l_length) %>%
  mutate(l_length_exp = ifelse(l_length > l_ref_length*0.9 & l_length < l_ref_length*1.1, T, F)) %>%
  filter(country != "NULL") %>%
  filter(l_length_exp == T)

z_sequences <- all_sequences %>%
  filter(!is.na(z_length)) %>%
  select(accession_lassa, defined_lassa, host, country, region, year, gene, z_seq, z_length) %>%
  mutate(z_length_exp = ifelse(z_length > z_ref_length*0.9 & z_length < z_ref_length*1.1, T, F)) %>%
  filter(country != "NULL") %>%
  filter(z_length_exp == T)

l_segment_sequence <- all_sequences %>%
  mutate(l_segment = ifelse(nuc_length > l_seq_ref_length*0.9 & nuc_length < l_seq_ref_length*1.1, T, F)) %>%
  filter(!is.na(l_segment)) %>%
  select(accession_lassa, defined_lassa, host, country, region, year, nuc_seq, nuc_length, name_nuc_seq, l_segment) %>%
  filter(country != "NULL") %>%
  filter(l_segment == T)

s_segment_sequence <- all_sequences %>%
  mutate(s_segment = ifelse(nuc_length > s_seq_ref_length*0.9 & nuc_length < s_seq_ref_length*1.1, T, F)) %>%
  filter(!is.na(s_segment)) %>%
  select(accession_lassa, defined_lassa, host, country, region, year, nuc_seq, nuc_length, name_nuc_seq, s_segment) %>%
  filter(country != "NULL") %>%
  filter(s_segment == T)

incomplete_sequences <- all_sequences %>%
  filter(country != "NULL") %>%
  filter(!accession_lassa %in% c(l_segment_sequence$accession_lassa, s_segment_sequence$accession_lassa))

np_gpc_sequences <- full_join(np_sequences,
                              gpc_sequences,
                              by = c("accession_lassa", "defined_lassa", "host", "country", "region", "year", "gene"))

l_z_sequences <- full_join(l_sequences,
                           z_sequences,
                           by = c("accession_lassa", "defined_lassa", "host", "country", "region", "year", "gene"))

write_rds(np_gpc_sequences, here("cleaned_data", paste("cleaned_np_gpc_", search_date, ".rds", sep = "")))
write_rds(l_z_sequences, here("cleaned_data", paste("cleaned_l_z_", search_date, ".rds", sep = "")))
write_rds(s_segment_sequence, here("cleaned_data", paste("cleaned_s_segment_", search_date, ".rds", sep = "")))
write_rds(l_segment_sequence, here("cleaned_data", paste("cleaned_l_segment_", search_date, ".rds", sep = "")))


# The list of sequences we have aligned previously
previously_aligned <- read_csv(here("cleaned_data", "aligned_sequences.csv"))

new_s <- s_segment_sequence %>%
  filter(!accession_lassa %in% previously_aligned$accession_s) %>%
  dplyr::select(accession_lassa, nuc_seq)
seqRFLP::dataframe2fas(new_s, file = here("cleaned_data", "new_s_sequences.fasta"))

new_l <- l_segment_sequence %>%
  filter(!accession_lassa %in% previously_aligned$accession_l) %>%
  dplyr::select(accession_lassa, nuc_seq)
seqRFLP::dataframe2fas(new_l, file = here("cleaned_data", "new_l_sequences.fasta"))

new_np <- np_sequences %>%
  filter(!accession_lassa %in% previously_aligned$accession_np) %>%
  dplyr::select(accession_lassa, np_seq)
seqRFLP::dataframe2fas(new_np, file = here("cleaned_data", "new_np_sequences.fasta"))

new_gpc <- gpc_sequences %>%
  filter(!accession_lassa %in% previously_aligned$accession_gp) %>%
  dplyr::select(accession_lassa, gpc_seq)
seqRFLP::dataframe2fas(new_gpc, file = here("cleaned_data", "new_gpc_sequences.fasta"))
