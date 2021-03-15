all_sequences <- read_rds(here("data","all_sequences_2021-03-15.rds")) %>%
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
                         "Tararba" = "Taraba"))

reference_sequences <- all_sequences %>%
  filter(accession_lassa %in% c("NC_004297", "NC_004296"))
# These are the genbank reference sequences

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

np_sequences <- all_sequences %>%
  filter(!is.na(np_length)) %>%
  select(accession_lassa, defined_lassa, host, country, region, year, gene, np_seq, np_length) %>%
  mutate(np_length_exp = ifelse(np_length > np_ref_length*0.9 & np_length < np_ref_length*1.1, T, F))

gpc_sequences <- all_sequences %>%
  filter(!is.na(gpc_length)) %>%
  select(accession_lassa, defined_lassa, host, country, region, year, gene, gpc_seq, gpc_length) %>%
  mutate(gpc_length_exp = ifelse(gpc_length > gpc_ref_length*0.9 & gpc_length < gpc_ref_length*1.1, T, F))

l_sequences <- all_sequences %>%
  filter(!is.na(l_length)) %>%
  select(accession_lassa, defined_lassa, host, country, region, year, gene, l_seq, l_length) %>%
  mutate(l_length_exp = ifelse(l_length > l_ref_length*0.9 & l_length < l_ref_length*1.1, T, F))

z_sequences <- all_sequences %>%
  filter(!is.na(z_length)) %>%
  select(accession_lassa, defined_lassa, host, country, region, year, gene, z_seq, z_length) %>%
  mutate(z_length_exp = ifelse(z_length > z_ref_length*0.9 & z_length < z_ref_length*1.1, T, F))

np_gpc_sequences <- full_join(np_sequences %>%
                                  filter(np_length_exp == T),
                                gpc_sequences %>%
                                  filter(gpc_length_exp == T),
                                by = c("accession_lassa", "defined_lassa", "host", "country", "region", "year", "gene")) %>%
  filter(country != "NULL")

write_rds(np_gpc_sequences, here("data", paste("cleaned_sequences_", search_date, ".rds", sep = "")))
