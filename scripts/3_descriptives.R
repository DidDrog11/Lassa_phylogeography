source(here::here("scripts", "library_api.r"))

search_date <- "2021-09-24"

np_gpc <- read_rds(here("cleaned_data", paste0("cleaned_np_gpc_", search_date, ".rds")))
l_z <- read_rds(here("cleaned_data", paste0("cleaned_l_z_", search_date, ".rds")))
s_seq <- read_rds(here("cleaned_data",  paste0("cleaned_s_segment_", search_date, ".rds")))
l_seq <- read_rds(here("cleaned_data", paste0("cleaned_l_segment_", search_date, ".rds")))

included_records <- unique(c(np_gpc$accession_lassa,
                             l_z$accession_lassa,
                             s_seq$accession_lassa,
                             l_seq$accession_lassa))

included_phylogenetics <- read_csv(here("data", "included_accession.csv"))
included_phylogenetics_september <- bind_rows(tibble(segment = "L", accession_lassa = included_phylogenetics$L_accession),
                                              tibble(segment = "S", accession_lassa = included_phylogenetics$S_accession))

included_records <- tibble(included_phylogenetics_september) %>%
  left_join(., tibble(bind_rows(np_gpc, l_z, s_seq, l_seq) %>%
                        select(accession_lassa, name_nuc_seq, host, country, region, year) %>%
                        mutate(accession_lassa = as.character(accession_lassa)) %>%
                        distinct(accession_lassa, .keep_all = T)),
            by = "accession_lassa")

all_sequences <- read_rds(here("data", paste0("all_sequences_", search_date, ".rds"))) %>%
  select(accession_lassa, name_nuc_seq, host, country, region, year) %>%
  mutate(accession_lassa = as.character(accession_lassa),
         country = case_when(country == "Cote d'Ivoire" ~ "Côte d'Ivoire",
                             TRUE ~ country)) %>%
  distinct(accession_lassa, .keep_all = T)

prop.table(table(all_sequences$country))

regional_data <- all_sequences %>%
  filter(!is.na(region))

prop.table(table(regional_data$country))

# Year of collection ----------------------------------------------------------
table(all_sequences$year)

annual_sequences <- all_sequences %>%
  group_by(country, year) %>%
  summarise(n = n()) %>%
  mutate(country = case_when(country == "NULL" ~ "Missing",
                             TRUE ~ country)) %>%
  ggplot() +
  geom_col(aes(x = year, y = n, fill = country)) +
  theme_minimal() +
  labs(x = "Year",
       y = "Number of sequences") +
  scale_fill_manual(name = "Country (n)",
                    labels = c("Benin (32)",
                               "Côte d'Ivoire (58)",
                               "Ghana (4)",
                               "Guinea (442)",
                               "Liberia (83)",
                               "Mali (25)",
                               "Missing (134)",
                               "Nigeria (1217)",
                               "Sierra Leone (297)",
                               "Togo (6)"),
                    values = c("#871a10", "#e6853b",
                               "#000000", "#be2a2c",
                               "#0a2864", "#50ab4a",
                               "#89CFF0", "#398053",
                               "#2f71c1", "#f7d147")) +
  theme(legend.position = "bottom")

# Records mapping ---------------------------------------------------------

# Attempt to enrich locations from citations
# Andersen et al. collected samples from Kenema and Irrua
andersen <- read_xlsx(here("data", "andersen_2015_accession.xlsx")) %>%
  select(country = "Country", s_accession = "S Acc. #", l_accession = "L Acc. #") %>%
  mutate(s_accession = case_when(s_accession == "N/A" ~ as.character(NA),
                                 TRUE ~ s_accession),
         l_accession = case_when(l_accession == "N/A" ~ as.character(NA),
                                 TRUE ~ l_accession)) %>%
  filter(country %in% c("Sierra Leone", "Nigeria")) %>%
  pivot_longer(cols = contains("accession"), values_to = "accession_lassa") %>%
  drop_na(accession_lassa) %>%
  mutate(accession_lassa = paste0(accession_lassa, ".1"),
         region = case_when(country == "Nigeria" ~ "Edo",
                            country == "Sierra Leone" ~ "Kenema")) %>%
  mutate(host = case_when(accession_lassa %in% c("KM822128.1", "KM822127.1") ~ "Homo sapiens")) %>%
  select(-name)

# Yadouleton et al. collected samples from an outbreak in Benin
yadouleton <- read_xlsx(here("data", "yadouleton_2020_accession.xlsx")) %>%
  filter(!is.na(region)) %>%
  select(country, region, gp_accession = "N° GP", np_accession = "N° NP", l_accession = "N° L") %>%
  mutate(gp_accession = case_when(gp_accession == "NA" ~ as.character(NA),
                                 TRUE ~ gp_accession),
         np_accession = case_when(np_accession == "NA" ~ as.character(NA),
                                 TRUE ~ np_accession),
         l_accession = case_when(l_accession == "NA" ~ as.character(NA),
                                  TRUE ~ l_accession)) %>%
  pivot_longer(cols = contains("accession"), values_to = "accession_lassa") %>%
  drop_na(accession_lassa) %>%
  mutate(accession_lassa = paste0(accession_lassa, ".1")) %>%
  select(-name)

# Mateo et al. produced a sequence from a fatality from Bangolo District, Cote d'Ivoire
mateo <- all_sequences %>%
  filter(accession_lassa %in% c("MK978784.1", "MK978785.1")) %>%
  mutate(region = "Bangolo") %>%
  select(accession_lassa, country, region)

# Baumann and Mandy samples were from Bamako
baumann <- all_sequences %>%
  filter(accession_lassa %in% c("MH473587.1", "MH473586.1")) %>%
  mutate(region = "Bamako") %>%
  select(accession_lassa, country, region)

# Welch et al. samples were from several regions in Liberia
welch <- read_xlsx(here("data", "welch_2019_accession.xlsx")) %>%
  select(l_accession, s_accession, country = Country, region = `County or prefecture`) %>%
  pivot_longer(cols = c("l_accession", "s_accession"), names_to = "segment", values_to = "accession_lassa") %>%
  mutate(accession_lassa = paste0(accession_lassa, ".1"),
         host = "Homo Sapiens")

# Gunther et al. sequence is from Bantou, Guinea
gunther <- all_sequences %>%
  filter(accession_lassa %in% c("MK044799.1")) %>%
  mutate(region = "Bantou") %>%
  select(accession_lassa, country, region)

# Siddle et al. 2018
siddle <- read_xlsx(here("data", "siddle_2018_accession.xlsx")) %>%
  mutate(region = str_split(ID, pattern = "-", simplify = TRUE)[,2],
         country = "Nigeria") %>%
  pivot_longer(cols = contains("Accession"), values_to = "accession_lassa") %>%
  drop_na(accession_lassa) %>%
  select(accession_lassa, country, region) %>%
  mutate(accession_lassa = paste0(accession_lassa, ".1"))

# Metsky at al. 2019
metsky <- read_xlsx(here("data", "metsky_2019_accession.xlsx")) %>%
  mutate(country = "Nigeria") %>%
  mutate(accession_lassa = paste0(accession_lassa, ".1"))

# Rossi et al.
rossi <- all_sequences %>%
  filter(accession_lassa %in% c("KU978807.1", "KU978808.1", "KU978811.1", "KU978812.1")) %>%
  select(accession_lassa, country) %>%
  mutate(region = case_when(country == "Guinea" ~ "Faranah",
                            country == "Nigeria" ~ "Osun"))

# Omilabu et al. 2016
omilabu <- all_sequences %>%
  filter(accession_lassa %in% c("DQ010031.1", "DQ010030.1", "KJ944261.1", "KJ944260.1", "KJ944259.1")) %>%
  select(accession_lassa, country) %>%
  mutate(region = case_when(accession_lassa %in% c("DQ010031.1", "DQ010030.1") ~ "Edo",
                            accession_lassa == "KJ944261.1" ~ "Lagos",
                            accession_lassa == "KJ944260.1" ~ "Ebonyi",
                            accession_lassa == "KJ944259.1" ~ "Enugu"),
         host = "Homo sapiens")

# Vieth et al. 2007 with missing countries were from Sierra Leone
vieth <- all_sequences %>%
  filter(str_starts(accession_lassa, "AY3")) %>%
  mutate(country = "Sierra Leone",
         host = "Homo sapiens") %>%
  select(accession_lassa, country, region, host) %>%
  bind_rows(all_sequences %>%
              filter(accession_lassa %in% c("AY870334.1", "AY693640.1", "AY693639.1",
                                            "AY693638.1", "AY693637.1")) %>%
              mutate(host = "Homo sapiens") %>%
            select(accession_lassa, country, region, host))

all_sequences$host[all_sequences$accession_lassa %in% c("AY179172.1", "AY179173.1")] <- "Homo sapiens"

# Ehichioya et al. 2011 and 2012
ehichioya <- read_xlsx(here("data", "ehichioya_2011_accession.xlsx")) %>%
  pivot_longer(cols = contains("accession"), values_to = "accession_lassa") %>%
  drop_na(accession_lassa) %>%
  select(accession_lassa, country, region) %>%
  mutate(accession_lassa = paste0(accession_lassa, ".1"))

# Safronetz et al. 2010
safronetz <- all_sequences %>%
  filter(str_starts(name_nuc_seq, "Soromba|Ouoma|Komina|Bamba")) %>%
  mutate(region = "Sikasso") %>%
  select(accession_lassa, country, region)

# Atkin et al. 2009
atkin <- all_sequences %>%
  filter(accession_lassa == "FJ824031.1") %>%
  mutate(region = "Sikasso") %>%
  select(accession_lassa, country, region)

# Kouadio et al. 2015
kouadio <- all_sequences %>%
  filter(accession_lassa %in% c("LN823985.1", "LN823984.1", "LN823983.1", "LN823982.1")) %>%
  mutate(region = "Komborodougou") %>%
  select(accession_lassa, country, region)

# Leski et al. 2014
leski <- read_xlsx(here("data", "leski_2014_accession.xlsx")) %>%
  pivot_longer(cols = contains("accession"), values_to = "accession_lassa") %>%
  drop_na(accession_lassa) %>%
  select(accession_lassa, region) %>%
  mutate(country = "Sierra Leone") %>%
  mutate(accession_lassa = paste0(accession_lassa, ".1"))

# Bonney et al. 2013
bonney <- all_sequences %>%
  filter(accession_lassa %in% c("KF425246.1")) %>%
  mutate(country = "Liberia",
         region = "Lofa")

# Bowen et al. 2000
bowen <- read_xlsx(here("data", "bowen_2000_accession.xlsx")) %>%
  select(accession_lassa, country, region, host = species)  %>%
  mutate(accession_lassa = paste0(accession_lassa, ".1"))

# Lecompte et al. 2006
lecompte <- all_sequences %>%
  filter(str_starts(accession_lassa, "DQ832")) %>%
  mutate(host = "Mastomys natalensis") %>%
  select(accession_lassa, country, region, host)

found_locations <- bind_rows(andersen, yadouleton, mateo, baumann, welch,
                             gunther, siddle, metsky, rossi, omilabu,
                             vieth, ehichioya, safronetz, atkin,
                             kouadio, leski, bonney, bowen, lecompte) %>%
  select(accession_lassa, country, region, host)

enriched <- all_sequences %>%
  mutate(country = case_when(country == "NULL" ~ as.character(NA),
                             TRUE ~ country),
         host = case_when(host == "NULL" ~ as.character(NA),
                          TRUE ~ host),
         region = case_when(region == "NULL" ~ as.character(NA),
                            is.na(region) ~ as.character(NA),
                            TRUE ~ region)) %>%
  left_join(found_locations, by = c("accession_lassa")) %>%
  group_by(accession_lassa) %>%
  mutate(country = coalesce(country.y, country.x),
         region = coalesce(region.y, region.x),
         host = coalesce(host.y, host.x)) %>%
  arrange(accession_lassa, region, country, host, year) %>%
  slice(1) %>%
  select(accession_lassa, name_nuc_seq, host, year, country, region) %>%
  distinct()

# Ondo itself identifies as a region in Osun state it should fix it by specifying the state
enriched$region[str_detect(enriched$region, "Ondo|ONDO")] <- "Ondo state"
enriched$region[str_detect(enriched$region, "Nasarawa State")] <- "Nasarawa"
enriched$region[str_detect(enriched$region, "Abuja")] <- "Federal Capital Territory"

geocoded <- enriched %>%
  mutate(region = ifelse(region == "Unknown", NA, region)) %>%
  unite("location", c(region, country), sep = ", ", na.rm = T, remove = F)

coordinate <- tibble(location = unique(geocoded$location)) %>%
  bind_cols(geocode(unique(geocoded$location))) %>% # This will pull coordinates for the location (built from region and country) for some reason Togo is listed as being in the US
  filter(location != "Togo" | location == "NULL")

country_coordinate <- na.omit(unique(geocoded$country[geocoded$country != "NULL"])) %>%
  bind_cols(geocode(na.omit(unique(geocoded$country[geocoded$country != "NULL"])))) %>%
  rename("country" = 1,
         "country_lon" = "lon",
         "country_lat" = "lat") %>%
  mutate(country_lon = ifelse(country == "Togo", 1.217, country_lon),
         country_lat = ifelse(country == "Togo", 6.133, country_lat))

geocoded <- geocoded %>%
  left_join(., coordinate, by = "location") %>%
  mutate(lon = case_when(location == "Togo" ~ 1.217,
                         TRUE ~ lon),
         lat = case_when(location == "Togo" ~ 6.133,
                         TRUE ~ lat))

missing_coords = tibble(region = c("Silimi", "Safrani", "Khoria", "Brissa", "Denguedou", "Songo Hospital, Sebgwema", "Bantou", "Tanganya", "Gbetaya", "Odo-akaba", "Damania",
                                   "Nyandeyama", "Worogui", "Jirandogo"),
                        lon = c(-10.65, -10.738367, -10.893033, -10.688767, -10.448611, -10.948721, -10.583333, -10.972783, -11.040150, 2.601698, -10.8667, -11.05305, 2.67307, -0.345833),
                        lat = c(9.966667,  10.057817, 9.942267, 10.216833, 8.495556, 8.004604, 10.066667, 10.000400, 9.841017, 8.774286, 9.8, 8.46485, 8.884305, 8.347167))

fix_coords <- geocoded %>%
  filter(!is.na(region)) %>%
  mutate(matched_lon = case_when(lon %in% country_coordinate$country_lon ~ FALSE,
                                 TRUE ~ TRUE),
         matched_lat = case_when(lat %in% country_coordinate$country_lat ~ FALSE,
                                 TRUE ~ TRUE)) %>%
  filter(matched_lon == FALSE | matched_lat == FALSE) %>%
  left_join(missing_coords, by = c("region")) %>%
  select(accession_lassa, name_nuc_seq, host, location, country, region, year, lon = lon.y, lat = lat.y)

geocoded <- geocoded %>%
  filter(!accession_lassa %in% fix_coords$accession_lassa) %>%
  bind_rows(fix_coords)


country_geocode <- geocoded %>%
  left_join(., country_coordinate, by = "country")

country_summary <- country_geocode %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  drop_na(country) %>%
  mutate(country = case_when(country == "Cote d'Ivoire" ~ "Côte d'Ivoire",
                             TRUE ~ country))

w_africa <- read_rds(here("data_download", "level_0_admin.rds")) %>%
  filter(!GID_0 %in% c("CMR", "DZA", "ESH", "MAR", "TCD"))

map_seqs <- w_africa %>%
  left_join(., country_summary,
            by = c("NAME_0" = "country")) %>%
  mutate(n = case_when(is.na(n) ~ as.numeric(0),
                       TRUE ~ as.numeric(n)))

write_rds(geocoded, here("data", "geocoded_raw.rds"))

complete_geocode = bind_rows(geocoded %>%
                               select(region, country, lon, lat) %>%
                               distinct(),
                             missing_coords) %>%
  mutate(region = str_to_sentence(region),
         country = str_to_sentence(country)) %>%
  drop_na(lon, lat) %>%
  group_by(region, country) %>%
  arrange(lon, lat) %>%
  slice(1)

# Human sequence locations
geocoded %>%
  filter(str_detect(host, "Human|Homo sapiens")) %>%
  select(lon, lat) %>%
  left_join()

map_regional <- geocoded %>%
  drop_na(region) %>%
  filter(!str_detect(region, "unknown|UNKNOWN")) %>%
  left_join(., missing_coords, by = "region") %>%
  mutate(lon = case_when(region %in% missing_coords$region ~ lon.y,
                         TRUE ~ lon.x),
         lat = case_when(region %in% missing_coords$region ~ lat.y,
                         TRUE ~ lat.x)) %>%
  drop_na(lon, lat) %>%
  group_by(lon, lat) %>%
  summarise(n = sqrt(n())) %>%
  st_as_sf(coords = c("lon", "lat")) %>%
  st_set_crs(value = st_crs(map_seqs))

map_sequences <- ggplot() +
  geom_sf(data = map_seqs, aes(fill = cut(n, c(1, 10, 40, 100, 500, 1300))), alpha = 0.8) +
  geom_sf(data = map_regional, aes(size = n), alpha = 0.4) +
  scale_fill_viridis_d(labels = c("1-10", "11-40", "41-100", "101-500", "501-1,300", "None")) +
  theme_minimal() +
  annotation_north_arrow(height = unit(1, "cm"),
                         style = north_arrow_minimal(text_size = 8)) +
  annotation_scale(height = unit(0.1, "cm"),
                   location = "tr") +
  guides(size = "none") +
  labs(fill = "Sequences")

# Human/rodent ------------------------------------------------------------
prop.table(table(all_sequences$host, useNA = "ifany"))

host_sequences <- all_sequences %>%
  filter(!host %in% c("Cavia porcellus (guinea pig)", "NULL")) %>%
  mutate(host = case_when(host %in% c("Mastomys", "Mastomys sp") ~ "Mastomys sp.",
                          TRUE ~ host),
         genus = case_when(host == "Homo sapiens" ~ "Human",
                           TRUE ~ "Rodent")) %>%
  group_by(country, host, year, genus) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_col(aes(x = year, y = n, fill = host)) +
  facet_wrap(~ genus, scales = "free_y") +
  theme_minimal() +
  labs(x = "Year",
       y = "Number of sequences",
       fill = "Host") +
  scale_fill_viridis_d(labels = c("Homo sapiens (1419)",
                                  "Hylomyscus pamfi (10)",
                                  "Mastomys erythroleucus (18)",
                                  "Mastomys natalensis (609)",
                                  "Mastomys sp. (55)",
                                  "Mus baoulei (9)"),
                       direction = -1) +
  theme(legend.position = "bottom")

first_row <- plot_grid(annual_sequences, host_sequences, labels = "AUTO")
second_row <- plot_grid(map_sequences, labels = "C")
save_plot(here("outputs", "combined_fig1.png"),
          plot_grid(first_row,
          second_row,
          nrow = 2,
          rel_heights = c(1, 2)),
          base_height = 18,
          base_width = 16)
