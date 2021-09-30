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
included_phylogenetics_september <- tibble(accession_lassa = c(included_phylogenetics$L_accession,
                                                               included_phylogenetics$S_accession)) %>%
  filter(!is.na(accession_lassa))


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

geocoded <- all_sequences %>%
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

missing_coords = tibble(region = c("Silimi", "Safrani", "Khoria", "Brissa", "Denguedou", "Songo Hospital, Sebgwema", "Bantou", "Tanganya", "Gbetaya", "Odo-akaba"),
                        lon = c(-10.65, -10.738367, -10.893033, -10.688767, -10.448611, -10.948721, -10.583333, -10.972783, -11.040150, 2.601698),
                        lat = c(9.966667,  10.057817, 9.942267, 10.216833, 8.495556, 8.004604, 10.066667, 10.000400, 9.841017, 8.774286))

map_regional <- geocoded %>%
  drop_na(region) %>%
  filter(region != "unknown") %>%
  left_join(., missing_coords, by = "region") %>%
  mutate(lon = case_when(is.na(lon.y) ~ lon.x,
                         TRUE ~ lon.y),
         lat = case_when(is.na(lat.y) ~ lat.x,
                         TRUE ~ lat.y)) %>%
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
