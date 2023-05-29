source(here::here("scripts", "library_api.r"))

search_date <- "2021-09-24"

geocoded <- read_rds(here("data", "geocoded_raw.rds"))

iso_3 <- getData("ISO3") %>%
  filter(NAME %in% unique(geocoded$country)) %>%
  pull(ISO3)

gadm_1 <- lapply(iso_3, function(x) getData(name = "GADM", country = x, level = 1, path = here("data_download")))

gadm_2 <- lapply(iso_3, function(x) getData(name = "GADM", country = x, level = 2, path = here("data_download")))

gadm_1 <- do.call(bind, gadm_1) %>%
  st_as_sf() %>%
  mutate(NAME_1 = str_replace_all(NAME_1, "é", "e"),
         NAME_1 = str_replace_all(NAME_1, "ô", "o"),
         NAME_1 = case_when(str_detect(NAME_1, "GrandBassa") ~ "grand bassa",
                            str_detect(NAME_1, "GrandGedeh") ~ "grand gedeh",
                            str_detect(NAME_1, "GrandKru") ~ "grand kru",
                            TRUE ~ str_to_lower(NAME_1)))

gadm_2 <- do.call(bind, gadm_2) %>%
  st_as_sf() %>%
  mutate(NAME_1 = str_replace_all(NAME_1, "é", "e"),
         NAME_1 = str_replace_all(NAME_1, "ô", "o"),
         NAME_1 = str_replace_all(NAME_1, "è", "e"),
         NAME_1 = str_to_lower(NAME_1),
         NAME_2 = str_replace_all(NAME_2, "é", "e"),
         NAME_2 = str_replace_all(NAME_2, "ô", "o"),
         NAME_2 = str_replace_all(NAME_2, "è", "e"),
         NAME_2 = str_to_lower(NAME_2))

# Associate cases with level 1 and level 2 administrative regions
epi <- read_csv(here("data", "confirmed_cases.csv")) %>%
  mutate(confirmed_cases = replace_na(confirmed_cases, 0),
         region = str_replace_all(region, "é", "e"),) %>%
  group_by(country, region) %>%
  summarise(regional_cases = sum(confirmed_cases)) %>%
  mutate(regional_cases = replace_na(regional_cases, 0),
         country = str_to_title(country),
         region = case_when(region == "fct" ~ "federal capital territory",
                            TRUE ~ str_to_lower(str_replace_all(region, "_", " ")))) %>%
  ungroup() %>%
  mutate(ID = row_number())

country_level_cases <- read_csv(here("data", "confirmed_cases.csv")) %>%
  filter((year >= 2012 & country == "nigeria" & region == "all")|
           (year < 2012 & country == "nigeria") |
           country != "nigeria") %>%
  group_by(country) %>%
  summarise(confirmed_cases  = sum(confirmed_cases, na.rm = TRUE))

epi_level_1 <- c("Benin", "Burkina Faso", "Ghana", "Liberia", "Mali", "Nigeria", "Togo")
epi_level_2 <- c("Guinea", "Sierra Leone")

cases_level_2 <- gadm_2 %>%
  select(GID_0, NAME_0, NAME_1, NAME_2) %>%
  filter(NAME_0 %in% epi_level_2) %>%
  left_join(epi,
            by = c("NAME_0" = "country",
                   "NAME_2" = "region")) %>%
  mutate(regional_cases = replace_na(regional_cases, 0)) %>%
  group_by(GID_0, NAME_1) %>%
  summarise(regional_cases = sum(regional_cases))

cases_level_1 <- gadm_1 %>%
  select(GID_0, NAME_0, NAME_1) %>%
  filter(NAME_0 %in% epi_level_1) %>%
  left_join(epi %>%
              select(country, region, regional_cases), by = c("NAME_0" = "country",
                                                              "NAME_1" = "region")) %>%
  mutate(regional_cases = replace_na(regional_cases, 0)) %>%
  group_by(GID_0, NAME_1) %>%
  summarise(regional_cases = sum(regional_cases))

cases_combined <- bind_rows(cases_level_1, cases_level_2) %>%
  tibble() %>%
  select(GID_0, NAME_1, regional_cases)

# All samples
clean_geocode <- geocoded %>%
  drop_na(region, lon, lat) %>%
  mutate(region = case_when(str_detect(region, "Abuja|FCT") ~ "federal capital territory",
                            str_detect(region, "Adamawa State") ~ "adamawa",
                            str_detect(region, "EBONYI|Ebonyi|Ebonyi State") ~ "ebonyi",
                            str_detect(region, "EDDO|EDO|Edo|Edo State") ~ "edo",
                            str_detect(region, "Enugu State") ~ "enugu",
                            str_detect(region, "Kissedougou") ~ "kissidougou",
                            str_detect(region, "N'Zerekore|Nzerekore") ~ "nzerekore",
                            str_detect(region, "Nasarawa|Nasarawa|Nasarawa State|NASARAWA") ~ "nassarawa",
                            str_detect(region, "Ondo state") ~ "ondo",
                            str_detect(region, "Plateau|PLATEAU|Plateau State") ~ "plateau",
                            str_detect(region, "Taraba|TARABA|Tararba") ~ "taraba",
                            str_detect(region, "unknown|UNKNOWN") ~ as.character(NA),
                            TRUE ~ str_to_lower(region)),
         level_1 = case_when(region %in% gadm_1$NAME_1 ~ TRUE,
                             TRUE ~ FALSE),
         level_2 = case_when(region %in% gadm_2$NAME_2 ~ TRUE,
                             TRUE ~ FALSE),
         exact = case_when(level_1 == FALSE & level_2 == FALSE ~ TRUE,
                           TRUE ~ FALSE))

sequences_level_1 <- clean_geocode %>%
  ungroup() %>%
  st_as_sf(coords = c("lon", "lat"), crs = crs(gadm_1)) %>%
  st_join(gadm_1, st_within) %>%
  tibble() %>%
  group_by(GID_0, NAME_1) %>%
  summarise(n_sequences = n())

human_level_1_sequences <- clean_geocode %>%
  ungroup() %>%
  filter(str_detect(host, "Human|Homo sapiens|Homo Sapiens")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = crs(gadm_1)) %>%
  st_join(gadm_1, st_within) %>%
  tibble() %>%
  group_by(GID_0, NAME_1) %>%
  summarise(n_sequences = n())

rodent_level_1_sequences <- clean_geocode %>%
  ungroup() %>%
  filter(str_detect(host, "Mast|Mus|Rodent")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = crs(gadm_1)) %>%
  st_join(gadm_1, st_within) %>%
  tibble() %>%
  group_by(GID_0, NAME_1) %>%
  summarise(n_sequences = n())

# aggregate the cases and number of sequences at level 1 administrative regions
# use the centroid of these for x and y coordinates and produce a smooth using cases and sequences
# check for population density effect too

level_1_sf <- gadm_1 %>%
  left_join(sequences_level_1, by = c("GID_0", "NAME_1")) %>%
  left_join(cases_combined, by = c("GID_0", "NAME_1")) %>%
  select(GID_0, NAME_1, n_sequences, regional_cases) %>%
  mutate(n_sequences = replace_na(n_sequences, 0),
         regional_cases = replace_na(regional_cases, 0))

if(!file.exists(here("data", "level_1_pop_count.rds"))) {
  # Rasters of population count have been downloaded from https://hub.worldpop.org/geodata/listing?id=78 and saved in data_download
  pop_rasters <- lapply(list.files(here("data_download"), pattern = "_2020.tif", full.names = TRUE), function(x) rast(x)) %>% sprc()
  single_pop_raster <- merge(pop_rasters)
  t <- aggregate(single_pop_raster, fact = 20, fun = "sum", na.rm = TRUE)
  pop_level_1 <- terra::extract(t, gadm_1, fun = sum, na.rm = TRUE)
  
  write_rds(pop_level_1, here("data", "level_1_pop_count.rds"))
} else {
  
  pop_level_1 <- read_rds(here("data", "level_1_pop_count.rds"))
  
}

df_1 <- level_1_sf %>%
  cbind(pop_level_1) %>%
  rename(population_count = ben_2020) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(lon = st_coordinates(st_centroid(geometry))[1],
         lat = st_coordinates(st_centroid(geometry))[2]) %>%
  tibble() %>%
  mutate(cases_per_100000 = (regional_cases/population_count) * 100000) %>%
  select(GID_0, NAME_1, n_sequences, regional_cases, population_count, cases_per_100000, lon, lat)



# Describe sequencing bias ------------------------------------------------

case_sequence_bias <- geocoded %>%
  filter(!is.na(country)) %>%
  group_by(country) %>%
  summarise(n_sequences = n()) %>%
  left_join(country_level_cases %>%
              mutate(country = case_when(country == "cote d'ivoire" ~ "Côte d'Ivoire",
                                         TRUE ~ str_to_title(country))))

ggplot(data = case_sequence_bias, aes(x = confirmed_cases, y = n_sequences)) +
  geom_point(aes(colour = country))

cor.test(case_sequence_bias$n_sequences, case_sequence_bias$confirmed_cases)

case_human_sequence_bias <- geocoded %>%
  filter(!is.na(country) & str_detect(host, "Human|Homo")) %>%
  group_by(country) %>%
  summarise(n_sequences = n()) %>%
  left_join(country_level_cases %>%
              mutate(country = case_when(country == "cote d'ivoire" ~ "Côte d'Ivoire",
                                         TRUE ~ str_to_title(country))))

cor.test(case_human_sequence_bias$n_sequences, case_human_sequence_bias$confirmed_cases)

case_rodent_sequence_bias <- geocoded %>%
  filter(!is.na(country) & !str_detect(host, "Human|Homo")) %>%
  group_by(country) %>%
  summarise(n_sequences = n()) %>%
  left_join(country_level_cases %>%
              mutate(country = case_when(country == "cote d'ivoire" ~ "Côte d'Ivoire",
                                         TRUE ~ str_to_title(country))))

cor.test(case_rodent_sequence_bias$n_sequences, case_rodent_sequence_bias$confirmed_cases)

# Model sequencing bias ---------------------------------------------------

m_1 <- gam(n_sequences ~ s(lon, lat, bs = "ts", k = 102),
    family = "tw",
    data = df_1,
    select = TRUE)

gam.check(m_1)

rect_border <- tibble(x = c(-17.64, -17.64, 15.995642, 15.995642),
                      y = c(4.2, 27.683098, 4.2, 27.683098)) %>%
  st_as_sf(coords = c("x", "y"), crs = crs(gadm_1)) %>%
  summarise() %>%
  st_cast("POLYGON")

w_africa <- read_rds(here("data_download", "level_0_admin.rds")) %>%
  filter(!GID_0 %in% c("CMR", "DZA", "ESH", "MAR", "TCD", "CPV"))

country_border <- w_africa %>%
  summarise()

inverse_country <- st_difference(rect_border, country_border)

# Reduce the opacity of the grid lines: Default is 255
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

# Diverging palette
brbg_hcl <- colorspace::diverging_hcl(11,
                                      h = c(180, 50), c = 80, l = c(20, 95), power = c(0.7, 1.3),
                                      register = "diverging_map")

viz_m1 <- getViz(m_1)

model_1_raster <- plot(sm(viz_m1, 1), n = 300, too.far = 0.08) +
  l_fitRaster(pTrans = zto1(0.05, 2, 0.1)) +
  l_fitContour() +
  coord_cartesian(expand = FALSE)
 
sequencing_effort <- model_1_raster +
  geom_sf(data = inverse_country, fill = "white", colour = "white", inherit.aes = FALSE) +
  geom_sf(data = w_africa, fill = NA, alpha = 1, lwd = 0.5, inherit.aes = FALSE) +
  scale_fill_continuous_diverging(palette = "diverging_map", na.value = "#ffffff00", limits = c(-12, 12)) +
  theme_minimal() +
  labs(title = element_blank(),
       x = element_blank(),
       y = element_blank(),
       fill = "Relative \nsequencing effort") +
  annotation_north_arrow(height = unit(2, "cm"),
                         style = north_arrow_minimal(text_size = 10)) +
  annotation_scale(height = unit(0.1, "cm"),
                   location = "tr") +
  theme(panel.ontop = TRUE,
        panel.grid = element_line(color = col_grid)) +
  theme_bw()

# save_plot(filename = here("outputs", "Figure_2.png"), plot = as_grob(sequencing_effort$ggObj), base_height = 7, base_width = 9)
# save_plot(filename = here("outputs", "Figure_2.pdf"), plot = as_grob(sequencing_effort$ggObj), base_height = 7, base_width = 9)

m_2 <- gam(n_sequences ~ s(lon, lat, bs = "ts", k = 102) + s(cases_per_100000),
           family = "tw",
           data = df_1,
           select = TRUE)

gam.check(m_2)
summary(m_2)

viz_m2 <- getViz(m_2)

model_2_raster <- plot(sm(viz_m2, 1), n = 300, too.far = 0.08) +
  l_fitRaster(pTrans = zto1(0.05, 2, 0.1)) +
  l_fitContour() +
  coord_cartesian(expand = FALSE)

sequencing_effort_2 <- model_2_raster +
  geom_sf(data = inverse_country, fill = "white", colour = "white", inherit.aes = FALSE) +
  geom_sf(data = w_africa, fill = NA, alpha = 1, lwd = 0.5, inherit.aes = FALSE) +
  scale_fill_continuous_diverging(palette = "diverging_map", na.value = "#ffffff00", limits = c(-12, 12)) +
  theme_minimal() +
  labs(title = element_blank(),
       x = element_blank(),
       y = element_blank(),
       fill = "Relative \nsequencing effort") +
  annotation_north_arrow(height = unit(2, "cm"),
                         style = north_arrow_minimal(text_size = 10)) +
  annotation_scale(height = unit(0.1, "cm"),
                   location = "tr") +
  theme(panel.ontop = TRUE,
        panel.grid = element_line(color = col_grid)) +
  theme_bw()

save_plot(filename = here("outputs", "Figure_2.png"), plot = as_grob(sequencing_effort_2$ggObj), base_height = 7, base_width = 9)
save_plot(filename = here("outputs", "Figure_2.pdf"), plot = as_grob(sequencing_effort_2$ggObj), base_height = 7, base_width = 9)

# By host
host_geo <- geocoded %>%
  drop_na(lon, lat, region) %>%
  mutate(host_clean = case_when(is.null(host) ~ as.character(NA),
                                host == "NULL" ~ as.character(NA),
                          str_detect(host, "Homo|Human") ~ "Human",
                          TRUE ~ "Rodent")) %>%
  filter(host != "NULL") %>%
  group_by(lon, lat, host_clean) %>%
  summarise(samples = n()) %>%
  ungroup() %>%
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(gadm_1))

level_1_host <- st_join(host_geo, gadm_1, join = st_within) %>%
  tibble() %>%
  select(samples, host_clean, NAME_0, GID_1) %>%
  group_by(NAME_0, GID_1, host_clean) %>%
  summarise(all_samples = sum(samples, na.rm = TRUE)) %>%
  right_join(gadm_1, by = c("NAME_0", "GID_1")) %>%
  st_as_sf(crs = st_crs(gadm_1))

level_1_host <- level_1_host %>%
  tibble() %>%
  select(NAME_0, GID_1, all_samples, host_clean) %>%
  drop_na(host_clean) %>%
  group_split(host_clean)

level_1_host <- lapply(level_1_host, function(x) 
  x %>%
    full_join(gadm_1 %>%
                filter(NAME_0 %in% x$NAME_0)))

plot_level_1 <- lapply(level_1_host, function(x)
  st_as_sf(x, crs = st_crs(gadm_1)) %>%
    mutate(log_all_samples = log10(all_samples)) %>%
    ggplot() +
    geom_sf(aes(fill = log_all_samples)) +
    geom_sf(data = w_africa, fill = NA, lwd = 0.6, colour = "black") +
    labs(fill = bquote('Number of sequences '(log[10])),
         title = paste(unique(x$host_clean))) +
    scale_fill_continuous(limits = c(0, 3)) +
    theme_bw() +
    theme(legend.direction = "horizontal",
          legend.position = "bottom"))

combined_plot <- plot_grid(plotlist = list(plot_level_1[[1]] +
                                             theme(legend.position = "none"),
                                           plot_level_1[[2]] +
                                             theme(legend.position = "none")),
                           ncol = 1)
legend_b <- get_legend(plot_level_1[[1]] + 
    guides(color = guide_legend(nrow = 1)))

plot_legend <- plot_grid(combined_plot, legend_b, ncol = 1, rel_heights = c(1, .1))

save_plot(filename = here("outputs", "Figure_1.png"), plot = plot_legend, base_height = 8, base_width = 6)
save_plot(filename = here("outputs", "Figure_1.pdf"), plot = plot_legend, base_height = 8, base_width = 6)

# Supplementary figure for country names
country_locations <- w_africa %>%
  ggplot() +
  geom_sf() +
  geom_sf_label_repel(aes(label = NAME_0)) +
  theme_bw() +
  labs(x = element_blank(),
       y = element_blank())

save_plot(filename = here("outputs", "country_locations.png"), plot = country_locations, base_height = 8, base_width = 6)
