source(here::here("scripts", "library_api.r"))

np_gpc <- read_rds(here("cleaned_data", "cleaned_np_gpc_2021-03-15.rds"))
l_z <- read_rds(here("cleaned_data", "cleaned_l_z_2021-03-15.rds"))
s_seq <- read_rds(here("cleaned_data", "cleaned_s_segment_2021-03-15.rds"))
l_seq <- read_rds(here("cleaned_data", "cleaned_l_segment_2021-03-15.rds"))

included_records <- unique(c(np_gpc$accession_lassa, l_z$accession_lassa, s_seq$accession_lassa, l_seq$accession_lassa))

included_records <- tibble(included_records) %>%
  rename("accession_lassa" = included_records) %>%
  left_join(., tibble(bind_rows(np_gpc, l_z, s_seq, l_seq) %>%
                        select(accession_lassa, host, country, region, year) %>%
                        mutate(accession_lassa = as.character(accession_lassa)) %>%
                        distinct(accession_lassa, .keep_all = T)),
            by = "accession_lassa")


# Records mapping ---------------------------------------------------------
geocoded <- included_records %>%
  mutate(region = ifelse(region == "Unknown", NA, region)) %>%
  unite("location", c(region, country), sep = ", ", na.rm = T, remove = F)

coordinate <- geocode(geocoded$location) # This will pull coordinates for the location (built from region and country) for some reason Togo is listed as being in the US

country_coordinate <- unique(geocoded$country) %>%
  bind_cols(geocode(unique(geocoded$country))) %>%
  rename("country" = 1,
         "country_lon" = "lon",
         "country_lat" = "lat") %>%
  mutate(country_lon = ifelse(country == "Togo", 1.217, country_lon),
         country_lat = ifelse(country == "Togo", 6.133, country_lat))

geocoded %<>%
  bind_cols(coordinate) %>%
  mutate(lon = ifelse(location == "Togo", 1.217, lon),
         lat = ifelse(location == "Togo", 6.133, lat))

country_geocode <- geocoded %>%
  left_join(., country_coordinate, by = "country")


country_summary <- country_geocode %>%
  group_by(country) %>%
  summarise(n = n())

west_afr <- c("Benin", "Burkina Faso", "Cote d'Ivoire", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Liberia", "Mali",
              "Mauritania", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo")

data("World")

afr <- World %>%
  filter(name %in% west_afr) %>%
  left_join(., country_summary %>%
              rename("name" = "country"),
            by = "name") %>%
  st_transform(crs = st_crs("+proj=merc"))
  
tm_shape(afr) + 
  tm_fill("n", title = "", style = "log10_pretty", as.count = T,
          textNA = "No records") +
  tm_layout(main.title = "Lassa virus sequences in GenBank",
            legend.outside = T) +
  tm_borders() +
  tm_grid(alpha = 0.2) +
  tm_xlab("Longitude") +
  tm_ylab("Latitude") +
  tm_scale_bar(color.dark = "gray60",
               position = c("right", "top")) + 
  tm_compass(type = "4star", 
             text.size = 0.5, # set size of the compass
             color.dark = "gray60", # color the compass
             text.color = "gray60",
             position = "left") + 
  tm_shape(afr %>% filter(n > 0)) +
  tm_text("name", size = 0.7) +
  tmap_save(filename = here("outputs", "records_map.png"))


# Nigeria mapping ---------------------------------------------------------
nigeria <- read_rds(here("data", "nigeria_sf.rds")) %>%
  st_as_sf()

region_summary <- geocoded %>%
  filter(!is.na(region)) %>%
  group_by(country, region) %>%
  summarise(n = n())

nigeria_region <- region_summary %>%
  filter(country == "Nigeria") %>%
  ungroup() %>%
  mutate(region = recode(region,
                         "Cross river" = "Cross River",
                         "Nasarawa" = "Nassarawa",
                         "Owo" = "Oyo")) %>%
  rename("NAME_1" = "region")

nigeria_region <- nigeria %>%
  full_join(., nigeria_region,
  by = "NAME_1")

tm_shape(nigeria_region) + 
  tm_fill("n", title = "", style = "pretty", as.count = T,
          textNA = "No records") +
  tm_layout(main.title = "Lassa virus sequences from Nigeria",
            legend.outside = T) +
  tm_borders() +
  tm_grid(alpha = 0.2) +
  tm_xlab("Longitude") +
  tm_ylab("Latitude") +
  tm_scale_bar(color.dark = "gray60",
               position = c("right", "bottom")) + 
  tm_compass(type = "4star", 
             text.size = 0.5, # set size of the compass
             color.dark = "gray60", # color the compass
             text.color = "gray60",
             position = "left") + 
  tm_shape(nigeria_region %>% filter(n > 0)) +
  tm_text("NAME_1", size = 0.7) +
  tmap_save(filename = here("outputs", "nigeria_map.png"))


# Human/rodent ------------------------------------------------------------
table(included_records$host)

# Year of collection ----------------------------------------------------------
table(included_records$year)

ggplot(included_records %>%
         group_by(country, year) %>%
         summarise(n = n())) +
  geom_col(aes(x = year, y = n, fill = country)) +
  theme_minimal() +
  labs(x = "Year",
       y = "Number") +
  scale_fill_discrete(name = "Country")
  
included_records %>%
  group_by(year) %>%
  summarise(n = n()) %>%
  filter(year < 2010) %>%
  summarise(n = sum(n))
