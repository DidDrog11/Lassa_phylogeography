source(here::here("scripts", "library_api.r"))

search_date <- "2021-09-24"

geocoded <- read_rds(here("data", "geocoded_raw.rds"))

iso_3 <- getData("ISO3") %>%
  filter(NAME %in% unique(geocoded$country)) %>%
  pull(ISO3)

gadm <- lapply(iso_3, function(x) getData(name = "GADM", country = x, level = 1, path = here("data_download")))

gadm <- do.call(bind, gadm) %>%
  st_as_sf()

table(!is.na(geocoded$region)) # 1761 sequences have region locations

# All samples

clean_geo <- geocoded %>%
  drop_na(lon, lat, region) %>%
  group_by(lon, lat) %>%
  summarise(samples = n()) %>%
  ungroup() %>%
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(gadm))

level_1 <- st_join(clean_geo, gadm, join = st_within) %>%
  tibble() %>%
  select(samples, NAME_0, GID_1) %>%
  group_by(NAME_0, GID_1) %>%
  summarise(all_samples = sum(samples)) %>%
  right_join(gadm, by = c("NAME_0", "GID_1")) %>%
  mutate(all_samples = replace_na(all_samples, 0)) %>%
  st_as_sf(crs = st_crs(gadm))

epi <- read_xlsx(here("data", "suspected_confirmed_deaths.xlsx")) %>%
  filter(country != "nigeria") %>% 
  group_by(country) %>%
  summarise(national_cases = sum(confirmed_cases)) %>%
  mutate(national_cases = replace_na(national_cases, 0),
         country = str_to_title(country))

nig_epi <- read_xlsx(here("data", "suspected_confirmed_deaths.xlsx")) %>%
  filter(country == "nigeria") %>% 
  group_by(country, region) %>%
  summarise(regional_cases = sum(confirmed_cases)) %>%
  mutate(regional_cases = replace_na(regional_cases, 0),
         country = str_to_title(country),
         region = case_when(region == "fct" ~ "Federal Capital Territory",
                            TRUE ~ str_to_title(str_replace_all(region, "_", " ")))) %>%
  filter(region != "All")

df_1 <- level_1 %>%
  left_join(epi, by = c("NAME_0" = "country")) %>%
  left_join(nig_epi, by = c("NAME_0" = "country", "NAME_1" = "region")) %>%
  mutate(confirmed_cases = coalesce(regional_cases, national_cases)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(confirmed_cases = replace_na(confirmed_cases, 0),
         lon = st_coordinates(st_centroid(geometry))[1],
         lat = st_coordinates(st_centroid(geometry))[2]) %>%
  tibble() %>%
  select(NAME_0, GID_1, all_samples, confirmed_cases, lon, lat)

m_1 <- gam(all_samples ~ s(lon, lat, bs = "ts", k = 102),
    family = "tw",
    data = df_1,
    select = TRUE)

gam.check(m_1)

rect_border <- tibble(x = c(-17.64, -17.64, 15.995642, 15.995642),
                      y = c(4.2, 27.683098, 4.2, 27.683098)) %>%
  st_as_sf(coords = c("x", "y"), crs = crs(gadm)) %>%
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

save_plot(filename = here("outputs", "modeled_effort.png"), plot = as_grob(sequencing_effort$ggObj), base_height = 7, base_width = 9)
save_plot(filename = here("outputs", "modeled_effort.pdf"), plot = as_grob(sequencing_effort$ggObj), base_height = 7, base_width = 9)

m_2 <- gam(all_samples ~ s(lon, lat, bs = "ts", k = 102) + s(confirmed_cases, k = 13),
           family = "tw",
           data = df_1,
           select = TRUE)

gam.check(m_2)

viz_m1 <- getViz(m_2)

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

save_plot(filename = here("outputs", "supplementary_figure_1.png"), plot = as_grob(sequencing_effort_2$ggObj), base_height = 7, base_width = 9)

# By host
host_geo <- geocoded %>%
  drop_na(lon, lat, region) %>%
  mutate(host_clean = case_when(is.null(host) ~ as.character(NA),
                                host == "NULL" ~ as.character(NA),
                          host == "Homo sapiens" ~ "Human",
                          TRUE ~ "Rodent")) %>%
  filter(host != "NULL") %>%
  group_by(lon, lat, host_clean) %>%
  summarise(samples = n()) %>%
  ungroup() %>%
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(gadm))

level_1_host <- st_join(host_geo, gadm, join = st_within) %>%
  tibble() %>%
  select(samples, host_clean, NAME_0, GID_1) %>%
  group_by(NAME_0, GID_1, host_clean) %>%
  summarise(all_samples = sum(samples)) %>%
  right_join(gadm, by = c("NAME_0", "GID_1")) %>%
  st_as_sf(crs = st_crs(gadm))

level_1_host <- level_1_host %>%
  tibble() %>%
  select(NAME_0, GID_1, all_samples, host_clean) %>%
  drop_na(host_clean) %>%
  group_split(host_clean)

level_1_host <- lapply(level_1_host, function(x) 
  x %>%
    full_join(gadm %>%
                filter(NAME_0 %in% x$NAME_0)))

plot_level_1 <- lapply(level_1_host, function(x)
  st_as_sf(x, crs = st_crs(gadm)) %>%
    mutate(log_all_samples = log(all_samples)) %>%
    ggplot() +
    geom_sf(aes(fill = log_all_samples)) +
    geom_sf(data = w_africa, fill = NA, lwd = 0.6, colour = "black") +
    labs(fill = "Number of sequences (log)",
         title = paste(unique(x$host_clean))) +
    scale_fill_continuous(limits = c(0, 7)) +
    theme_bw() +
    theme(legend.direction = "horizontal",
          legend.position = "bottom"))

combined_plot <- plot_grid(plotlist = plot_level_1, ncol = 1)

save_plot(filename = here("outputs", "sequence_locations.png"), plot = combined_plot, base_height = 8, base_width = 6)
save_plot(filename = here("outputs", "sequence_locations.pdf"), plot = combined_plot, base_height = 8, base_width = 6)

# Supplementary figure for country names
country_locations <- w_africa %>%
  ggplot() +
  geom_sf() +
  geom_sf_label_repel(aes(label = NAME_0)) +
  theme_bw()

save_plot(filename = here("outputs", "country_locations.png"), plot = country_locations, base_height = 8, base_width = 6)
