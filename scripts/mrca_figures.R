source(here::here("scripts", "library_api.r"))
# Run 3_descriptives.R first

mrca_l <- read_xlsx(here("data", "MRCA_L_segment.xlsx"))
mrca_s <- read_xlsx(here("data", "MRCA_S_segment.xlsx"))

mrca_combined <- bind_rows(mrca_s %>%
                             rename("accession_number" = S_segment_accession,
                                    "tmrca" = MRCA) %>%
                             mutate(segment = "S segment"),
                           mrca_l %>%
                             rename("accession_number" = L_segment_accession) %>%
                             mutate(segment = "L segment",
                                    accession_number = paste0(accession_number, ".1"))) %>%
  mutate(segment = factor(segment, levels = c("S segment", "L segment")))


mrca_enriched <- mrca_combined %>%
  left_join(., all_sequences, by = c("accession_number" = "accession_lassa")) %>%
  mutate(region = coalesce(region.x, region.y),
         region = str_to_sentence(region)) %>%
  left_join(., complete_geocode, by = c("region")) %>%
  drop_na(region) %>%
  select(accession_number, segment, host, region, country, year, tmrca, lon, lat) %>%
  drop_na(lon, lat) %>%
  st_as_sf(coords = c("lon", "lat")) %>%
  st_set_crs(value = "EPSG:4326") %>%
  mutate(region = case_when(str_detect(region, "Ebonyi") ~ "Ebonyi",
                            str_detect(region, "Ed") ~ "Edo",
                            str_detect(region, "Plat") ~ "Plateau",
                            str_detect(region, "Tar") ~ "Taraba",
                            TRUE ~ region))

nigeria_shapefile <- read_rds(here("data_download", "gadm36_NGA_1_sp.rds")) %>%
  st_as_sf() %>%
  st_set_crs(value = "EPSG:4326")
nigeria_bbox <- st_bbox(nigeria_shapefile)

mrca_summary <- tibble(mrca_enriched) %>%
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) %>%
  group_by(segment, x, y, region, tmrca) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(tmrca) %>%
  mutate(tmrca_discrete = cut(tmrca, c(1400, 1600, 1700, 1800, 1900, 1950, 1980, 1990, 2000, 2010, 2020), dig.lab = 4, include.lowest = TRUE))

S_segment_mrca <- ggplot() +
  geom_sf(data = nigeria_shapefile, fill = "white") +
  coord_sf(xlim = c(nigeria_bbox[1], nigeria_bbox[3]), ylim = c(nigeria_bbox[2], nigeria_bbox[4])) +
  geom_point(data = mrca_summary %>%
               filter(segment == "S segment"),
             aes(x = x,
                 y = y,
                 colour = tmrca_discrete),
             alpha = 1,
             size = 0.6,
             position = position_jitter(width = 0.25, height = 0.25)) +
  scale_colour_viridis_d(guide = guide_coloursteps(show.limits = TRUE), drop = FALSE) +
  theme_bw() +
  labs(colour = "MRCA",
       x = element_blank(),
       y = element_blank()) +
  annotation_north_arrow(height = unit(1, "cm"),
                         style = north_arrow_minimal(text_size = 8)) +
  annotation_scale(height = unit(0.1, "cm"),
                   location = "tr")

L_segment_mrca <- ggplot() +
  geom_sf(data = nigeria_shapefile, fill = "white") +
  coord_sf(xlim = c(nigeria_bbox[1], nigeria_bbox[3]), ylim = c(nigeria_bbox[2], nigeria_bbox[4])) +
  geom_point(data = mrca_summary %>%
               filter(segment == "L segment"),
             aes(x = x,
                 y = y,
                 colour = tmrca_discrete),
             alpha = 1,
             size = 0.6,
             position = position_jitter(width = 0.25, height = 0.25)) +
  scale_colour_viridis_d(guide = guide_coloursteps(show.limits = TRUE), drop = FALSE) +
  theme_bw() +
  labs(colour = "MRCA",
       x = element_blank(),
       y = element_blank()) +
  annotation_north_arrow(height = unit(1, "cm"),
                         style = north_arrow_minimal(text_size = 8)) +
  annotation_scale(height = unit(0.1, "cm"),
                   location = "tr")

legend <- get_legend(L_segment_mrca +
  theme(legend.direction = "horizontal",
        legend.key.width = unit(2, "cm")))

figure_3 <- plot_grid(plotlist = list(S_segment_mrca +
                                        theme(legend.position = "none"),
                                      L_segment_mrca +
                                        theme(legend.position = "none"),
                                      legend),
                      ncol = 1, 
                      labels = c("A", "B", " "),
                      rel_heights = c(1, 1, 0.2))

save_plot(figure_3, filename = here("outputs", "Figure_3.png"), base_width = 5, base_height = 8)
save_plot(figure_3, filename = here("outputs", "Figure_3.pdf"), base_width = 5, base_height = 8)
