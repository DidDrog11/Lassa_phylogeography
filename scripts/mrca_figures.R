source(here::here("scripts", "library_api.r"))

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
  st_as_sf(coords = c("lon", "lat")) %>%
  st_set_crs(value = "EPSG:4326") %>%
  mutate(region = case_when(str_detect(region, "Ebonyi") ~ "Ebonyi",
                            str_detect(region, "Ed") ~ "Edo",
                            str_detect(region, "Plat") ~ "Plateau",
                            str_detect(region, "Tar") ~ "Taraba",
                            TRUE ~ region))

nigeria_shapefile <- read_rds(here("data", "nigeria_sf.rds")) %>%
  st_set_crs(value = "EPSG:4326")
nigeria_bbox <- st_bbox(nigeria_shapefile)

mrca_summary <- tibble(mrca_enriched) %>%
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) %>%
  group_by(segment, x, y, region, tmrca) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(tmrca) %>%
  pivot_wider(names_from = tmrca,
              values_from = n) %>%
  mutate(across(.cols = 5:62, ~ replace_na(.x, 0)))
mrca_summary$radius <- 0.35

A <- ggplot() +
  geom_sf(data = nigeria_shapefile, alpha = 0.5) +
  coord_sf(xlim = c(nigeria_bbox[1], nigeria_bbox[3]), ylim = c(nigeria_bbox[2], nigeria_bbox[4])) +
  geom_scatterpie(data = mrca_summary, aes(x = x, y = y, group = region, r = radius), cols = colnames(mrca_summary[,c(5:62)]), color = NA) +
  scale_fill_viridis_d() +
  facet_wrap(~ segment) +
  theme_minimal() +
  annotation_north_arrow(height = unit(1, "cm"),
                         style = north_arrow_minimal(text_size = 8)) +
  annotation_scale(height = unit(0.1, "cm"),
                   location = "tr") +
  labs(fill = "MRCA") +
  theme(strip.text = element_text(hjust = 0, size = 14))

B <- ggplot(mrca_enriched) +
  geom_histogram(aes(x = tmrca), fill = "#398053") +
  facet_wrap(~ segment) +
  labs(x = "MRCA",
       y = element_blank(),
       fill = "Region") +
  theme_minimal() +
  theme(strip.text = element_text(hjust = 0))

A_B <- plot_grid(plotlist = list(A, B), align = "v", nrow = 2, rel_heights = c(4, 2), labels = c("A", "B"))

save_plot(A_B, filename = here("outputs", "combined_fig5.png"), base_width = 16, base_height = 18)
