source(here::here("scripts", "library_api.r"))

np_gpc <- read_rds(here("cleaned_data", "cleaned_np_gpc_2021-03-15.rds"))
l_z <- read_rds(here("cleaned_data", "cleaned_l_z_2021-03-15.rds"))

geocoded_np_gpc <- np_gpc %>%
  mutate(region = recode(region,
                         "Unknown" = NA)) %>%
  unite("location", c(region, country), sep = ", ", na.rm = T, remove = F) %>%
  mutate(coordinates = geocode(geocoded_np_gpc$location),
         coordinates.lat = if_else(country == "Togo", coordinates.lat = 6.133333, coordinates.lat),
         coordinates.lon = if_else(country == "Togo", coordinates.lon = 1.216667, coordinates.lon))

missing_coords <- tibble(location = c("Togo"),
                         coordinates.lon = c(1.217),
                         coordinates.lat = c(6.133))

         
geocoded_np_gpc %>%
  filter(!is.na(geocoded_np_gpc$coordinates.lon))
