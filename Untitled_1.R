
library(tidyverse)
library(stars)
library(terra)
library(furrr)

plan(multisession)


# *****************************************************************************

download.file("https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/KAGRFI/GO4YTX",
              "/mnt/pers_disk/ag_risk_data/cassava.tif",
              method = "wget",
              quiet = T)

download.file("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_bio.zip",
              "/mnt/pers_disk/ag_risk_data/bio_hist.zip",
              "wget",
              quiet = T)

dir.create("/mnt/pers_disk/ag_risk_data/bio_hist")
unzip("/mnt/pers_disk/ag_risk_data/bio_hist.zip",
      exdir = "/mnt/pers_disk/ag_risk_data/bio_hist/")

future_walk2(
  c(
    "https://geodata.ucdavis.edu/cmip6/2.5m/MPI-ESM1-2-HR/ssp585/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2021-2040.tif",
    "https://geodata.ucdavis.edu/cmip6/2.5m/MPI-ESM1-2-HR/ssp585/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2041-2060.tif",
    "https://geodata.ucdavis.edu/cmip6/2.5m/MPI-ESM1-2-HR/ssp585/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2061-2080.tif",
    "https://geodata.ucdavis.edu/cmip6/2.5m/MPI-ESM1-2-HR/ssp585/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp585_2081-2100.tif"
  ),
  
  c(
    "/mnt/pers_disk/ag_risk_data/bio_2021-2040.tif",
    "/mnt/pers_disk/ag_risk_data/bio_2041-2060.tif",
    "/mnt/pers_disk/ag_risk_data/bio_2061-2080.tif",
    "/mnt/pers_disk/ag_risk_data/bio_2081-2100.tif"
  ),
  
  function(f, n){
    
    download.file(f, n, "wget", quiet = T)
    
  }
  
)



# *****************************************************************************

bio_hist <- 
  tibble(f =
           "/mnt/pers_disk/ag_risk_data/bio_hist/" %>% 
           list.files(full.names = T),
         
         v = f %>% str_sub(start = 50) %>% parse_number()) %>% 
  arrange(v) %>% 
  pull(f) %>% 
  read_stars(RasterIO = list(nXOff = 3890, nYOff = 1720, nXSize = 220, nYSize = 180))

nn <- 
  names(bio_hist) %>% 
  str_remove("wc2.1_2.5m_") %>% 
  str_remove(".tif")

bio_hist <- 
  bio_hist %>% 
  setNames(nn)

cassava <- 
  "/mnt/pers_disk/ag_risk_data/cassava.tif" %>% 
  read_stars() %>% 
  st_warp(bio_hist) %>% 
  setNames("area")

tb <- 
  c(cassava,
  bio_hist) %>% 
  as_tibble()


bio_future <- 
  "/mnt/pers_disk/ag_risk_data/" %>% 
  list.files(full.names = T) %>% 
  str_subset("20") %>% 
  read_stars(proxy = F,
             RasterIO = list(nXOff = 3890, nYOff = 1720, nXSize = 220, nYSize = 180)) %>% 
  setNames(c("t-2021-2040", "t-2041-2060", "t-2061-2080", "t-2081-2100")) %>% 
  merge(name = "time") %>% 
  split("band") %>% 
  setNames(nn)





# *****************************************************************************

library(tidymodels)

model_spec <- 
  rand_forest(trees = 1000) %>% 
  set_engine("ranger", num.threads = 8) %>% 
  set_mode("regression")

model_rec <- 
  recipe(area ~ ., data = tb) %>% 
  update_role(x, y, new_role = "ID") %>% 
  step_naomit(everything())
  # prep() %>% 
  # bake(new_data = NULL)

model_wf <- 
  workflow() %>% 
  add_model(model_spec) %>% 
  add_recipe(model_rec)

model_fit <- 
  fit(model_wf, data = tb)


tb_hist <- tb %>% as_tibble() %>% na.omit()

tb_pred_hist <- 
  tb_hist %>% 
  select(x, y, area) %>% 
  bind_cols(predict(model_fit, tb_hist))


tb_pred_hist %>% 
  pivot_longer(c(area, .pred)) %>% 
  mutate(name = factor(name, levels = c("area", ".pred")),
         value = raster::clamp(value, 0, 0.4)) %>% 
  
  ggplot(aes(x,y, fill = value)) +
  geom_raster() +
  coord_equal() +
  facet_wrap(~name, nrow = 1) +
  colorspace::scale_fill_continuous_sequential("plasma",
                                               rev = F,
                                               na.value = "transparent",
                                               name = "ha*1000",
                                               breaks = c(0, 0.08, 0.16, 0.24, 0.32, 0.4),
                                               labels = c("", 
                                                          "0.08 (unsuitable)",
                                                          "0.16 (marginal)",
                                                          "0.24 (suitable)",
                                                          "0.32 (optimal)",
                                                          "")) +
  theme(axis.title = element_blank())
  


# *****************************************************************************

l_tb <- 
  bio_future %>% 
  as_tibble() %>% 
  group_split(time)

l_fut <- 
  c("2021-2040", "2041-2060", "2061-2080", "2081-2100") %>% 
  map2_dfr(l_tb, function(t, d){
    
    tb <- 
      d %>% 
      select(-time) %>% 
      na.omit()
    
    tb %>% 
      select(x,y) %>% 
      bind_cols(predict(model_fit, tb)) %>% 
      mutate(time = t)
    
  })


l_fut %>% 
  mutate(.pred = raster::clamp(.pred, 0, 0.4)) %>% 
  
  ggplot(aes(x,y, fill = .pred)) +
  geom_raster() +
  coord_equal() +
  facet_wrap(~time, nrow = 1) +
  colorspace::scale_fill_continuous_sequential("plasma",
                                               rev = F,
                                               na.value = "transparent",
                                               name = "ha*1000",
                                               breaks = c(0, 0.08, 0.16, 0.24, 0.32, 0.4),
                                               labels = c("", 
                                                          "0.08 (unsuitable)",
                                                          "0.16 (marginal)",
                                                          "0.24 (suitable)",
                                                          "0.32 (optimal)",
                                                          "")) +
  theme(axis.title = element_blank())

