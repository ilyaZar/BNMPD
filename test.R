
<!--header-includes:
  - \usepackage{bbm}
- \usepackage{my_styles}
bibliography: my_bibliography.bib
-->
  <!--./doc/summary_results/TO-DOS-AND-TESTS/doc/00_intro.Rmd-->


  <!--```{r 01_child_documents, child = "./doc/01_child_documents.Rmd"}
```

```{r 02_bibliography, child = "./doc/02_bibliography.Rmd"}
```
-->



# Aim: create animation showing shifting US boundaries
# depends on 17 MB USAboundariesData package
# link to script file that shows chaning state boundaries
# install.packages("USAboundaries")
library(USAboundaries)
library(tidyverse)
library(tmap)
library(sf)
dates = paste(historydata::us_state_populations$year, "01", "01", sep = "-")
dates_unique = unique(dates)
usb1 = USAboundaries::us_boundaries(map_date = dates[1])
usb1$year = lubridate::year(dates_unique[1])
plot(usb1$geometry)
usbl = map(dates, ~USAboundaries::us_boundaries(map_date = .))
# usb = do.call(rbind, usbl)s
statepop = historydata::us_state_populations %>%
  dplyr::select(-GISJOIN) %>% rename(name = state)
sel = usb1$name %in% statepop$name
summary(sel)
usb1$name[!sel]
usbj = left_join(usb1, statepop)
plot(usbj["population"])
i = 2

dates_unique[dates_unique > "2000-12-31"] = "2000-12-31"

for(i in 2:length(dates_unique)) {
  usbi = USAboundaries::us_boundaries(map_date = dates_unique[i])
  print(st_crs(usbi))
  usbi$year = lubridate::year(dates_unique[i])
  if(dates_unique[i] == "2000-12-31") usbi$year = 2010
  plot(usbi$geometry)
  usbji = left_join(usbi, statepop)
  plot(usbji["population"])
  usbj = rbind(usbj, usbji)
}

summary(usbj)
usa_contig = usbji[!grepl(pattern = "Alaska|Haw", usbji$name), ]
plot(usa_contig["population"])
usa_union = st_union(usa_contig) %>%
  st_transform(2163)
plot(usa_union)
bb_contig = st_bbox(usa_union)

# map_dbl(statepop_sf, ~sum(is.na(.)))  # looks about right
# usbj = st_transform(usbj, 4269)
pal = viridis::viridis(n = 7, direction = -1)
pb = c(0, 1, 2, 5, 10, 20, 30, 40) * 1e6
facet_anim = tm_shape(usbj, bbox = bb_contig, projection = 2163) +
  tm_polygons("population", colorNA = NULL, palette = pal, breaks = pb) +
  tm_facets(free.scales.fill = FALSE, ncol = 1, nrow = 1, along = "year") +
  tm_shape(usa_union) + tm_borders(lwd = 2) +
  tm_layout(legend.position = c("left", "bottom"))
tmap_animation(tm = facet_anim, filename = "09-us_pop.gif", width=800, delay=40)


data(NLD_prov)
m1 <- tm_shape(NLD_prov) +
  tm_polygons("yellow") +
  tm_facets(along = "name")

tmap_animation(m1, filename="Dutch_provinces.gif", width=800, delay=40)

data(World, metro)

m2 <- tm_shape(World, simplify = 0.5) +
  tm_fill() +
  tm_shape(metro) +
  tm_bubbles(size = paste0("pop", seq(1970, 2030, by=10)),
             col = "purple",
             border.col = "black", border.alpha = .5,
             scale = 2) +
  tm_facets(free.scales.symbol.size = FALSE, nrow=1,ncol=1) +
  tm_format("World", scale=.5)

tmap_animation(m2, filename="World population.gif", width=1200, delay=100)
