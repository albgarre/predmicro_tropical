
library(bioinactivation2)
library(tidyverse)
library(readxl)
library(ggsci)
library(ggpubr)

## Load dynamic data

d <- excel_sheets("./data/Resultados dinamicos Janellys y Paola UTP revisados.xlsx") %>%
  set_names(., .) |>
  map(
    ~ read_excel("./data/Resultados dinamicos Janellys y Paola UTP revisados.xlsx",
                 sheet = .)
  ) |>
  map(
    ~ select(., time = Time, temp = temperature, logN = LogN)
  ) |>
  map(
    ~ mutate(., temp = as.numeric(gsub(",", ".", temp)))
  )

## Predictions based on primary models

bioinactivation::predict_inactivation(
  "Bigelow",
  seq(0, 12.5, length = 100),
  c(N0 = 1e5, Bigelow_1step$`TSB Janellys`$parameters),
  d$`Janellys 2C TSB` |> select(time, temperature = temp)
) |>
  plot(plot_temp  = TRUE) +
  geom_point(aes(x = time, y = logN),
             data = d$`Janellys 2C TSB`)

bioinactivation::predict_inactivation(
  "Bigelow",
  seq(0, 2.5, length = 100),
  c(N0 = 1e5, Bigelow_1step$`TSB Janellys`$parameters),
  d$`Janellys 10C TSB` |> select(time, temperature = temp)
) |>
  plot(plot_temp  = TRUE) +
  geom_point(aes(x = time, y = logN),
             data = d$`Janellys 10C TSB`)

bioinactivation::predict_inactivation(
  "Bigelow",
  seq(0, 12.5, length = 100),
  c(N0 = 1e5, Bigelow_1step$`coco Janellys`$parameters),
  d$`Janellys 2C Coco` |> select(time, temperature = temp)
) |>
  plot(plot_temp  = TRUE) +
  geom_point(aes(x = time, y = logN),
             data = d$`Janellys 2C Coco`)

bioinactivation::predict_inactivation(
  "Bigelow",
  seq(0, 2.5, length = 100),
  c(N0 = 1e5, Bigelow_1step$`coco Janellys`$parameters),
  d$`Janellys 10C coco` |> select(time, temperature = temp)
) |>
  plot(plot_temp  = TRUE) +
  geom_point(aes(x = time, y = logN),
             data = d$`Janellys 10C coco`)

## Fitting dynamic models

secondary_models <- list(
  list(
    par = "D",
    model = "Bigelow",
    depends_on = c("temp")
  )
)

dynamic_guess <- c(
  logN0 = 5,
  D_temp_logz = log10(7),
  # D_temp_z = 7,
  D_logref = 1
)

dynamic_known <- c(
  D_temp_xref = 55
)

check_inactivation_guess(
  "dynamic",
  d[[1]],
  "Bigelow",
  c(dynamic_guess, dynamic_known),
  secondary_models,
  d[[1]]
)

all_fits <- d[1:4] |>
  map(
    ~ fit_inactivation(
      approach = "dynamic",
      fit_data = select(., time, logN),
      primary_model_name = "Bigelow",
      guess = dynamic_guess,
      known = dynamic_known,
      secondary_models = secondary_models,
      env_conditions = select(., time, temp),
      algorithm = "regression"
    )
  )

## Plots of the fits

p1 <- bioinactivation::predict_inactivation(
  "Bigelow",
  seq(0, 12.5, length = 100),
  c(N0 = 1e5, Bigelow_1step$`TSB Janellys`$parameters),
  d$`Janellys 2C TSB` |> select(time, temperature = temp)
) |>
  plot(plot_temp  = TRUE,
       label_y1 = "Microbial concentration (log CFU/g)",
       label_y2 = "Temperature (ºC)",
       ylims = c(0, 6)) +
  geom_point(aes(x = time, y = logN),
             data = d$`Janellys 2C TSB`) +
  geom_line(aes(x = time, y = logN),
            data = all_fits$`Janellys 2C TSB`$best_prediction$simulation,
            linetype = 3,
            linewidth = 1,
            colour = ggsci::pal_bmj()(2)[2]
  ) +
  theme_pubr(base_size = 14)

p2 <- bioinactivation::predict_inactivation(
  "Bigelow",
  seq(0, 2.5, length = 100),
  c(N0 = 1e5, Bigelow_1step$`TSB Janellys`$parameters),
  d$`Janellys 10C TSB` |> select(time, temperature = temp)
) |>
  plot(plot_temp  = TRUE,
       label_y1 = "Microbial concentration (log CFU/g)",
       label_y2 = "Temperature (ºC)",
       ylims = c(0, 6)) +
  geom_point(aes(x = time, y = logN),
             data = d$`Janellys 10C TSB`) +
  geom_line(aes(x = time, y = logN),
            data = all_fits$`Janellys 10C TSB`$best_prediction$simulation,
            linetype = 3,
            linewidth = 1,
            colour = ggsci::pal_bmj()(2)[2]
  ) +
  theme_pubr(base_size = 14)

p3 <- bioinactivation::predict_inactivation(
  "Bigelow",
  seq(0, 12.5, length = 100),
  c(N0 = 1e5, Bigelow_1step$`coco Janellys`$parameters),
  d$`Janellys 2C Coco` |> select(time, temperature = temp)
) |>
  plot(plot_temp  = TRUE,
       label_y1 = "Microbial concentration (log CFU/g)",
       label_y2 = "Temperature (ºC)",
       ylims = c(0, 6)) +
  geom_point(aes(x = time, y = logN),
             data = d$`Janellys 2C Coco`) +
  geom_line(aes(x = time, y = logN),
            data = all_fits$`Janellys 2C Coco`$best_prediction$simulation,
            linetype = 3,
            linewidth = 1,
            colour = ggsci::pal_bmj()(2)[2]
  ) +
  theme_pubr(base_size = 14)

p4 <- bioinactivation::predict_inactivation(
  "Bigelow",
  seq(0, 2.5, length = 100),
  c(N0 = 1e5, Bigelow_1step$`coco Janellys`$parameters),
  d$`Janellys 10C coco` |> select(time, temperature = temp)
) |>
  plot(plot_temp  = TRUE,
       label_y1 = "Microbial concentration (log CFU/g)",
       label_y2 = "Temperature (ºC)",
       ylims = c(0, 6)) +
  geom_point(aes(x = time, y = logN),
             data = d$`Janellys 10C coco`) +
  geom_line(aes(x = time, y = logN),
            data = all_fits$`Janellys 10C coco`$best_prediction$simulation,
            linetype = 3,
            linewidth = 1,
            colour = ggsci::pal_bmj()(2)[2]
  ) +
  theme_pubr(base_size = 14)

p <- ggarrange(p1, p2, p3, p4,
               labels = "AUTO")

p

ggsave(p, filename = "dynamic_curves_janellis.png", width = 12, height = 8)

## Model parameters

p <- all_fits |>
  map(~ summary(.)$par) |>
  map(~ as_tibble(., rownames = "par")) |>
  imap_dfr(~ mutate(.x, cond = .y)) |>
  separate(cond, into = c("foo", "profile", "media")) |>
  filter(par != "logN0") |>
  mutate(par = if_else(par == "D_logref", "D-value at 55ºC (log min)",
                       "z-value (log ºC)")) |>
  mutate(media = if_else(media == "TSB", "TSB", "Coconut water")) |>
  mutate(profile = if_else(profile == "10C", "10ºC/min", "2ºC/min")) |>
  ggplot(aes(x = media, y = Estimate, colour = profile)) +
  geom_point(position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`),
                width = .2,
                position = position_dodge(.3)) +
  facet_wrap("par", scales = "free") +
  labs(x = "", y = "", colour = "Heating rate") +
  theme_pubr() +
  scale_color_bmj()

p

ggsave(p, filename = "dynamic_pars_janellis.png", width = 8, height = 4)

## Export dynamic parameters

all_fits |>
  map(~ summary(.)$par) |>
  map(~ as_tibble(., rownames = "par")) |>
  imap_dfr(~ mutate(.x, cond = .y)) |>
  separate(cond, into = c("foo", "profile", "media")) |>
  # filter(par != "logN0") |>
  # mutate(par = if_else(par == "D_logref", "D-value at 55ºC (log min)",
  #                      "z-value (log ºC)")) |>
  # mutate(media = if_else(media == "TSB", "TSB", "Coconut water")) |>
  # mutate(profile = if_else(profile == "10C", "10ºC/min", "2ºC/min")) |>
  write_excel_csv2(file = "dynamic_pars_janellys.csv")



all_fits |>
  map(~ summary(.)$par) |>
  map(~ as_tibble(., rownames = "par")) |>
  imap_dfr(~ mutate(.x, cond = .y)) |>
  separate(cond, into = c("foo", "profile", "media")) |>
  filter(par == "D_logref") |>
  mutate(aa = 10^Estimate)





























