
library(bioinactivation2)
library(tidyverse)
library(readxl)
library(ggsci)
library(ggpubr)

## Load data

all_d <- excel_sheets("./data/datos isotermicos Janellys y Paola.xlsx") %>%
  set_names(., .) |>
  map(
    ~ read_excel("./data/datos isotermicos Janellys y Paola.xlsx",
                 sheet = .)
  ) |>
  map(
    ~ set_names(., c("time", "temp", "logN"))
  ) |>
  map(~ filter(., temp < 61))

## Plot

all_d |>
  imap_dfr(~ mutate(.x, cond = .y)) |>
  separate(cond, into = c("media", "who")) |>
  ggplot(aes(x = time, y = logN, colour = media)) +
  geom_point() +
  geom_line() +
  facet_grid(who~temp,
             scales = "free")

## Two-steps fitting

Bigelow_2step <- all_d |>
  map(
    ~ bioinactivation2:::fit_two_step(
      .,
      "Bigelow",
      known = c(D_temp_xref = 55),
      secondary_models = list(
        list(
          par = "D",
          model = "Bigelow",
          depends_on = c("temp")
        )
      )
    )
  )

Mafart_2step <- all_d |>
  map(
    ~ bioinactivation2:::fit_two_step(
      .,
      "Mafart",
      known = c(delta_temp_xref = 55),
      secondary_models = list(
        list(
          par = "delta",
          model = "Bigelow",
          depends_on = c("temp")
        )
      )
    )
  )

##

list(
  `Bigelow` = Bigelow_2step |> map(~ .$primary_pars) |> imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart` = Mafart_2step |> map(~ .$primary_pars) |> imap_dfr(~ mutate(.x, cond = .y))
) |>
  imap_dfr(~ mutate(.x, model = .y)) |>
  mutate(D = if_else(is.na(D), delta, D)) |>
  separate(cond, into = c("medium", "who")) |>
  filter(who == "Janellys") |>
  ggplot(aes(x = temp, y = D, colour = model)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # facet_grid("model") +
  scale_y_log10()

##

Bigelow_1step <- all_d |>
  map(
    ~ mutate(.,
             logN0 = mean(na.rm = TRUE, if_else(time == 0, logN, NA)),
             .by = c("temp")
    ) |>
      mutate(log_diff = logN - logN0)
  ) |>
  map(
    ~ bioinactivation::fit_isothermal_inactivation("Bigelow",
                                                   .,
                                                   list(D_R = 1, z = 5),
                                                   list(temp_ref = 55)
    )
  )

Mafart_1step <- all_d |>
  map(
    ~ mutate(.,
             logN0 = mean(na.rm = TRUE, if_else(time == 0, logN, NA)),
             .by = c("temp")
    ) |>
      mutate(log_diff = logN - logN0)
  ) |>
  map(
    ~ bioinactivation::fit_isothermal_inactivation("Mafart",
                                                   .,
                                                   list(delta_ref = 1, z = 5, p = 1),
                                                   list(temp_ref = 55)
    )
  )

## All estimates

list(
  `Bigelow_1step` = Bigelow_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_1step` = Mafart_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Bigelow_2step` = Bigelow_2step |>
    map(~ summary(.)$D$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_2step` = Mafart_2step |>
    map(~ summary(.)$delta$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y))
) |>
  imap_dfr(~ mutate(.x, model = .y)) |>
  separate(model, into = c("model", "approach"), sep = "_") |>
  separate(cond, into = c("media", "who"), sep = " ") |>
  mutate(par = if_else(par == "ref", "D_R", par),
         par = if_else(par == "delta_ref", "D_R", par),
         par = if_else(par == "temp_z", "z", par),
         ) %>%
  filter(par != "p") %>%
  split(.$who) |>
  map(
    ~ ggplot(., aes(x = media, colour = approach, y = Estimate)) +
      geom_point(aes(shape = model),
                 position = position_dodge(width = .2)
                 ) +
      geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                        ymax = Estimate + `Std. Error`,
                        linetype = model
                        ),
                    position = position_dodge(width = .2),
                    width = .2) +
      facet_wrap("par", scales = "free")
  ) %>%
  ggarrange(plotlist = ., labels = c("Janellys", "Paola"))


##

Bigelow_1step |> map(summary)
Mafart_1step |> map(summary)
Mafart_1step |> map(plot, make_gg = TRUE)
Bigelow_1step |> map(plot, make_gg = TRUE)

## Fitted curves

p <- list(Bigelow_1step$`TSB Paola`,
          Bigelow_1step$`granad Paola`) |>
  map(plot, make_gg = TRUE) |>
  map(
    ~ . + theme_pubr() +
      scale_color_bmj() +
      labs(x = "Time (min)", y = "Log-reductions (log CFU/mL)")
  ) %>%
  ggarrange(plotlist = .,
            ncol = 1,
            labels = "AUTO")

ggsave("iso_curves_paola.png", width = 8, height = 6)

p <- list(Bigelow_1step$`TSB Janellys`,
          Bigelow_1step$`coco Janellys`) |>
  map(plot, make_gg = TRUE) |>
  map(
    ~ . + theme_pubr() +
      scale_color_bmj() +
      labs(x = "Time (min)", y = "Log-reductions (log CFU/mL)")
  ) %>%
  ggarrange(plotlist = .,
            ncol = 1,
            labels = "AUTO")

p

ggsave("iso_curves_janellis.png", width = 8, height = 6)


## Plot of parameters for Paola

p <- list(
  `Bigelow_1step` = Bigelow_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_1step` = Mafart_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Bigelow_2step` = Bigelow_2step |>
    map(~ summary(.)$D$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_2step` = Mafart_2step |>
    map(~ summary(.)$delta$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y))
) |>
  imap_dfr(~ mutate(.x, model = .y)) |>
  separate(model, into = c("model", "approach"), sep = "_") |>
  separate(cond, into = c("media", "who"), sep = " ") |>
  mutate(par = if_else(par == "ref", "D_R", par),
         par = if_else(par == "delta_ref", "D_R", par),
         par = if_else(par == "temp_z", "z", par),
  ) %>%
  filter(par != "p") %>%
  filter(who == "Paola") |>
  mutate(media = if_else(media == "granad", "Granadilla", "TSB")) |>
  mutate(par = if_else(par == "D_R", "D-value at 55ºC (min)", "z-value (ºC)")) |>
  ggplot(aes(x = media, colour = approach, y = Estimate)) +
  geom_point(aes(shape = model),
             position = position_dodge(width = .2)
  ) +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`,
                    linetype = model
  ),
  position = position_dodge(width = .2),
  width = .2) +
  facet_wrap("par", scales = "free") +
  theme_pubr() +
  scale_color_bmj() +
  labs(x = "", y = "", linetype = "Model", shape = "Model", colour = "Approach")

p

ggsave(p, filename = "iso_pars_paola.png", width = 6, height = 4)

## Plot of parameters for janellis

p <- list(
  `Bigelow_1step` = Bigelow_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_1step` = Mafart_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Bigelow_2step` = Bigelow_2step |>
    map(~ summary(.)$D$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_2step` = Mafart_2step |>
    map(~ summary(.)$delta$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y))
) |>
  imap_dfr(~ mutate(.x, model = .y)) |>
  separate(model, into = c("model", "approach"), sep = "_") |>
  separate(cond, into = c("media", "who"), sep = " ") |>
  mutate(par = if_else(par == "ref", "D_R", par),
         par = if_else(par == "delta_ref", "D_R", par),
         par = if_else(par == "temp_z", "z", par),
  ) %>%
  filter(par != "p") %>%
  filter(who == "Janellys") |>
  mutate(media = if_else(media == "coco", "Coconut water", "TSB")) |>
  mutate(par = if_else(par == "D_R", "D-value at 55ºC (min)", "z-value (ºC)")) |>
  ggplot(aes(x = media, colour = approach, y = Estimate)) +
  geom_point(aes(shape = model),
             position = position_dodge(width = .2)
  ) +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                    ymax = Estimate + `Std. Error`,
                    linetype = model
  ),
  position = position_dodge(width = .2),
  width = .2) +
  facet_wrap("par", scales = "free") +
  theme_pubr() +
  scale_color_bmj() +
  labs(x = "", y = "", linetype = "Model", shape = "Model", colour = "Approach")

p

ggsave(p, filename = "iso_pars_janellis.png", width = 6, height = 4)

## Export parameter table

list(
  `Bigelow_1step` = Bigelow_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_1step` = Mafart_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Bigelow_2step` = Bigelow_2step |>
    map(~ summary(.)$D$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_2step` = Mafart_2step |>
    map(~ summary(.)$delta$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y))
) |>
  imap_dfr(~ mutate(.x, model = .y)) |>
  separate(model, into = c("model", "approach"), sep = "_") |>
  separate(cond, into = c("media", "who"), sep = " ") |>
  write_excel_csv2("pars_isothermal.csv")

list(
  `Bigelow_1step` = Bigelow_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_1step` = Mafart_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Bigelow_2step` = Bigelow_2step |>
    map(~ summary(.)$D$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_2step` = Mafart_2step |>
    map(~ summary(.)$delta$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y))
) |>
  imap_dfr(~ mutate(.x, model = .y)) |>
  separate(model, into = c("model", "approach"), sep = "_") |>
  separate(cond, into = c("media", "who"), sep = " ") |>
  filter(who == "Paola") |>
  mutate(e = paste(round(Estimate, 2), " +- ", round(`Std. Error`, 2))) |>
  select(media, par, e, model, approach) |>
  mutate(par = if_else(par == "temp_z", "z", par)) |>
  mutate(par = if_else(par == "delta_ref", "D_R", par)) |>
  mutate(par = if_else(par == "ref", "D_R", par)) |>
  pivot_wider(names_from = "par", values_from = "e") |>
  arrange(media, model, approach) |>
  write_excel_csv2("pars_paola.csv")

list(
  `Bigelow_1step` = Bigelow_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_1step` = Mafart_1step |>
    map(~ summary(.)$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Bigelow_2step` = Bigelow_2step |>
    map(~ summary(.)$D$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y)),
  `Mafart_2step` = Mafart_2step |>
    map(~ summary(.)$delta$par) |>
    map(~ as_tibble(., rownames = "par")) |>
    imap_dfr(~ mutate(.x, cond = .y))
) |>
  imap_dfr(~ mutate(.x, model = .y)) |>
  separate(model, into = c("model", "approach"), sep = "_") |>
  separate(cond, into = c("media", "who"), sep = " ") |>
  filter(who != "Paola") |>
  mutate(e = paste(round(Estimate, 2), " +- ", round(`Std. Error`, 2))) |>
  select(media, par, e, model, approach) |>
  mutate(par = if_else(par == "temp_z", "z", par)) |>
  mutate(par = if_else(par == "delta_ref", "D_R", par)) |>
  mutate(par = if_else(par == "ref", "D_R", par)) |>
  pivot_wider(names_from = "par", values_from = "e") |>
  arrange(media, model, approach) |>
  write_excel_csv2("pars_janellis.csv")









