library(tidyverse)
library(lubridate)
last_day <- function(date) {
  ceiling_date(date, "month", change_on_boundary = TRUE) - days(1)
}

log_diff <- function(x) {
  100*(log(x) - lag(log(x), 1))
}
scale_diff <- function(x) {
  x <- (x-10)/10
  x <- x - dplyr::lag(x)
}

data(full_tbl)

full_tbl2 <- unnest(full_tbl)   %>%
  filter(fcst_date > "2000-01-01" & date > "1996-07-01" & fcst_date < "2016-01-01")

monthly_tbl <- full_tbl2 %>%
  group_by(fcst_date) %>%
  mutate(infl = log_diff(infl), ip = log_diff(ip), eti = scale_diff(eti)) %>%
  nest()

quarterly_tbl <- unnest(monthly_tbl) %>%
  mutate(quarter = quarter(.$date, with_year = TRUE)) %>%
  group_by(quarter, fcst_date) %>%
  summarise(unemp = mean(unemp), infl = mean(infl), ip = mean(ip), eti = mean(eti), gdp = mean(gdp, na.rm = TRUE),
            interest = mean(interest), n_months = n()) %>%
  filter(n_months == 3) %>%
  ungroup() %>%
  select(-n_months) %>%
  mutate(date = paste0(floor(quarter), "-", round((quarter-floor(quarter))*30), "-01") %>% ymd %>% last_day)  %>%
  select(date, everything()) %>%
  replace_na(replace = list(gdp = NA)) %>%
  select(-quarter) %>%
  group_by(fcst_date) %>%
  nest()

data_list <- as.list(inner_join(monthly_tbl, quarterly_tbl, by = "fcst_date") %>%
  rename(mf = data.x, qf = data.y))

data_list$mf <- lapply(data_list$mf, function(x) {
  x <- as.data.frame(x)
  rownames(x) <- x$date
  x <- x[, -1]
})
data_list$qf <- lapply(data_list$qf, function(x) {
  x <- as.data.frame(x)
  rownames(x) <- x$date
  x <- x[, -1]
})
