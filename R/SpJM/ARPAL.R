library(ARPALData)
dat=get_ARPA_Lombardia_AQ_municipal_data(
  Date_begin = "2019-01-01",
  Date_end = "2019-12-31"
)

dat_w=get_ARPA_Lombardia_W_data(
  ID_station = NULL,
  Date_begin = "2021-01-01",
  Date_end = "2021-12-31",
  Frequency = "monthly",
  Var_vec = c("Rainfall","Temperature","Wind_speed",
              "Wind_direction"),
  Fns_vec=c("sum","mean","mean","mean"),
  parallel=T
)
