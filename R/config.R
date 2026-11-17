# ==============================================================================
# config.R -- shared paths and data-access helpers for every script in R/
#
# Every script starts with:   source(here::here("R", "config.R"))
# `here` locates the repository root (marked by the empty file `.here`), so the
# scripts run from any working directory on any machine -- no setwd(), no
# absolute paths.
#
# Folder layout (created automatically):
#   Data/    input data + intermediate .rds files (downloaded on demand)
#   Stan/    .stan model files (versioned with the code)
#   Output/  fitted stanfit objects (.rds); large, not versioned
#   Figure/  figures written by the scripts
#   Tables/  LaTeX / CSV tables written by the scripts
# ==============================================================================
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

PROJ_DIR  <- here::here()
DATA_DIR  <- file.path(PROJ_DIR, "Data")
STAN_DIR  <- file.path(PROJ_DIR, "Stan")
MODEL_DIR <- file.path(PROJ_DIR, "Output")
FIG_DIR   <- file.path(PROJ_DIR, "Figure")
TABLE_DIR <- file.path(PROJ_DIR, "Tables")
# PRISM archive (several GB). Defaults to Data/prism; to keep it elsewhere set the
# environment variable PRISM_DIR (e.g. in ~/.Renviron) instead of editing scripts.
PRISM_DIR <- Sys.getenv("PRISM_DIR", unset = file.path(DATA_DIR, "prism"))
for (d in c(DATA_DIR, MODEL_DIR, FIG_DIR, TABLE_DIR, PRISM_DIR))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ------------------------------------------------------------------------------
# Data archive.
# Source: Zenodo record https://doi.org/10.5281/zenodo.23605794 (record id below).
# Dropbox links are kept only as a fallback. Files are downloaded once into Data/
# (or Output/ for fitted models) and reused thereafter.
# ------------------------------------------------------------------------------
ZENODO_RECORD <- "23605794"

DROPBOX_URLS <- c(
  "Study_site.csv" = "https://www.dropbox.com/scl/fi/si346imz380lpdgo9yekr/Study_site.csv?rlkey=yiue42npkzzu9w8dggjr4fent&dl=1",
  "Initialdata.csv" = "https://www.dropbox.com/scl/fi/b93bvocqltadc36xirak2/Initialdata.csv?rlkey=8hd3z4th35lqvtfvam83kb972&dl=1",
  "census2023.csv" = "https://www.dropbox.com/scl/fi/fkwm0dan6nx2eaeyxjrjw/census2023.csv?rlkey=hy9209t53j9n7vxhta7axl5jk&dl=1",
  "census2024.csv" = "https://www.dropbox.com/scl/fi/52c1hzv97cml698kb74tq/census2024.csv?rlkey=pqiz8g0jgnhxen08j2450w7a8&dl=1",
  "data_ini_2024.csv" = "https://www.dropbox.com/scl/fi/zncghunh1p9ull9j1jhkp/data_ini_2024.csv?rlkey=fkhbp3sdvb65va0gg4rwp4ra4&dl=1",
  "census_2025.csv" = "https://www.dropbox.com/scl/fi/oeqdgik07lyzxbkeiwpfp/census_2025.csv?rlkey=0midqalrvaaqu6i8v4h2z1vpw&dl=1",
  "aghy_occ_raw.rds" = "https://www.dropbox.com/scl/fi/ijl1i7964qxzcxvm7blz8/aghy_occ_raw.rds?rlkey=msbs3xzkjc719uld8yw7hb6cj&dl=1",
  "elvi_occ_raw.rds" = "https://www.dropbox.com/scl/fi/0ssa5gepxyz28b7ykw1x8/elvi_occ_raw.rds?rlkey=4dx0q4lw2112droh73hmh7xte&dl=1",
  "poa_occ_raw.rds" = "https://www.dropbox.com/scl/fi/oip7ndyf0d99rqxcqxb0q/poa_occ_raw.rds?rlkey=920uql1gd4gahnh8utw9fz96l&dl=1",
  "climate_census_years.rds" = "https://www.dropbox.com/scl/fi/fj6aqhej58k9fjt9v02ap/climate_census_years.rds?rlkey=28dsn6b9o3x06wd1c2ehs5ama&dl=1",
  "source_pop.csv" = "https://www.dropbox.com/scl/fi/a6cycehor6k6jp7xvbw3c/source_pop.csv?rlkey=chftk4ym5rsc7edjug90ie6o7&dl=1",
  "herbivory.csv" = "https://www.dropbox.com/scl/fi/2gnlfozxpd2u9gprzp9oi/herbivory.csv?rlkey=sz2cloqxtbc6ou29j97l3t10f&dl=1",
  "fit_surv_abio_bio_endo.rds" = "https://www.dropbox.com/scl/fi/vu750s78j6zjnmbohp8a0/fit_surv_abio_bio_endo.rds?rlkey=52n3g1naikyqzj0ix57ve6e23&dl=1",
  "fit_grow_abio_bio_endo.rds" = "https://www.dropbox.com/scl/fi/wdxva2182fcd2i6jzsjcr/fit_grow_abio_bio_endo.rds?rlkey=v79tfa2br3l0sd6mblas8rto2&dl=1",
  "fit_inf_abio_bio_endo.rds" = "https://www.dropbox.com/scl/fi/857ar8r838n5wauqo3dtz/fit_inf_abio_bio_endo.rds?rlkey=2mot44dp52xlxj19bov5jndly&dl=1",
  "fit_spik_abio_bio_endo.rds" = "https://www.dropbox.com/scl/fi/a6o21wspoofthisda1aye/fit_spik_abio_bio_endo.rds?rlkey=sxsxzzfjjqopqnh5km3s92xpa&dl=1",
  "fit_surv_abio_bio_endo_linear.rds" = "https://www.dropbox.com/scl/fi/mh5es9xqo4t608h12zg4q/fit_surv_abio_bio_endo_linear.rds?rlkey=akzrlhtqbrx3sut9h58aidp0v&dl=1",
  "fit_grow_abio_bio_endo_linear.rds" = "https://www.dropbox.com/scl/fi/mu3gpry42ad7fkfbtlvne/fit_grow_abio_bio_endo_linear.rds?rlkey=pz8uiqevdm5ogy7qvm369qyvb&dl=1",
  "fit_inf_abio_bio_endo_linear.rds" = "https://www.dropbox.com/scl/fi/5lfkgq6d5a2vzx5h2t9yr/fit_inf_abio_bio_endo_linear.rds?rlkey=pfqvm7sg8un5c14slypn1aqa9&dl=1",
  "fit_spik_abio_bio_endo_linear.rds" = "https://www.dropbox.com/scl/fi/7ivmicuigz1pahg4vxa7q/fit_spik_abio_bio_endo_linear.rds?rlkey=8h3js8dnaue95ojom8evrkmkl&dl=1",
  "fit_inf_abio_bio_endo_hurdle_linear.rds" = "https://www.dropbox.com/scl/fi/v5kcirie79uae1dex6l1a/fit_inf_abio_bio_endo_hurdle_linear.rds?rlkey=pxp0pdptn95j5ubkjs819ou3x&dl=1",
  "demography_climate.rds" = "https://www.dropbox.com/scl/fi/b7s8xk3131vpubcqq0413/demography_climate.rds?rlkey=ak5b5dl6t18fhiehv3mgapyfk&dl=1",
  "prism_edge_yr_means.rds" = "https://www.dropbox.com/scl/fi/h8m15t64sxkkghgsnsypp/prism_edge_yr_means.rds?rlkey=ljjf16kqvvoe15qmhrqh8qfah&dl=1",
  "climate_site_scaled.rds" = "https://www.dropbox.com/scl/fi/irecsnoh3xrq6g8d5cysa/climate_site_scaled.rds?rlkey=63r7ugrtkuo5ncmqfoywbidps&dl=1",
  "herbivory.rds" = "https://www.dropbox.com/scl/fi/s5ferd5ha1hjrinvory2f/herbivory.rds?rlkey=vel87d2imoxriu74q3fzuioos&dl=1",
  "reproduction_and_biomass-endo-soil.csv" = "https://www.dropbox.com/scl/fi/tujihmutv87jufadtvbwj/reproduction_and_biomass-endo-soil.csv?rlkey=pa04j3jfae5cldlbc880w3xr3&dl=1",
  "Agrostis_hyemalis_precipitation.csv" = "https://www.dropbox.com/scl/fi/bst61ph31v0ewz6dyn79k/Agrostis_hyemalis_precipitation.csv?rlkey=xp1hq45eu7h8u6wnuznb9zlh0&dl=1",
  "Poa_autumnalis_precipitation.csv" = "https://www.dropbox.com/scl/fi/f6kkq559n795a3pgbdvsi/Poa_autumnalis_precipitation.csv?rlkey=adf5hwslj7uo74ui8c663y9hq&dl=1",
  "Elymus_virginicus_precipitation.csv" = "https://www.dropbox.com/scl/fi/71ospl61twndy42ffz21m/Elymus_virginicus_precipitation.csv?rlkey=qihct3x5kcit7584sp78evcdy&dl=1",
  "fit_surv_sens_noKERSON.rds" = "https://www.dropbox.com/scl/fi/nu2azmsok9mztqqdjjj99/fit_surv_sens_noKERSON.rds?rlkey=xbmxiz0j7woo71fv954jybz2u&dl=1",
  "fit_grow_sens_noKERSON.rds" = "https://www.dropbox.com/scl/fi/e0verxpbqgq0l4ovb5a63/fit_grow_sens_noKERSON.rds?rlkey=0to4b7zkpsog6svnq86kdc4ds&dl=1",
  "fit_inf_sens_noKERSON.rds" = "https://www.dropbox.com/scl/fi/oufjkibr5uk2x1y43ppzh/fit_inf_sens_noKERSON.rds?rlkey=yjrve44p1saofezvgsn29w46x&dl=1",
  "fit_spik_sens_noKERSON.rds" = "https://www.dropbox.com/scl/fi/16lhmj4t0x8ku5ut35nx8/fit_spik_sens_noKERSON.rds?rlkey=51x4r2rm9ya28wr328xy8c2by&dl=1"
)

# Return a local path to `name`, downloading it first if it is not already there.
# Fitted Stan objects (fit_*.rds, herbivory.rds) live in Output/, everything else in Data/.
data_path <- function(name) {
  dir   <- if (grepl("^(fit_|herbivory\\.rds)", name)) MODEL_DIR else DATA_DIR
  local <- file.path(dir, name)
  if (file.exists(local) && file.size(local) > 0) return(local)

  urls <- character(0)
  if (nzchar(ZENODO_RECORD))
    urls <- c(urls, sprintf("https://zenodo.org/records/%s/files/%s?download=1",
                            ZENODO_RECORD, utils::URLencode(name, reserved = TRUE)))
  if (name %in% names(DROPBOX_URLS)) urls <- c(urls, DROPBOX_URLS[[name]])
  if (length(urls) == 0)
    stop("'", name, "' not found in ", dir, " and no download URL is registered.\n",
         "Download the Data/ folder from the Zenodo record ",
         "and place the file in ", dir, ".", call. = FALSE)

  old <- options(timeout = max(3600, getOption("timeout"))); on.exit(options(old))
  for (u in urls) {
    ok <- tryCatch({ utils::download.file(u, local, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (ok && file.exists(local) && file.size(local) > 0) return(local)
    unlink(local)
  }
  stop("Could not download '", name, "'. Tried:\n  ", paste(urls, collapse = "\n  "), call. = FALSE)
}
