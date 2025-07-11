# CHELSA BIOCLIM+ Data Download ----
# Purpose: Automate downloading of monthly CHELSA climate variables (e.g., hurs, clt, vpd, pet, etc.)
# Source: https://envicloud.wsl.ch/#/?bucket=https://os.zhdk.cloud.switch.ch/chelsav2/&prefix=GLOBAL/monthly/
# Note: This section assumes you are working with monthly rasters in CHELSA V2.1, organized by variable and year
# Recommendation: Run this only once, then comment it out to avoid unnecessary re-downloads

# Define base URL and file path
# base_url <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/monthly/"
# base_dir <- "/Users/jm200/Documents/chelsa_BioclimPlus/monthly"  # Update to your preferred directory

# List of variables to download
# vars <- c("hurs", "clt", "sfcWind", "vpd", "rsds", "pet", "cmi")
# months <- sprintf("%02d", 1:12)
# years <- 1990:2019  # 30-year period ending in 2019

# for (var in vars) {
#   var_dir <- file.path(base_dir, var)
#   dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
#   
#   for (year in years) {
#     for (month in months) {
#       filename <- sprintf("CHELSA_%s_%s_%s_V.2.1.tif", var, month, year)
#       url <- paste0(base_url, var, "/", filename)
#       destfile <- file.path(var_dir, filename)
#       
#       if (!file.exists(destfile)) {
#         cat("⬇️ Downloading:", filename, "\n")
#         tryCatch({
#           download.file(url, destfile, mode = "wb", quiet = TRUE)
#         }, error = function(e) {
#           message("❌ Failed: ", url)
#         })
#       } else {
#         cat("✅ Already exists:", filename, "\n")
#       }
#     }
#   }
# }
