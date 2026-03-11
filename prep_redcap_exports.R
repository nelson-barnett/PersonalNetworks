# Social Connectedness Project collects demographic and persnet surveys at different times
# This combines the individual exports into the expected form for analysis_als.R

### Setup
root_dir <- ""

persnet_path <- file.path(root_dir, "")
demographic_path <- file.path(root_dir, "")
out_name <- "soccon_persnet_demographic"

###### Code
out_meta <- gsub(".*PERSNET(\\d+)_DATA_(.*)_\\d+.csv", "\\1_\\2", persnet_path)
out_dir <- dirname(persnet_path)
out_path <- paste0(file.path(out_dir, paste(out_name, out_meta, sep = "_")), ".csv")

df_persnet <- read.csv(
    persnet_path,
    stringsAsFactors = FALSE,
    colClasses = c(zip_code = "character", zip_code_v2 = "character")
)
df_demographic <- read.csv(
    demographic_path,
    stringsAsFactors = FALSE,
    colClasses = c(zip_code = "character", zip_code_v2 = "character")
)

df_merged <- merge(
    df_persnet,
    df_demographic,
    by = "record_id"
)

write.csv(
    df_merged,
    out_path,
    row.names = FALSE
)

print(paste("File saved to:", out_path))
