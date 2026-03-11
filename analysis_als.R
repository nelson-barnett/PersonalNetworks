# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PROJECT: Process Persnet
# PURPOSE: Organize and analyze egonetwork data from the persnet instrument
# DIR:     ""
# INPUTS:  Fake data ("fake_data.csv")
# OUTPUTS: A .csv file of network measures, a grid of sociograms, and individual
# sociograms. Built on the legacy code of persnet analyses.
# Requirements: tidyverse (version 2.0.0), tidygraph (1.3.1), stats (4.1.1) and
# purrr (1.0.2); igraph (2.0.3); patchwork (1.3.0), ggraph, patchwork, cowplot,
# rstudioapi
# AUTHORS: Zachary Wehrwein, Liam McCafferty, Amar Dhand
# UPDATED AND ADAPTED BY NELSON BARNETT
# CREATED: 03/26/2025
# LATEST:  03/03/2026
# PSERIES: NA
# NSERIES: NA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This file will not run if the following occurs:
# 1. igraph, tidyverse, tidygraph go through a major update
# 2. the way in which persnet data is loaded changes (e.g. character vs numeric)
# 3. different variable names, but instructions to adjust are included below
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################# Code Setup ##################################

#Attaching packages
library(tidyverse)
library(igraph)
library(cowplot)
library(tidygraph)
library(purrr)
library(stats)
library(ggraph)
library(patchwork)
library(rstudioapi)

#Clearing Working Directory
rm(list = ls())

### Path info
data_file <- ""
out_dir <- ""
dict_file <- "PersonalNetworkSurvey_RedcapCodebook.csv"

# Settings/flags
MAKE_FIGS <- TRUE
use_relat_names_in_figs <- TRUE
replace_ids_in_figs <- TRUE
repel_labels <- FALSE
is_v2 <- TRUE

out_suffix <- as.character(Sys.Date())
metrics_out_name <- "persnet_metrics"

# Get funcs
source("utils.R")
source("metrics.R")
source("viz.R")

############# ANALYSIS #############
#Importing data
df_input <- read.csv(
    data_file,
    stringsAsFactors = FALSE,
    colClasses = c(zip_code = "character", zip_code_v2 = "character")
) #this reads in zip retains 0s

df_dict <- read.csv(dict_file)

############### EGO ###############

################# Extracting values #################
### record_id is always generated
record_id = df_input$record_id

# Set values based on v1/v2 flag
age_col <- ifelse(is_v2, "age_v2", "age")
race_col <- ifelse(is_v2, "race_v2", "race")
employ_col <- ifelse(is_v2, "employment_v2", "employment")
smoke_col <- ifelse(is_v2, "smoke_v2", "smoke")

# Data only collected during demographic v1
if (!is_v2) {
    ### sex
    sex_map <- get_codebook_mapping(df_dict, "sex")
    sex <- names(sex_map[pull(df_input, "sex")])

    ### Education
    edu_map <- get_codebook_mapping(df_dict, "education")
    education <- factor(
        pull(df_input, "education"),
        levels = edu_map,
        labels = names(edu_map)
    )

    # Ethnicity
    eth_map <- get_codebook_mapping(df_dict, "ethnicity")
    ethnicity <- names(eth_map[pull(df_input, "ethnicity")])

    # Region of onset
    region_map <- get_codebook_mapping(df_dict, "region")
    region <- extract_ego_multi_attributes(df_input, "region", region_map)
}

### age
age = pull(df_input, age_col)

### race
race_map <- get_codebook_mapping(df_dict, "race")
ego_races <- extract_ego_multi_attributes(df_input, race_col, race_map)

### Employment
employ_map <- get_codebook_mapping(df_dict, "employment")
ego_employ <- extract_ego_multi_attributes(df_input, employ_col, employ_map)


### Smoking
smoke_map <- get_codebook_mapping(df_dict, "smoke")
smoke <- factor(
    pull(df_input, smoke_col),
    levels = smoke_map,
    labels = names(smoke_map)
)

### Zip
zip = as.character(df_input$zip_code)

################# Network Stats #################
network_size = n_unique_alters(df_input, use_more_names = TRUE)
density = egoless_density_df(df_input)

gra_list <- organize_list_tidygraphs(df_input)
constraint = sapply(gra_list, node_constraint)
effsize = sapply(gra_list, node_ens)
mean_degree = sapply(gra_list, egoless_degree, mean)
max_degree = sapply(gra_list, egoless_degree, max)

############### ALTERS ###############
age_sd <- numeric_attr_sd(df_input, "age")

################# IQV #################
gender_map <- get_codebook_mapping(df_dict, "name1sex")
iqv_sex_datanorm <- attribute_iqv(
    df_input,
    "sex",
    gender_map,
    normalize_by = "data"
)
iqv_sex_mapnorm <- attribute_iqv(
    df_input,
    "sex",
    gender_map,
    normalize_by = "mapping"
)

educ_map <- list(
    only_high_school = c(1, 2),
    some_college = c(3, 4),
    college_grad = c(5, 6),
    dont_know = 99
)
iqv_educ_datanorm <- attribute_iqv(
    df_input,
    "educ",
    educ_map,
    normalize_by = "data"
)
iqv_educ_mapnorm <- attribute_iqv(
    df_input,
    "educ",
    educ_map,
    normalize_by = "mapping"
)

################# Proportions #################
prop_q1_in_network <- prop_of_qnames_in_network(df_input, 1)
prop_q2_in_network <- prop_of_qnames_in_network(df_input, 2)
prop_q3_in_network <- prop_of_qnames_in_network(df_input, 3)
num_unique_alters_namegen <- n_unique_alters(df_input, use_more_names = FALSE)
num_unique_alters_total <- n_unique_alters(df_input, use_more_names = TRUE)

# All mappings will be the same, so just take the first one
relat_map <- get_codebook_mapping(df_dict, "name1relat")
kin_prop <- prop_alters_multians(
    df_input,
    c("Spouse", "Family"),
    relat_map,
    "relat"
)

health_map <- get_codebook_mapping(df_dict, "name1health")
als_illness_prop <- prop_alters_multians(
    df_input,
    c("ALS", "Serious Illness"),
    health_map,
    "health"
)

speak_map <- get_codebook_mapping(df_dict, "name1speak")
weak_freq_prop <- prop_alters_singleans(
    df_input,
    c("Monthly", "Less often"),
    speak_map,
    "speak"
)

len_map <- get_codebook_mapping(df_dict, "name1length")
weak_dur_prop <- prop_alters_singleans(
    df_input,
    c("Less than three", "Three to six"),
    len_map,
    "length"
)

dist_map <- get_codebook_mapping(df_dict, "name1dist")
far_dist_prop <- prop_alters_singleans(
    df_input,
    c("16-50 miles", "50+ miles"),
    dist_map,
    "dist"
)

als_map <- get_codebook_mapping(df_dict, "name1als")
met_through_als_prop <- prop_alters_singleans(
    df_input,
    "Yes",
    als_map,
    "als"
)

################# Blau #################
blau_gender <- blau_alter_heterophily(df_input, "sex", gender_map)
blau_educ <- blau_alter_heterophily(df_input, "educ", educ_map)
blau_dist <- blau_alter_heterophily(df_input, "dist", dist_map)
blau_speak <- blau_alter_heterophily(df_input, "speak", speak_map)

length_map <- get_codebook_mapping(df_dict, "name1length")
blau_length <- blau_alter_heterophily(df_input, "length", length_map)


############# MAKE OUTPUTS #############
############### Data Table ###############

if (is_v2) {
    v1_vars <- c(
        "sex",
        paste0("race", 1:2),
        "education",
        "ethnicity",
        paste0("region", 1:5)
    )

    for (x in v1_vars) {
        assign(x, "see v1")
    }
} else {
    race1 <- sapply(ego_races, "[", 1)
    race2 <- sapply(ego_races, "[", 2)
    region1 <- sapply(region, "[", 1)
    region2 <- sapply(region, "[", 2)
}


df_clean = tibble::tibble(
    record_id = record_id,
    age = age,
    sex = sex,
    race1 = race1,
    race2 = race2,
    region1 = region1,
    region2 = region2,
    ethnicity = ethnicity,
    education = education,
    employment = sapply(ego_employ, "[", 1),
    zip = zip,
    network_size = network_size,
    density = density,
    constraint = constraint,
    effsize = effsize,
    max_degree = max_degree,
    mean_degree = mean_degree,
    kin_prop = kin_prop,
    als_or_illness_prop = als_illness_prop,
    prop_of_network_in_q1 = prop_q1_in_network,
    prop_of_network_in_q2 = prop_q2_in_network,
    prop_of_network_in_q3 = prop_q3_in_network,
    n_unique_alters_namegen = num_unique_alters_namegen,
    n_unique_alters_total = num_unique_alters_total,
    prop_through_als = met_through_als_prop,
    age_sd = round(age_sd, digits = 2),
    iqv_sex_datanorm = round(iqv_sex_datanorm, digits = 2),
    IQV_sex_all_options_norm = round(iqv_sex_mapnorm, digits = 2),
    iqv_educ_datanorm = round(iqv_educ_datanorm, digits = 2),
    iqv_educ_all_options_norm = round(iqv_educ_mapnorm, digits = 2),
    weak_freq_prop = round(weak_freq_prop, digits = 2),
    weak_dur_prop = round(weak_dur_prop, digits = 2),
    far_dist_prop = round(far_dist_prop, digits = 2),
    blau_gender = round(blau_gender, digits = 2),
    blau_educ = round(blau_educ, digits = 2),
    blau_distance = round(blau_dist, digits = 2),
    blau_length = round(blau_length, digits = 2),
    blau_speak = round(blau_speak, digits = 2)
)

df_clean$zip <- as.character(df_clean$zip) #double check zip as string


### Writing files
write.csv(
    df_clean,
    paste(
        out_dir,
        metrics_out_name,
        ifelse(out_suffix != "", paste("_", out_suffix, sep = ""), ""),
        ".csv",
        sep = ""
    ),
    row.names = FALSE
)

############### FIGURES ###############
if (MAKE_FIGS) {
    if (replace_ids_in_figs) {
        df_input <- df_input %>%
            mutate(
                record_id = gsub(".*-", "", df_input$record_id) %>%
                    as.numeric() %>%
                    paste0("P", .)
            )
    }

    df_for_figs <- df_input %>%
        mutate(n_alts_with_data = n_alters_with_data(.)) %>%
        arrange(desc(n_alts_with_data))

    if (use_relat_names_in_figs) {
        names(relat_map) <- gsub("\\s+", "\n", names(relat_map))
        df_for_figs <- names_to_relat(df_for_figs, relat_map)
    }

    list_tidygras <- organize_list_tidygraphs(df_for_figs)

    list_network_plots_labels <- mapply(
        plot_single_network_node_labels,
        tidygra = list_tidygras,
        ego_name = df_for_figs$record_id,
        repel = repel_labels,
        fig_title = paste(
            "ID:",
            df_for_figs$record_id,
            "\nnum alters shown:",
            df_for_figs$n_alts_with_data
        ),
        SIMPLIFY = FALSE
    )

    ############# SAVE OUTPUTS #############

    # This section creates a montage of network graphs in a grid orientation from top
    #  left to bottom right. These graphs are purposefully stripped of alter names
    #  and simplified for a better visualization at smaller sizes. Output is a PDF by default.

    # Default sizing/scaling for this section is set for the fake datasets of n = 6. You may need to
    #  make adjustments for larger/smaller datasets and differently sized outputs. This section
    #  only needs data import and functions in the "Establishing Network Graph Functions" section to work.

    # To customize the montage easily, modify these three variables to create a desired output.
    #  Note that at smaller output sizes and/or larger datasets you may need to reduce the
    #  graph_scale, otherwise the edges/nodes will be too large for their respective graphs.

    ############### Setup ###############

    # Width in inches of the output PDF
    output_width <- 9
    # Height in inches of the output PDF
    output_height <- 5
    # Change the size of the edges/nodes for graphs within montage
    graph_scale <- 4
    ## for ~100 networks, we recommend 8.5, 11, and 2

    # Gives how many graphs the montage will contain
    graph_count <- length(list_tidygras)

    # Area of the output PDF
    workspace_area <- output_width * output_height
    # Estimated area of each graph
    est_graph_size <- sqrt(workspace_area / graph_count)
    # Estimated number of columns to most evenly fill grid
    est_column_count <- ceiling(output_width / est_graph_size)
    # Changing scale of graphs by how large graphs will be
    est_graph_scale <- graph_scale / est_graph_size

    # Constructing graphs, note that this is where our edge/node size scaling actually occurs
    #  and this section determines the ratio of edge size to node size.
    list_network_plots <- lapply(
        list_tidygras,
        plot_single_network,
        edge_size = 0.1 * est_graph_scale,
        node_size = 0.5 * est_graph_scale
    )

    # Constructing grid of plots, note this is where we set how many columns the grid will contain,
    #  therefore setting the size and aspect ratio of each plot.
    wrap_list_network_plots <- wrap_plots(
        plotlist = list_network_plots,
        ncol = est_column_count,
    )

    ############### WRITE ###############
    # Creating output PDF. Note that if you'd like to change the output to something other than a PDF,
    #  change the extention within the output file name and uncomment the last two lines of code.
    ggsave(
        paste(
            out_dir,
            "Social_Network_Grid",
            ifelse(out_suffix != "", paste("_", out_suffix, sep = ""), ""),
            ".pdf",
            sep = ""
        ),
        wrap_list_network_plots,
        width = output_width,
        height = output_height,
        # This code is only needed in case the user wants to change the output to a raster output (not PDF)
        # units = "in",
        # dpi = 300
    )

    i = 0
    single_networks_out_dir <- paste0(
        out_dir,
        "single_networks",
        "_",
        out_suffix,
        "/"
    )
    dir.create(single_networks_out_dir, showWarnings = FALSE)
    for (plt in list_network_plots_labels) {
        i = i + 1
        ggsave(
            paste(
                single_networks_out_dir,
                "Single_Network",
                df_for_figs$record_id[[i]],
                ifelse(out_suffix != "", paste("_", out_suffix, sep = ""), ""),
                ".png",
                sep = ""
            ),
            plt,
            dpi = 150,
            width = 12,
            height = 10
        )
    }

    pdf(
        paste(
            out_dir,
            "Single_Networks",
            ifelse(out_suffix != "", paste("_", out_suffix, sep = ""), ""),
            ".pdf",
            sep = ""
        ),
        width = 12,
        height = 10
    )
    # each print() automatically starts a new page when mfrow = c(1,1)
    par(mfrow = c(1, 1))

    for (plt in list_network_plots_labels) {
        print(plt)
    }

    dev.off()
}
