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
# CREATED: 03/26/2025
# LATEST:  06/08/2025
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


### Locate data
data_file <- "fake_data.csv"
dict_file <- "codebook.csv"

out_dir <- ""
out_suffix <- as.character(Sys.Date())

#Automatically set working directory, requires rstudioapi package. If not working set WD manually
#to do so manually go to Session -> Set Working Directory -> To Source File location
#copy paste result from console over following commented line of code:
#Example:
# setwd("~/Desktop/network")
setwd(dirname(getActiveDocumentContext()$path))

# Get funcs
source("utils.R")
source("metrics.R")
source("viz.R")

# Flags
MAKE_FIGS <- FALSE

############# ANALYSIS #############
#Importing data
df_input <- read.csv(
    data_file,
    stringsAsFactors = FALSE,
    colClasses = c(zip = "character")
) #this reads in zip retains 0s

df_dict <- read.csv(dict_file)

############### EGO ###############

################# Extracting values #################
### record_id
record_id = df_input$record_id

### age
age = df_input$age

### sex
sex_map <- get_codebook_mapping(df_dict, "sex")
sex <- factor(df_input$sex, levels = sex_map, labels = names(sex_map))

### race
race_map <- get_codebook_mapping(df_dict, "race")
ego_races <- extract_ego_multi_attributes(df_input, "race", race_map)

### Employment
employ_map <- get_codebook_mapping(df_dict, "employment")
employment <- factor(
    df_input$employment,
    levels = employ_map,
    labels = names(employ_map)
)

### Education
edu_map <- get_codebook_mapping(df_dict, "edu")
education <- factor(df_input$edu, levels = edu_map, labels = names(edu_map))

### Zip
zip = as.character(df_input$zip)

### Occupation
occ_map <- get_codebook_mapping(df_dict, "occupation")
occupation <- factor(
    df_input$occupation,
    levels = occ_map,
    labels = names(occ_map)
)

### Income
income_map <- get_codebook_mapping(df_dict, "income")
income <- factor(
    df_input$income,
    levels = income_map,
    labels = names(income_map)
)

### Married
marry_map <- get_codebook_mapping(df_dict, "married")
married <- factor(
    df_input$married,
    levels = marry_map,
    labels = names(marry_map)
)

## Live_alone
live_map <- get_codebook_mapping(df_dict, "live_alone")
live_alone <- factor(
    df_input$live_alone,
    levels = live_map,
    labels = names(live_map)
)

### household_number
household_number <- df_input$household_number

### Smoking
smoke_map <- get_codebook_mapping(df_dict, "smoke")
smoke <- factor(df_input$smoke, levels = smoke_map, labels = names(smoke_map))

### Alcohol
alc_map <- get_codebook_mapping(df_dict, "alcohol")
alcohol <- factor(df_input$alcohol, levels = alc_map, labels = names(alc_map))

### Exercise
exercise_map <- get_codebook_mapping(df_dict, "exercise")
exercise <- factor(
    df_input$exercise,
    levels = exercise_map,
    labels = names(exercise_map)
)

### Diet
diet_map <- get_codebook_mapping(df_dict, "diet")
diet <- factor(df_input$diet, levels = diet_map, labels = names(diet_map))

### Health probs
health_map <- get_codebook_mapping(df_dict, "health")
ego_health_probs <- extract_ego_multi_attributes(df_input, "health", health_map)

################# Network Stats #################
network_size = calc_total_alters_df(df_input)
density = calc_egoless_density_df(df_input)

gra_list <- organize_list_tidygraphs(df_input)
constraint = sapply(gra_list, calc_node_constraint)
effsize = sapply(gra_list, calc_node_ens)
mean_degree = sapply(gra_list, calc_egoless_mean_degree)
max_degree = sapply(gra_list, calc_egoless_max_degree)

############### ALTERS ###############
age_sd <- calc_numeric_attr_sd(df_input, "age")

################# IQV #################
gender_map <- get_codebook_mapping(df_dict, "name1sex")
iqv_sex <- calc_attribute_iqv(df_input, "sex", gender_map)

educ_map <- list(
    only_high_school = c(1, 2),
    some_college = c(3, 4),
    college_grad = c(5, 6),
    dont_know = 99
)
iqv_educ <- calc_attribute_iqv(df_input, "educ", educ_map)

race_map <- get_codebook_mapping(df_dict, "name1race")
iqv_race <- calc_attribute_iqv(df_input, "race", race_map)

################# Proportions #################
prop_q1_in_network <- prop_of_qnames_in_network(df_input, 1)
prop_q2_in_network <- prop_of_qnames_in_network(df_input, 2)
prop_q3_in_network <- prop_of_qnames_in_network(df_input, 3)

# All mappings will be the same, so just take the first one
relat_map <- get_codebook_mapping(df_dict, "name1relat")
kin_prop <- calc_prop_alters_multians(
    df_input,
    c("Spouse", "Family"),
    relat_map,
    "relat"
)

speak_map <- get_codebook_mapping(df_dict, "name1speak")
weak_freq_prop <- calc_prop_alters_singleans(
    df_input,
    c("Monthly", "Less often"),
    speak_map,
    "speak"
)

len_map <- get_codebook_mapping(df_dict, "name1length")
weak_dur_prop <- calc_prop_alters_singleans(
    df_input,
    c("Less than three", "Three to six"),
    len_map,
    "length"
)

dist_map <- get_codebook_mapping(df_dict, "name1dist")
far_dist_prop <- calc_prop_alters_singleans(
    df_input,
    c("16-50 miles", "50+ miles"),
    dist_map,
    "dist"
)

alc_map <- get_codebook_mapping(df_dict, "name1alcohol")
heavy_drinkers_prop <- calc_prop_alters_singleans(
    df_input,
    c("Yes", "No"),
    alc_map,
    "alcohol"
)

smoke_map <- get_codebook_mapping(df_dict, "name1smoke")
smoking_prop <- calc_prop_alters_singleans(
    df_input,
    c("Yes", "No"),
    smoke_map,
    "smoke"
)

exer_map <- get_codebook_mapping(df_dict, "name1exer")
no_exercise_prop <- calc_prop_alters_singleans(
    df_input,
    "No",
    exer_map,
    "exer"
)

diet_map <- get_codebook_mapping(df_dict, "name1diet")
bad_diet_prop <- calc_prop_alters_singleans(
    df_input,
    "No",
    diet_map,
    "diet"
)

health_map <- get_codebook_mapping(df_dict, "name1health")
health_prob_prop <- calc_prop_alters_multians(
    df_input,
    names(health_map[health_map != 0 & health_map != 99]),
    health_map,
    "health"
)


################# Blau #################
blau_gender <- calc_blau_alter_heterophily(df_input, "sex", gender_map)
blau_educ <- calc_blau_alter_heterophily(df_input, "educ", educ_map)
blau_dist <- calc_blau_alter_heterophily(df_input, "dist", dist_map)
blau_speak <- calc_blau_alter_heterophily(df_input, "speak", speak_map)

length_map <- get_codebook_mapping(df_dict, "name1length")
blau_length <- calc_blau_alter_heterophily(df_input, "length", length_map)

############# MAKE OUTPUTS #############
############### Data Table ###############
my_df_clean = tibble::tibble(
    record_id = record_id,
    age = age,
    sex = sex,
    race1 = sapply(ego_races, "[", 1),
    race2 = sapply(ego_races, "[", 2),
    zip = zip,
    education = education,
    employment = employment,
    occupation = occupation,
    income = income,
    married = married,
    live_alone = live_alone,
    household_number = household_number,
    ego_alcohol = alcohol,
    ego_smoke = smoke,
    ego_exercise = exercise,
    ego_healthy_diet = diet,
    health_problems1 = sapply(ego_health_probs, "[", 1),
    health_problems2 = sapply(ego_health_probs, "[", 2),
    health_problems3 = sapply(ego_health_probs, "[", 3),
    health_problems4 = sapply(ego_health_probs, "[", 4),
    network_size = network_size,
    density = density,
    constraint = constraint,
    effsize = effsize,
    max_degree = max_degree,
    mean_degree = mean_degree,
    kin_prop = round(kin_prop, digits = 2),
    age_sd = round(age_sd, digits = 2),
    IQV_sex = round(iqv_sex, digits = 2),
    IQV_educ = round(iqv_educ, digits = 2),
    weak_freq_prop = round(weak_freq_prop, digits = 2),
    weak_dur_prop = round(weak_dur_prop, digits = 2),
    far_dist_prop = round(far_dist_prop, digits = 2),
    heavy_drinkers_prop = round(heavy_drinkers_prop, digits = 2),
    smoking_prop = round(smoking_prop, digits = 2),
    no_exercise_prop = round(no_exercise_prop, digits = 2),
    bad_diet_prop = round(bad_diet_prop, digits = 2),
    health_prob_prop = round(health_prob_prop, digits = 2),
    blau_gender = round(blau_gender, digits = 2),
    blau_educ = round(blau_educ, digits = 2),
    blau_distance = round(blau_dist, digits = 2),
    blau_length = round(blau_length, digits = 2),
    blau_speak = round(blau_speak, digits = 2),
    prop_q1 = prop_q1_in_network,
    prop_q2 = prop_q2_in_network,
    prop_q3 = prop_q3_in_network
)

df_clean$zip <- as.character(df_clean$zip) #double check zip as string


### Writing files
write.csv(
    df_clean,
    paste(
        out_dir,
        "Clean_Data",
        ifelse(out_suffix != "", paste("_", out_suffix, sep = ""), ""),
        ".csv",
        sep = ""
    ),
    row.names = FALSE
)

if (MAKE_FIGS) {
    ############### FIGURES ###############
    df_relat_names <- names_to_relat(df_input, relat_map)
    list_tidygras <- organize_list_tidygraphs(df_relat_names)
    list_network_plots_labels <- mapply(
        plot_single_network_node_labels,
        tidygra = list_tidygras,
        ego_name = record_id,
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

    # Creating list of graphs
    list_tidygras <- organize_list_tidygraphs(df_input)

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
        edge_size = 0.25 * est_graph_scale,
        node_size = 1 * est_graph_scale
    )
    # Constructing grid of plots, note this is where we set how many columns the grid will contain,
    #  therefore setting the size and aspect ratio of each plot.
    wrap_list_network_plots <- wrap_plots(
        plotlist = list_network_plots,
        ncol = est_column_count
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
